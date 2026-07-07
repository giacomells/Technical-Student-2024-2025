import logging
from dataclasses import dataclass

import numpy as np
import xcoll as xc
import xtrack as xt
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
from matplotlib import ticker

from crystal_extraction import plotters

logger = logging.getLogger(__name__)

@dataclass
class BeamArgs:
    emit_rms: float = 10e-6/426
    dpp_fullwidth: float = 2e-3
    nparticles: int = 1000
    p0c: float = 400e9
    mass0: float = xt.PROTON_MASS_EV

@dataclass
class SteinArgs:
    npoints: int = 100
    crystal_angular_acceptance: float = 20e-6   
    cscale: str = 'log'
    septum_loc: str = 'mst.21774'
    crystal_loc: str = 'lsf.43605'
    approach_side: int = - 1 # - 1: beam approaches crystal from "negative x" side 
    nsigma_for_delta0: float = 1.
    nsigma_emit: float = 3
    nstopbands_dpp: int = 3

sps_crystal = xc.EverestCrystal(material=xc.materials.SiliconCrystal,
                                 tilt = 0,
                                 bending_angle=174e-6,  # bending radius [m]
                                 jaw = 0,  # jaw position [m] 
                                 width=2e-3,  # horizontal dimension of crystal [m]
                                 height=30e-3,  # vertical dimension of crystal [m]
                                 length=2.5e-3,  # the crystal length [m]
                                 miscut=0.0,  # just set to 0 if crystal is assumed to be manufactured well
                                 lattice='strip',  # means 110
                                 side='+',  # will only keep positive x side jaw, and will ignore negative x side jaw
                                 )



class Steinbach:
    def __init__(self, line: xt.Line, crystal: xc.EverestCrystal = sps_crystal,
                  beam_args: BeamArgs = BeamArgs(), stein_args: SteinArgs = SteinArgs()):
        self.line = line
        self.crystal = crystal
        self.beam_args = beam_args
        self.stein_args = stein_args

    def parse_xsuite_twiss(self, tw: xt.Table) -> pd.DataFrame:
        tw_pandas = tw.to_pandas()
        tw_pandas.set_index('name', inplace=True)
        tw_pandas['dpp_sign'] = np.sign(self.stein_args.approach_side*tw_pandas.dx)
        tw_pandas['gammx'] = (1 + tw_pandas.alfx**2) / tw_pandas.betx
        tw_pandas['mux_to_sep_deg'] = self._mux_deg_to_ele(tw_pandas, mod=True)
        tw_pandas['mux_to_sep_rad'] = tw_pandas.mux_to_sep_deg * np.pi / 180
        tw_pandas['dpp_stopband_rms'] = tw_pandas.dpp_sign * np.sqrt(tw_pandas.betx) / tw_pandas.dx * self.beam_args.emit_rms**0.5
        tw_pandas['dpp_offset_0'] = tw_pandas.dpp_sign * self.crystal.width / 2 / tw_pandas.dx
        tw_pandas['dpp_offset_max'] = tw_pandas.dpp_offset_0 + tw_pandas.dpp_stopband_rms * self.stein_args.nsigma_emit
        tw_pandas['dxp_hardt_term'] = np.abs(tw_pandas.betx**0.5 * tw_pandas.dpx/tw_pandas.dx + tw_pandas.alfx/tw_pandas.betx**0.5)
        tw_pandas['dxp_stopband'] = self.stein_args.nsigma_emit * np.sqrt(self.beam_args.emit_rms) * tw_pandas.dxp_hardt_term
        tw_pandas['mm_per_mrad_to_sep'] = np.sqrt(tw_pandas.betx * tw_pandas.betx.loc[self.stein_args.septum_loc]) * np.sin(tw_pandas.mux_to_sep_rad)
        tw_pandas['gap_to_sep_mm'] = tw_pandas.mm_per_mrad_to_sep * self.crystal.bending_angle * 1e3
        tw_pandas['delta0'] = tw_pandas.dpx * tw_pandas.dpp_offset_0
        tw_pandas['delta0_aligned'] = tw_pandas.delta0 + tw_pandas.dpp_sign * (
                                     self.stein_args.nsigma_for_delta0*self.beam_args.emit_rms*
                                     (-tw_pandas.alfx / tw_pandas.betx**0.5) 
                                    + self.stein_args.crystal_angular_acceptance
                                    )
        tw_pandas['dpp_stopband_rms_rel'] = tw_pandas.dpp_stopband_rms / self.beam_args.dpp_fullwidth
        tw_pandas['gap_to_sep_in_sigma'] = tw_pandas.gap_to_sep_mm*1e-3 / np.sqrt(tw_pandas.betx * self.beam_args.emit_rms)
        tw_pandas['dxp_stopband_rel'] = tw_pandas.dxp_stopband / self.stein_args.crystal_angular_acceptance

        return tw_pandas

    def _mux_deg_to_ele(self, tw_df, mod=True):
        small_mux = tw_df.loc[self.stein_args.septum_loc].mux - tw_df.mux
        big_mux = tw_df.loc[self.stein_args.septum_loc].mux + (tw_df.mux.max() - tw_df.mux)
        mux_filter = small_mux >= 0
        mux = mux_filter * small_mux + ~mux_filter * big_mux
        if mod:
            return np.mod(mux * 360, 360)
        else:
            return mux * 360

    @staticmethod
    def _points_inside(x, xp, d_crystal, delta_crystal, d0=0, delta0=0):
        return np.sum((x > -d_crystal+d0) & (x-d0 < d_crystal+d0) & 
                      (xp > -delta_crystal + delta0) & 
                      (xp < delta_crystal + delta0)) / x.size * 100

    def _stein_boundaries(self, dpps, tw_crystal, d_crystal, delta_crystal, d0_crystal=0, delta0_crystal=0):
        positive_pos_amp = 1 / tw_crystal['betx']**0.5 * (dpps * tw_crystal['dx'] - d_crystal - d0_crystal)
        negative_pos_amp = -1 / tw_crystal['betx']**0.5 * (dpps * tw_crystal['dx'] + d_crystal - d0_crystal)
        positive_angle_amp = 1 / tw_crystal['gammx']**0.5 * (dpps * tw_crystal['dpx'] - delta_crystal - delta0_crystal)
        negative_angle_amp = -1 / tw_crystal['gammx']**0.5 * (dpps * tw_crystal['dpx'] + delta_crystal - delta0_crystal)

        return positive_pos_amp, negative_pos_amp, positive_angle_amp, negative_angle_amp
    
    
    def _plot_boundaries(self, ax, dpps, tw_crystal, d_crystal, delta_crystal, d0_crystal=0, delta0_crystal=0):
        positive_pos_amp, negative_pos_amp, positive_angle_amp, negative_angle_amp = self._stein_boundaries(
            dpps, tw_crystal, d_crystal, delta_crystal, d0_crystal=d0_crystal, delta0_crystal=delta0_crystal)

        # Plot boundaries
        ax.plot(dpps, positive_pos_amp, color='blue', label='Boundary (pos)')
        ax.plot(dpps, negative_pos_amp, color='blue')
        ax.plot(dpps, positive_angle_amp, color='fuchsia', label='Boundary (angle)')
        ax.plot(dpps, negative_angle_amp, color='fuchsia')

        return ax
    
    def create_xsuite_particles(self, tw=None, plot=False):
        if tw is None:
            tw = self.line.twiss(method='4d')

        tw_pandas = self.parse_xsuite_twiss(tw)
        tw_crystal = tw_pandas.loc[self.stein_args.crystal_loc]

        x_norms = stats.truncnorm.rvs(-self.stein_args.nsigma_emit,
                                      self.stein_args.nsigma_emit,
                                      size=self.beam_args.nparticles)
        px_norms = stats.truncnorm.rvs(-self.stein_args.nsigma_emit,
                                       self.stein_args.nsigma_emit,
                                       size=self.beam_args.nparticles)
        
        amps = np.sqrt(x_norms ** 2 + px_norms ** 2) * np.sqrt(self.beam_args.emit_rms)
        deltas_raw = tw_crystal['dpp_stopband_rms'] * stats.uniform.rvs(0,
                                                                    self.stein_args.nstopbands_dpp,
                                                                    size=self.beam_args.nparticles)
        ddeltas = tw_crystal['dpp_stopband_rms'] * amps / self.beam_args.emit_rms**0.5 + tw_crystal['dpp_offset_0']
        deltas = deltas_raw + ddeltas

        particles = self.line.build_particles(
            method="4d",
            zeta=0,
            delta=deltas,
            x_norm=x_norms,
            px_norm=px_norms,
            y_norm=0,
            py_norm=0,
            nemitt_x=self.beam_args.emit_rms*self.line.particle_ref.gamma0[0]*self.line.particle_ref.beta0[0], 
            nemitt_y=0,
            at_element=self.stein_args.crystal_loc
        )

        self.xsuite_particles = particles

        if plot:
            figs, axes = self.plot(tw)
            axes[0].scatter(deltas, amps, c='red', marker='.')
            axes[1].scatter(particles.x, particles.px, c='red', marker='.')
            for fig in figs:
                fig.tight_layout()

        return self.xsuite_particles
    


    def plot(self, tw=None):
        if tw is None:
            logger.warning('Twiss not provided, computing from line...')
            tw = self.line.twiss(method='4d')

        tw_pandas = self.parse_xsuite_twiss(tw)
        tw_crystal = tw_pandas.loc[self.stein_args.crystal_loc]
        tw_septum = tw_pandas.loc[self.stein_args.septum_loc]

        d_crystal = self.crystal.width / 2
        delta_crystal = self.stein_args.crystal_angular_acceptance / 2
        d0_crystal = 0
        delta0_crystal = tw_crystal['delta0']
        kick_crystal = self.crystal.bending_angle
        mux_rad = tw_crystal['mux_to_sep_rad']
        
        dpp_beam = self.beam_args.dpp_fullwidth
        amp_beam = self.stein_args.nsigma_emit*self.beam_args.emit_rms ** 0.5
        dpp_offset = tw_crystal['dpp_offset_max']
        

        amps, dpps = (np.linspace(0, 2 * amp_beam, self.stein_args.npoints), 
                        np.linspace(-2 * dpp_beam + dpp_offset,
                            2 * dpp_beam + dpp_offset, self.stein_args.npoints))
        

        # Create grid of amplitudes and dpps and compute those inside the crystal
        aa, dd = np.meshgrid(amps, dpps)
        aa_flat, dd_flat = aa.flatten(), dd.flatten()

        # basic circle
        angles = np.linspace(0, 2 * np.pi)
        Xs_basic, Xps_basic = np.cos(angles), np.sin(angles)
        xs_basic = tw_crystal['betx'] ** 0.5 * Xs_basic
        xps_basic = 1 / tw_crystal['betx'] ** 0.5 * Xps_basic - tw_crystal['alfx'] / tw_crystal['betx'] ** 0.5 * Xs_basic

        inside_pct = [
            self._points_inside(
            a * xs_basic + d * tw_crystal['dx'],
            a * xps_basic + d * tw_crystal['dpx'],
            d_crystal, delta_crystal, 
            d0=0,
            delta0=delta0_crystal
            )
            for a, d in zip(aa_flat, dd_flat)
        ]

        # Track the example particles that fall within the crystal

        aaa, xxx = np.meshgrid(aa_flat, xs_basic)
        ddd, _ = np.meshgrid(dd_flat, xs_basic)
        _, xxxp = np.meshgrid(aa_flat, xps_basic)
        aaa, xxx, ddd, xxxp, = aaa.flatten(), xxx.flatten(), ddd.flatten(), xxxp.flatten()

        particles = pd.DataFrame({'x': xxx, 'xp': xxxp, 'a': aaa, 'dpp': ddd})
        particles['x'] = particles.x * particles.a + particles.dpp * tw_crystal['dx']
        particles['xp'] = particles.xp * particles.a + particles.dpp * tw_crystal['dpx']
        pinside = particles[(particles.x > -d_crystal) &
                    (particles.x < d_crystal) &
                    (particles.xp > -delta_crystal + delta0_crystal) &
                    (particles.xp < delta_crystal + delta0_crystal) &
                    (particles.a <= amp_beam)]  # Limit only to where the beam is

        m11 = np.sqrt(tw_septum['betx'] / tw_crystal['betx']) * (np.cos(mux_rad) + tw_crystal['alfx'] * np.sin(mux_rad))
        m12 = np.sqrt(tw_crystal['betx'] * tw_septum['betx']) * np.sin(mux_rad)
        m21 = 1 / np.sqrt(tw_crystal['betx'] * tw_septum['betx']) * ((tw_crystal['alfx'] - tw_septum['alfx']) * np.cos(mux_rad) -
                                         (1 + tw_crystal['alfx'] * tw_septum['alfx']) * np.sin(mux_rad))
        m22 = np.sqrt(tw_crystal['betx'] / tw_septum['betx']) * (np.cos(mux_rad) - tw_septum['alfx'] * np.sin(mux_rad))

        # kicked particles
        xs_centred, xps_centred = pinside.x - pinside.dpp * tw_crystal['dx'], pinside.xp - pinside.dpp * tw_crystal['dpx'] + kick_crystal
        xs_centred_f, xps_centred_f = m11 * xs_centred + m12 * xps_centred, m21 * xs_centred + m22 * xps_centred
        xs_f, xps_f = xs_centred_f + pinside.dpp * tw_septum['dx'], xps_centred_f + pinside.dpp * tw_septum['dpx']

        # No kick version
        nxs_centred, nxps_centred = pinside.x - pinside.dpp * tw_crystal['dx'], pinside.xp - pinside.dpp * tw_crystal['dpx']
        nxs_centred_f, nxps_centred_f = m11 * nxs_centred + m12 * nxps_centred, m21 * nxs_centred + m22 * nxps_centred
        nxs_f, nxps_f = nxs_centred_f + pinside.dpp * tw_septum['dx'], nxps_centred_f + pinside.dpp * tw_septum['dpx']


        # Plotting

        f_stein, ax_stein = plt.subplots(figsize=(14, 4))

        self._plot_boundaries(ax_stein, dpps, tw_crystal, d_crystal, delta_crystal, d0_crystal, delta0_crystal)
        if self.stein_args.cscale == 'log':
            m = ax_stein.tricontourf(dd_flat, aa_flat, np.log10(np.clip(inside_pct, a_min=1, a_max=100)), vmin=0, vmax=2)
            cax = f_stein.colorbar(m, label='Fraction \n inside crystal (%)', ticks=[0, 1, 2],
                 format=ticker.FixedFormatter(['1', '10', '100']))
        else:
            m = ax_stein.tricontourf(dd_flat, aa_flat, inside_pct)
            cax = f_stein.colorbar(m, label='Fraction \n inside crystal (%)')

        ax_stein.set_ylim(0, amps.max())
        ax_stein.set_xlim(-2 * dpp_beam + dpp_offset, 2 * dpp_beam + dpp_offset)

        print(dpp_offset)

        ax_stein.set_ylabel(r'$A_x$ ($\sqrt{m}$)')
        ax_stein.set_xlabel(r'$\delta p / p$ (1)')
        ax_stein.fill_between((dpp_beam * np.sign(dpp_offset) + dpp_offset, dpp_offset), 0,
               amp_beam,
               label=f'beam, {self.stein_args.nsigma_emit}' + '$A_{rms}$', color='grey')

        # Plot points in the steinbach
        amps_plot = [0, 0.5, 1]
        dpps_plot = [0, 1]
        for amp_plot in amps_plot:
            for dpp_plot in dpps_plot:
                ax_stein.scatter(dpp_offset + np.sign(dpp_offset) * dpp_beam * dpp_plot,
                    amp_beam * amp_plot, color='black')

        ax_stein.set_title(f's={tw_crystal["s"]} m')
        ax_stein.legend(loc='upper right')

        f_crystal, ax_crystal = plt.subplots(figsize=(14, 4))

        ax_crystal.plot([-d_crystal, d_crystal, d_crystal, -d_crystal, -d_crystal],
               [-delta_crystal + delta0_crystal, -delta_crystal + delta0_crystal, delta_crystal + delta0_crystal,
            delta_crystal + delta0_crystal, -delta_crystal + delta0_crystal],
               color='fuchsia', label='crystal')

        ax_crystal.plot([0, 0], [delta0_crystal, delta0_crystal + self.crystal.bending_angle], color='pink', marker='s')

        # Plot points in phase space
        for amp_plot in amps_plot:
            for dpp_plot in dpps_plot:
                plotters.draw_ellipse(tw_crystal['alfx'], tw_crystal['betx'], (amp_plot*amp_beam) ** 2, ax=ax_crystal,
                        x0=((dpp_offset + np.sign(dpp_offset) * dpp_beam * dpp_plot) * tw_crystal['dx'],
                        (dpp_offset + np.sign(dpp_offset) * dpp_beam * dpp_plot) * tw_crystal['dpx']))

        ax_crystal.axvline(d_crystal, color='fuchsia', linestyle='dotted')
        ax_crystal.axvline(-d_crystal, color='fuchsia', linestyle='dotted')
        ax_crystal.set_xlabel('x (m)')
        ax_crystal.set_ylabel('xp (rad)')
        ax_crystal.legend()
        ax_crystal.set_title('at crystal')
        
        c = ax_crystal.scatter(pinside.x, pinside.xp, c=pinside.dpp, marker='.')
        f_crystal.colorbar(c, ax=ax_crystal, label='dpp')
        f_crystal.tight_layout()

        f0, (ax0, ax0_tfs) = plt.subplots(1, 2, figsize=(15, 4))
        for amp_plot in amps_plot:
            for dpp_plot in dpps_plot:
                plotters.draw_ellipse(tw_septum['alfx'], tw_septum['betx'], (amp_beam * amp_plot) ** 2, ax=ax0,
                        x0=((dpp_offset + np.sign(dpp_offset) * dpp_beam * dpp_plot) * tw_septum['dx'],
                        (dpp_offset + np.sign(dpp_offset) * dpp_beam * dpp_plot) * tw_septum['dpx']))

        c = ax0.scatter(xs_f, xps_f, c=pinside.dpp, marker='.')
        ax0.scatter(nxs_f, nxps_f, c=pinside.dpp, marker='x')
        f0.colorbar(c, ax=ax0, label='dpp')

        ax0.set_xlabel('x (m)')
        ax0.set_ylabel('xp (rad)')
        ax0.set_title(f'at {self.stein_args.septum_loc}')

        df_plot = tw_pandas[(tw_pandas.s > tw_crystal['s'] - 150) & (tw_pandas.s < tw_crystal['s'] + 150)]
        ax0_tfs.plot(df_plot.s, df_plot.dx, alpha=0.5, color='blue')
        ax_betx = ax0_tfs.twinx()
        ax_betx.plot(df_plot.s, df_plot.betx, color='red', alpha=0.5)
        ax0_tfs.set_ylabel(r'$D_x$ (m)', color='blue')
        ax_betx.set_ylabel(r'$\beta_x$ (m)', color='red')
        ax0_tfs.axvline(tw_crystal['s'], color='fuchsia')
        ax0_tfs.set_xlabel('s (m)')
        f0.tight_layout()

        return (f_stein, f_crystal, f0), (ax_stein, ax_crystal, ax0, ax0_tfs)