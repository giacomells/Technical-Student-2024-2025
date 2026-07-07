'''Useful utilities for the project'''
from math import ceil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm, ticker
from matplotlib.colors import LogNorm

import crystal_extraction.plotters as plotters


def track_with_mat(particles, mat):
    '''Tracks particles with a matrix'''
    ins = particles.state>=0
    (particles.x[ins],
        particles.px[ins],
        particles.y[ins],
        particles.py[ins],
        particles.zeta[ins],
        particles.delta[ins]) = mat@(particles.x[ins],
                                     particles.px[ins],
                                     particles.y[ins],
                                     particles.py[ins],
                                     particles.zeta[ins],
                                     particles.delta[ins])
# Analysis tool



def mux_deg_to_ele(df, ele, mod=True):
    small_mux = df.loc[ele].mux - df.mux
    big_mux = df.loc[ele].mux + (df.mux.max() - df.mux)
    mux_filter = small_mux >= 0
    mux = mux_filter * small_mux + ~mux_filter * big_mux
    if mod:
        return np.mod(mux * 360, 360)
    else:
        return mux * 360


def points_inside(x, xp, d_crystal, Delta_crystal, Delta0=0):
    return np.sum((x > -d_crystal) & (x < d_crystal) & (xp > -Delta_crystal + Delta0) & (
                xp < Delta_crystal + Delta0)) / x.size * 100


def plot_steinbach_with_tracks(s, df_, d_crystal=1e-3 / 2, Delta_crystal=20e-6 / 2, emit_beam=10e-6 / 426,
                               dpp_beam=2e-3, sigma_amp=3, npoints=100, Delta0=None, kick_crystal=160e-6,
                               cscale='log', sep='mst.21774', sign=1):
    '''
    Notes
    -----
    - Chromatic effect on the tune is ignored
    -
    '''

    df = df_.copy()
    if 'name' in df.columns:
        df.index = df.name

    df['mux_to_sep'] = mux_deg_to_ele(df, sep)
    df['gammx'] = (1 + df.alfx ** 2) / df.betx

    amp_beam = np.sqrt(emit_beam)
    p = df.iloc[np.argmin(np.abs(df.s - s))]
    psep = df.loc[sep]
    dpp_offset = sign * (- p['betx'] ** 0.5 * sigma_amp * amp_beam - d_crystal) / p['dx']
    dpp_stopband = sign * np.sqrt(p.betx) / np.abs(p.dx) * sign * sigma_amp * amp_beam

    angular_offset = p['dpx'] * dpp_offset
    if Delta0 is None:
        Delta0 = angular_offset + sign * (sigma_amp * amp_beam * (-p['alfx'] / p['betx'] ** 0.5) + Delta_crystal)

    amps = np.linspace(0, 2 * sigma_amp * amp_beam, npoints)
    dpps = np.linspace(-2 * dpp_beam + dpp_offset,
                       2 * dpp_beam + dpp_offset, npoints)

    positive_pos_amp = 1 / p['betx'] ** 0.5 * (dpps * p['dx'] - d_crystal)
    negative_pos_amp = -1 / p['betx'] ** 0.5 * (dpps * p['dx'] + d_crystal)
    positive_angle_amp = 1 / p['gammx'] ** 0.5 * (dpps * p['dpx'] - Delta_crystal - Delta0)
    negative_angle_amp = -1 / p['gammx'] ** 0.5 * (dpps * p['dpx'] + Delta_crystal - Delta0)

    # basic circle
    angles = np.linspace(0, 2 * np.pi)
    Xs_basic, Xps_basic = np.cos(angles), np.sin(angles)
    xs_basic = p['betx'] ** 0.5 * Xs_basic
    xps_basic = 1 / p['betx'] ** 0.5 * Xps_basic - p['alfx'] / p['betx'] ** 0.5 * Xs_basic

    aa, dd = np.meshgrid(amps, dpps)
    aa_flat, dd_flat = aa.flatten(), dd.flatten()

    inside_pct = [points_inside(a * xs_basic + d * p['dx'], a * xps_basic + d * p['dpx'],
                                d_crystal, Delta_crystal, Delta0=Delta0)
                  for a, d in zip(aa_flat, dd_flat)]

    f, ax = plt.subplots(1, 2, figsize=(14, 4))
    if cscale == 'log':
        m = ax[0].tricontourf(dd_flat, aa_flat, np.log10(np.clip(inside_pct, a_min=1, a_max=100)), vmin=0, vmax=2)
        cax = f.colorbar(m, label='Fraction \n inside crystal (%)', ticks=[0, 1, 2],
                         format=ticker.FixedFormatter(['1', '10', '100']))
    else:
        m = ax[0].tricontourf(dd_flat, aa_flat, inside_pct)
        cax = f.colorbar(m, label='Fraction \n inside crystal (%)')

        # ax.scatter(dd_flat, aa_flat, marker='.', s=1)
    ax[0].plot(dpps, positive_pos_amp, color='blue')
    ax[0].plot(dpps, negative_pos_amp, color='blue', label='Boundary (offset)')
    ax[0].plot(dpps, positive_angle_amp, color='fuchsia')
    ax[0].plot(dpps, negative_angle_amp, color='fuchsia', label='Boundary (angle)')
    ax[0].set_ylim(0, amps.max())
    ax[0].set_ylabel(r'$A_x$ ($\sqrt{m}$)')
    ax[0].set_xlabel(r'$\delta p / p$ (1)')
    ax[0].fill_between((dpp_beam * np.sign(dpp_offset) + dpp_offset, dpp_offset), 0,
                       sigma_amp * amp_beam,
                       label=f'beam, {sigma_amp}' + '$A_{rms}$', color='grey')
    ax[0].set_xlim(-2 * dpp_beam + dpp_offset, 2 * dpp_beam + dpp_offset)

    # Plot points in the steinbach
    amps_plot = [0, 0.5, 1]
    dpps_plot = [0, 1]
    for amp_plot in amps_plot:
        for dpp_plot in dpps_plot:
            ax[0].scatter(dpp_offset + np.sign(dpp_offset) * dpp_beam * dpp_plot,
                          sigma_amp * amp_beam * amp_plot, color='black')

    ax[0].set_title(f's={s} m')
    ax[0].legend(loc='upper right')

    ax[1].plot([-d_crystal, d_crystal, d_crystal, -d_crystal, -d_crystal],
               [-Delta_crystal + Delta0, -Delta_crystal + Delta0, Delta_crystal + Delta0,
                Delta_crystal + Delta0, -Delta_crystal + Delta0],
               color='fuchsia', label='crystal')

    ax[1].plot([0, 0], [Delta0, Delta0 + kick_crystal], color='pink', marker='s')

    # Plot points in phase space
    for amp_plot in amps_plot:
        for dpp_plot in dpps_plot:
            plotters.draw_ellipse(p['alfx'], p['betx'], (sigma_amp * amp_plot) ** 2 * emit_beam, ax=ax[1],
                                  x0=((dpp_offset + np.sign(dpp_offset) * dpp_beam * dpp_plot) * p['dx'],
                                      (dpp_offset + np.sign(dpp_offset) * dpp_beam * dpp_plot) * p['dpx']))

    ax[1].axvline(d_crystal, color='fuchsia', linestyle='dotted')
    ax[1].axvline(-d_crystal, color='fuchsia', linestyle='dotted')
    ax[1].set_xlabel('x (m)')
    ax[1].set_ylabel('xp (rad)')
    ax[1].legend()
    ax[1].set_title('at crystal')
    f.tight_layout()

    f0, (ax0, ax0_norm) = plt.subplots(1, 2, figsize=(15, 4))
    for amp_plot in amps_plot:
        for dpp_plot in dpps_plot:
            plotters.draw_ellipse(psep['alfx'], psep['betx'], (sigma_amp * amp_plot) ** 2 * emit_beam, ax=ax0,
                                  x0=((dpp_offset + np.sign(dpp_offset) * dpp_beam * dpp_plot) * psep['dx'],
                                      (dpp_offset + np.sign(dpp_offset) * dpp_beam * dpp_plot) * psep['dpx']))

    # Tracking
    # ----------

    # Initial conditions:

    aaa, xxx = np.meshgrid(aa_flat, xs_basic)
    ddd, _ = np.meshgrid(dd_flat, xs_basic)
    _, xxxp = np.meshgrid(aa_flat, xps_basic)
    aaa, xxx, ddd, xxxp, = aaa.flatten(), xxx.flatten(), ddd.flatten(), xxxp.flatten()

    particles = pd.DataFrame({'x': xxx, 'xp': xxxp, 'a': aaa, 'dpp': ddd})
    particles['x'] = particles.x * particles.a + particles.dpp * p['dx']
    particles['xp'] = particles.xp * particles.a + particles.dpp * p['dpx']
    pinside = particles[(particles.x > -d_crystal) &
                        (particles.x < d_crystal) &
                        (particles.xp > -Delta_crystal + Delta0) &
                        (particles.xp < Delta_crystal + Delta0) &
                        (particles.a <= amp_beam * sigma_amp)]  # Limit only to where the beam is

    c = ax[1].scatter(pinside.x, pinside.xp, c=pinside.dpp, marker='.')
    f.colorbar(c, ax=ax[1], label='dpp')

    mux_rad = p['mux_to_sep'] * np.pi / 180
    m11 = np.sqrt(psep['betx'] / p['betx']) * (np.cos(mux_rad) + p['alfx'] * np.sin(mux_rad))
    m12 = np.sqrt(p['betx'] * psep['betx']) * np.sin(mux_rad)
    m21 = 1 / np.sqrt(p['betx'] * psep['betx']) * ((p['alfx'] - psep['alfx']) * np.cos(mux_rad) -
                                                   (1 + p['alfx'] * psep['alfx']) * np.sin(mux_rad))
    m22 = np.sqrt(p['betx'] / psep['betx']) * (np.cos(mux_rad) - psep['alfx'] * np.sin(mux_rad))

    # kicked particles
    xs_centred, xps_centred = pinside.x - pinside.dpp * p['dx'], pinside.xp - pinside.dpp * p['dpx'] + kick_crystal
    xs_centred_f, xps_centred_f = m11 * xs_centred + m12 * xps_centred, m21 * xs_centred + m22 * xps_centred
    xs_f, xps_f = xs_centred_f + pinside.dpp * psep['dx'], xps_centred_f + pinside.dpp * psep['dpx']

    # No kick version
    nxs_centred, nxps_centred = pinside.x - pinside.dpp * p['dx'], pinside.xp - pinside.dpp * p['dpx']
    nxs_centred_f, nxps_centred_f = m11 * nxs_centred + m12 * nxps_centred, m21 * nxs_centred + m22 * nxps_centred
    nxs_f, nxps_f = nxs_centred_f + pinside.dpp * psep['dx'], nxps_centred_f + pinside.dpp * psep['dpx']

    c = ax0.scatter(xs_f, xps_f, c=pinside.dpp, marker='.')
    ax0.scatter(nxs_f, nxps_f, c=pinside.dpp, marker='x')
    f0.colorbar(c, ax=ax0, label='dpp')

    ax0.set_xlabel('x (m)')
    ax0.set_ylabel('xp (rad)')
    ax0.set_title(f'at {sep}')

    f1, ax1 = plt.subplots(figsize=(10, 4))

    df_plot = df[(df.s > s - 150) & (df.s < s + 150)]
    ax1.plot(df_plot.s, df_plot.dx, alpha=0.5, color='blue')
    ax_betx = ax1.twinx()
    ax_betx.plot(df_plot.s, df_plot.betx, color='red', alpha=0.5)
    ax1.set_ylabel(r'$D_x$ (m)', color='blue')
    ax_betx.set_ylabel(r'$\beta_x$ (m)', color='red')
    ax1.axvline(s, color='fuchsia')
    ax1.set_xlabel('s (m)')

    return {'Delta0': Delta0, 'angular_offset': angular_offset,
            'dpp_offset': dpp_offset, 'dpp_stopband': dpp_stopband,
            'dd': dd_flat, 'aa': aa_flat, 'x': xs_basic, 'xp': xps_basic,
            'p': p, 'm11': m11, 'm12': m12, 'm21': m21, 'm22': m22,
            'fig': f, 'ax': ax, 'fig0': f0, 'ax0': ax0, 'fig1': f1, 'ax1': ax1}


