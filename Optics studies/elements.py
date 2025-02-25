import numpy as np
import pandas as pd
import scipy.stats as stats
import typing as t

import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from shapely.geometry import Polygon

import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc

from cpymad.madx import Madx

p = 400.0  # beam momentum (GeV/c)
momentum = p  # beam momentum (GeV/c)
Brho = p * 3.3356  # beam rigidity ???

N_EX = 10e-6
N_EY = 5e-6
BANDWIDTH = 1e-4

proton_mass_GeV = xt.PROTON_MASS_EV * 1e-9
beam_energy_GeV = 400
gamma = beam_energy_GeV / proton_mass_GeV
EX = N_EX / gamma
EY = N_EY / gamma
DPP = 1.5e-3

class SeptumInteraction:
    def __init__(self, blade_position: float = 68e-3, thickness:float = 0.3e-3, kick:float = 1e-3) -> None:
        self.blade_position = blade_position
        self.thickness = thickness
        self.kick = kick
        pass

    def interact(self, particles: xp.Particles) -> t.Optional[t.Dict]:
        n_part = particles._num_active_particles

        # lose the particles on the blade
        particles.state[:n_part] = np.where((particles.x[:n_part]>self.blade_position) & (particles.x[:n_part]<(self.blade_position+self.thickness)) & (particles.state[:n_part]==1), -1, 1)

        # kick the particles beyond the blade
        if self.kick != 0:
            particles.px[:n_part] +=  np.where((particles.x[:n_part] > (self.blade_position+self.thickness)) & (particles.state[:n_part]==1), 1, 0) * self.kick

        return None


# install the 5 zs
def install_septa(line, install_zs=True, septum_aperture_size=68e-3):
    septa_names_with_apertures = []

    if install_zs:
        septum_names = ["zs.21633", "zs.21639", "zs.21655", "zs.21671", "zs.21676"]
        
        for septum_name in septum_names:
            zs = xt.BeamInteraction(
                length=0.0,
                interaction_process=SeptumInteraction(
                    blade_position=septum_aperture_size, thickness=0.3e-3, kick=0.440e-3 / 5
                ),
            )
            line.insert_element(septum_name + ".sep", zs, index=septum_name)
            septa_names_with_apertures.append(septum_name + ".sep")

        zs = xt.BeamInteraction(
            length=0.0,
            interaction_process=SeptumInteraction(
                blade_position=40e-3, thickness=5.2e-3, kick=0
            ),
        )
        line.insert_element("tpst.21760_entry" + ".sep", zs, index="tpst.21760_entry")
        septa_names_with_apertures.append("tpst.21760_entry" + ".sep")


    for mst in ["mst.21774", "mst.21779", "mst.21794"]:
        zs = xt.BeamInteraction(
            length=0.0,
            interaction_process=SeptumInteraction(
                blade_position=40e-3, thickness=5.2e-3, kick=1.69520713e-3 / 3
            ),
        )
        line.insert_element(mst + ".sep", zs, index=mst)
        septa_names_with_apertures.append(mst + ".sep")

    for mse in ["mse.21832", "mse.21837", "mse.21852", "mse.21857", "mse.21872"]:
        zs = xt.BeamInteraction(
            length=0.0,
            interaction_process=SeptumInteraction(
                blade_position=40e-3, thickness=20e-3, kick=9.74519477e-3 / 5
            ),
        )
        line.insert_element(mse + ".sep", zs, index=mse)
        septa_names_with_apertures.append(mse + ".sep")

    septum = xt.LimitRect(min_x=-1.0, max_x=septum_aperture_size, min_y=-1.0, max_y=1.0)


    line.insert_element(
        name="tt20.extraction",
        element=xt.LimitRect(min_x=-1.0, max_x=70e-3, min_y=-1.0, max_y=1.0),
        index="ap.do.mse21872",
    )
    return septa_names_with_apertures



def draw_synoptic(ax, line, line_df):
    
    compound_names = line_df["name"].unique()
    for compound_name in compound_names:
        compound_df = line_df[line_df["name"] == compound_name]
        if compound_df["element_type"].isin(["Quadrupole"]).any():
            k1 = (
                compound_df[compound_df["element_type"] == "Quadrupole"]["element"]
                .squeeze()
                .k1
            )
            s1, s2 = compound_df["s"].min(), compound_df["s"].max()
            _ = ax.add_patch(
                mpl.patches.Rectangle(
                    (s1, 0), s2 - s1, np.sign(k1), facecolor="k", edgecolor="k"
                )
            )
        elif compound_df["element_type"].isin(["Sextupole", 'Bend', 'Multipole']).any():
            s1, s2 = compound_df["s"].min(), compound_df["s"].max()
            _ = ax.add_patch(
                mpl.patches.Rectangle(
                    (s1, -1), s2 - s1, 2, facecolor="k", edgecolor="k"
                )
            )

def plot_twiss(fig, twiss, line):

    gs = mpl.gridspec.GridSpec(3, 1, height_ratios=[1, 4, 4])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2], sharex=ax1)     #UNUSED, IMPLEMENT TO HAVE THE VISUAL REPRESENTATION OF THE ENVELOPE THOUGH THE LATTICE

    plt.setp(ax2.get_xticklabels(), visible=False)

    # top plot is synoptic
    ax1.axis('off')
    ax1.set_ylim(-1.2, 1)
    ax1.plot([0, twiss['s'].max()], [0, 0], 'k-')


    #2nd plot is beta functions
    ax2.set_ylabel(r'Twiss (m)')
    ax2.plot(twiss['s'], twiss['betx'], 'r-', label=r'$\beta_x$')
    ax2.plot(twiss['s'], twiss['bety'], 'b-', label=r'$\beta_y$')
    ax2.plot(twiss['s'], twiss['dx']*10, 'g-', label=r'$D_x x10$')
    ax2.set_xlim(twiss['s'][0], twiss['s'][-1])
    ax2.legend(loc='upper right')

    line_df = line.to_pandas()
    line_df = line_df[(line_df['s'] >= twiss.s[0]) & (line_df['s'] <= twiss.s[-1])]

    draw_synoptic(ax1, line, line_df)
    
    axnames = ax1.twiny()
    axnames.spines['top'].set_visible(False)
    axnames.spines['left'].set_visible(False)
    axnames.spines['right'].set_visible(False)
    ax1.get_shared_x_axes()._grouper.join(ax1, axnames)
    
    ticks, ticks_labels = list(), list()
    
    for keyword in ['Quadrupole', 'Sextupole', 'Bend', 'Multipole']:
        sub_line = line_df[line_df['element_type'] == keyword]
        ticks += list(sub_line['s'])
        ticks_labels += list(sub_line['name'])

    axnames.set_xticks(ticks)
    axnames.set_xticklabels(ticks_labels, rotation=60)




TECA = xc.EverestCrystal(
    length=2e-3, 
    material=xc.materials.SiliconCrystal, 
    bending_angle = 174e-6 ,
    side="left",
    lattice="strip",
    jaw = - 52.4e-3,
    tilt = - 1.23e-3,
    width = 1.8e-3,
    height = 50e-3
    )