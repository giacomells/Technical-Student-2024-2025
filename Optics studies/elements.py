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

# Proton mass in GeV/c^2
proton_mass_GeV = xt.PROTON_MASS_EV * 1e-9

# Beam energy in GeV
beam_energy_GeV = 400

# Compute gamma
gamma = beam_energy_GeV / proton_mass_GeV

EX = N_EX / gamma
EY = N_EY / gamma

deltaP_P = 1.5e-3

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
def install_septa(line, install_zs=True, septum_aperture_size=68e-3):#
    septa_names_with_apertures = []

    if install_zs:
        septum_names = ["zs.21633", "zs.21639", "zs.21655", "zs.21671", "zs.21676"]
        
        for septum_name in septum_names:
            zs = xt.BeamInteraction(
                length=0.0,
                interaction_process=SeptumInteraction(
                    blade_position=septum_aperture_size, thickness=0.6e-3, kick=0.440e-3 / 5
                ),
            )
            line.insert_element(septum_name + ".sep", zs, index=septum_name)
            septa_names_with_apertures.append(septum_name + ".sep")
        
    tpst = xt.BeamInteraction(
                length=0.0,
                interaction_process=SeptumInteraction(
                    blade_position=0.0394, thickness=0.0045, kick=0
                ))
    line.insert_element("tpst.21760_entry" + ".sep", tpst, index="tpst.21760_entry")
    septa_names_with_apertures.append("tpst.21760_entry" + ".sep")

    i = 0
    for mst in ["mst.21774", "mst.21779", "mst.21794"]:
        mst_positions = [0.04079,0.0424, 0.04401]
        zs = xt.BeamInteraction(
            length=0.0,
            interaction_process=SeptumInteraction(
                blade_position=mst_positions[i], thickness=4.2e-3, kick=1.69520713e-3 / 3 * 2  ## REMOVE THE *2 TO HAVE THE REAL KICK
            ),
        )
        line.insert_element(mst + ".sep", zs, index=mst)
        septa_names_with_apertures.append(mst + ".sep")
        i =+ 1
    
    i = 0
    for mse in ["mse.21832", "mse.21837", "mse.21852", "mse.21857", "mse.21872"]:
        mse_positions = [0.05270, 0.05184, 0.05511, 0.06250, 0.074]
        zs = xt.BeamInteraction(
            length=0.0,
            interaction_process=SeptumInteraction(
                blade_position=mse_positions[i], thickness=17.25e-3, kick=9.74519477e-3 / 5
            ),
        )
        line.insert_element(mse + ".sep", zs, index=mse)
        septa_names_with_apertures.append(mse + ".sep")
        i =+ 1

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


def print_optics_features(line):
    tw_init = line.twiss(method='4d')
    print(tw_init.qx)
    print(tw_init.qy)

    tw0 = tw_init.to_pandas()
    tw0.index = tw0.name  

    betx_teca = tw0.loc['TECA.entry'].betx
    dx_teca = tw0.loc['TECA.entry'].dx

    # COMPUTING THE BEAM SIZE 
    tw0['sx_mm'] = np.sqrt(tw0['betx'] * EX + (tw0['dx'] * deltaP_P)**2) * 1e3
    tw0['sy_mm'] = np.sqrt(tw0['bety'] * EY + (tw0['dy'] * deltaP_P)**2) * 1e3
    
    B_term0 = betx_teca * EX
    D_term0 = (dx_teca * deltaP_P)**2
    ratio0 = D_term0 / B_term0

    print(f'D term / B term = {ratio0:.2f}')
    print(f"sigma_x max = {tw0.sx_mm.max()} mm")
    print(f"sigma_y max = {tw0.sy_mm.max()} mm")

    # DISPERSION AND BETX AT THE TECA LOCATION
    betx_teca = tw0.loc['TECA.entry'].betx
    dx_teca = tw0.loc['TECA.entry'].dx

    # COMPUTING THE BEAM SIZE 
    tw0['sx_mm'] = np.sqrt(tw0['betx'] * EX + (tw0['dx'] * deltaP_P)**2) 
    tw0['sy_mm'] = np.sqrt(tw0['bety'] * EY + (tw0['dy'] * deltaP_P)**2)
    
    B_term0 = betx_teca * EX
    D_term0 = (dx_teca * deltaP_P)**2
    ratio0 = D_term0 / B_term0

    sigma_xMAXQ22 = tw0.sx_mm.max()
    sigma_yMAXQ22 = tw0.sy_mm.max()


    print(f'D term / B term = {ratio0:.2f}')
    print(f"sigma_x max = {tw0.sx_mm.max()} m")
    print(f"sigma_y max = {tw0.sy_mm.max()} m")

    mu_x_teca = tw0.loc['TECA.entry'].mux
    mu_x_tpst = tw0.loc['tpst.21760_entry'].mux
    mu_x_tcsm = tw0.loc['tcsm.51932.'].mux

    phaseAdvanceTecaTpst = mu_x_tpst - mu_x_teca
    phaseAdvanceTecaTcsm = mu_x_tcsm - mu_x_teca
    phaseAdvanceTpstTcsm = mu_x_tcsm - mu_x_tpst

    print(f"Phase advance Teca - Tpst: {phaseAdvanceTecaTpst:.2f}")
    print(f"Phase advance Teca - Tcsm: {phaseAdvanceTecaTcsm:.2f}")
    print(f"Phase advance Tpst - Tcsm: {phaseAdvanceTpstTcsm:.2f}")

"""
TECA = xc.EverestCrystal(
    length=2e-3, 
    material=xc.materials.SiliconCrystal, 
    bending_angle = 174e-6 ,
    side="left",
    lattice="strip",
    jaw = - 51.4e-3,  #original setting
    tilt = - 1.23e-3,
    width = 1.8e-3,
    height = 50e-3
    )
 """
TECA = xc.EverestCrystal(
    length=2e-3, 
    material=xc.materials.SiliconCrystal, 
    bending_angle = 174e-6 ,
    side="left",
    lattice="strip",
    jaw = - 35e-3,  #original setting
    tilt = - 1.67e-3,
    width = 1.8e-3,
    height = 50e-3
    )

TECS = xc.EverestCrystal(
    length=2e-3, 
    material=xc.materials.SiliconCrystal, 
    bending_angle = 174e-6 ,
    side="left",
    lattice="strip",
    jaw = -0.0315,
    tilt = 10e-6,
    width = 2e-3,
    height = 50e-3
    )


# Given values
x_teca = TECA.jaw  
delta_x_teca_prime = TECA.tilt + TECA.bending_angle 



import numpy as np
import xtrack as xt


# Function that opens the apertures that ae blocking the particles
def open_blocking_apertures(line, TECA, deltaP_P):
    """
    Tracks both non-channeled and channeled particles through the given line 
    and adjusts apertures if needed.

    Parameters:
    line: xtrack.Line
        The beamline through which the particles are tracked.
    TECA: object
        An object containing jaw and tilt parameters for the initial particle.
    deltaP_P: float
        Momentum deviation of the particle.
    
    Returns:
    blocking_elements: list
        A list of elements that blocked the particles, along with the x positions.
    """
    
    # Storage for elements blocking the particle
    blocking_elements = []
    
    # Function to track particles and update apertures
    def track_particle(particle_type):
        while True:
            # Build tracker
            line.build_tracker()

            # Initialize particle based on type
            if particle_type == "non-channeled":
                particles = line.build_particles(
                    method='4d',
                    x=TECA.jaw,
                    px=TECA.tilt,
                    y=0.0,
                    py=0.0,
                    zeta=0.0,
                    delta=deltaP_P,
                    mass0=xt.PROTON_MASS_EV,
                    p0c=400e9
                )
            else:  # Channeled particle
                particles = line.build_particles(
                    method='4d',
                    x=TECA.jaw,
                    px=TECA.tilt + TECA.bending_angle,
                    y=0.0,
                    py=0.0,
                    zeta=0.0,
                    delta=deltaP_P,
                    mass0=xt.PROTON_MASS_EV,
                    p0c=400e9
                )

            # Track the particle
            line.track(particles, num_turns=1, turn_by_turn_monitor="ONE_TURN_EBE")

            # Check if the particle reached the end
            if particles.state[0] >= 1 and particles.at_element[0] == line.get_table()['name', '_end_point']:
                print(f"{particle_type.capitalize()} particle successfully reached the end.")
                break  # Exit loop

            # Check if particle is lost
            if np.any(particles.state < 1):
                lost_turn = np.where(particles.state < 1)[0][0]
                lost_element_index = particles.at_element[lost_turn]
                lost_element_name = line.element_names[lost_element_index]
                print(f"{particle_type.capitalize()} particle lost at element {lost_element_name} (index {lost_element_index}) on turn {lost_turn}")

                # Store blocking element details
                blocking_elements.append((lost_element_name, particles.x[lost_turn]))
                element = line[lost_element_name]

                # Adjust aperture if possible
                if isinstance(element, xt.LimitEllipse):
                    print(f"Adjusting {lost_element_name}: expanding LimitEllipse aperture.")
                    element.a *= 2  # Increase semi-major axis
                    element.b *= 2  # Increase semi-minor axis
                elif hasattr(element, 'min_x') and hasattr(element, 'max_x'):
                    print(f"Adjusting {lost_element_name}: increasing rectangular aperture.")
                    element.min_x *= 2  # Increase minimum x aperture
                    element.max_x *= 2  # Increase maximum x aperture
                else:
                    print(f"Element {lost_element_name} does not have adjustable apertures.")
                    break  # Exit loop to avoid infinite loop
            else:
                print(f"{particle_type.capitalize()} particle not lost during tracking.")
                break
    
    # Track both types of particles
    track_particle("non-channeled")
    track_particle("channeled")

    # Print and return blocking elements
    print("Blocking elements encountered:")
    for elem, x_pos in blocking_elements:
        print(f"Element: {elem}, X position: {x_pos}")

    return blocking_elements


def remove_ZS_apertures(line):
    """
    Modifies the apertures of elements within the specified range in the beamline.
    
    Parameters:
    line: xtrack.Line
        The beamline containing the elements to be modified.
    """
    # Extract the list of elements within the given range
    inside_range = False  # Flag to track when we are in the range

    # Iterate through all elements in the line
    for element_name in line.element_names:
        
        # Start modifying elements only when we reach "zs.21633"
        if element_name == "ap.up.zs21633_aper":
            inside_range = True
        
        # Stop modifying elements after "mse.21872"
        if inside_range:
            element = line[element_name]

            # Check if the element has adjustable apertures
            if hasattr(element, 'min_x') and hasattr(element, 'max_x'):
                print(f"Opening aperture for element: {element_name}")
                element.min_x = -1.0  # Set minimum x aperture to -1
                element.max_x = 1.0   # Set maximum x aperture to 1
            else:
                print(f"Element {element_name} does not have adjustable apertures.")
        
        # Exit once we pass "mse.21872"
        if element_name == "qda.21910":
            break  # Stop iterating beyond this point
        
        
def remove_inner_sideLimits_closeTECA(line):
    # Start modifying elements only when we reach "zs.21633"
        # Extract the list of elements within the given range
    inside_range = False  # Flag to track when we are in the range

        # Iterate through all elements in the line
    for element_name in line.element_names:
    
        if element_name == "qecd.31402_aper":
            inside_range = True
        
        
        if inside_range:
            element = line[element_name]

            # Check if the element has adjustable apertures
            if isinstance(element, xt.LimitEllipse):
                element.a = 1 # Increase semi-major axis
            elif hasattr(element, 'min_x') and hasattr(element, 'max_x'):
                element.min_x = -1  # Increase minimum x aperture
        
        if element_name == "ap.do.mse41891_aper":
            break  # Stop iterating beyond this point


def save_df_Limit_elements_features(line):
    # Initialize lists to store the data
    element_names = []
    max_x_values = []
    min_x_values = []
    positions = []

    # Iterate through all elements in the line
    for element_name, element in line.element_dict.items():
        if element_name != "tt20.extraction" and hasattr(element, 'max_x') and hasattr(element, 'min_x'):
            element_names.append(element_name)
            max_x_values.append(element.max_x)
            min_x_values.append(element.min_x)
            positions.append(line.get_table()['s', element_name])
        elif element_name != "tt20.extraction" and isinstance(element, xt.LimitEllipse):
            element_names.append(element_name)
            max_x_values.append(element.a)
            min_x_values.append(-element.a)
            positions.append(line.get_table()['s', element_name])

    # Create a DataFrame
    df_elements = pd.DataFrame({
        'Position': positions,
        'Element Name': element_names,
        'max_x': max_x_values,
        'min_x': min_x_values
    })
    # Display the DataFrame
    print(df_elements)
    return df_elements
    
def save_horizontal_positions_at_septa(recordNONCH, line, septa_names_with_apertures):
    # Extract recorded element indices and x positions
    recorded_elements = recordNONCH.at_element[0, :]  # Convert to 1D array
    recorded_x = recordNONCH.x[0, :]  # Convert to 1D array

    # Dictionary to store results
    x_positions = {}

    # Iterate through each septum element and find its x position
    for sep in septa_names_with_apertures:
        if sep != "tt20.extraction":
            try:
                # Find the index of the element in the recorded elements
                idx = np.where(recorded_elements == line.element_names.index(sep))[0]

                if len(idx) > 0:
                    x_positions[sep] = recorded_x[idx[0]]  # Extract first occurrence
                else:
                    x_positions[sep] = np.nan  # If not found, return NaN
            except ValueError:
                x_positions[sep] = np.nan  # If element is missing, return NaN

    # Convert to DataFrame for readability
    df_x_positions_at_Septum = pd.DataFrame(list(x_positions.items()), columns=["Element", "X Position"])
    # Display the results
    print(df_x_positions_at_Septum)
    
    
