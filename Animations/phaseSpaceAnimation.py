# Import necessary modules
import typing as t  # Module for type hinting
import matplotlib.pyplot as plt  # For plotting data
import numpy as np  # For numerical operations
import xobjects as xo  # For handling objects related to particle tracking
import xpart as xp  # For particle definition and manipulation
import xtrack as xt  # For particle tracking and simulation
import matplotlib.animation as animation  # For creating animations

# Define constants for the simulation
N_EX = 10e-6  # Normalized emittance in the x direction
N_EY = 5e-6  # Normalized emittance in the y direction
DPP = 1e-4  # Momentum spread

# Set the septum aperture size
septum_aperture_size = 68e-3  # Size of the septum aperture in meters

# Load the simulation line from a JSON file
line = xt.Line.from_json("/sps_for_sx.json") #Replace with the pattern of the .json file if not in the same directory

# Define slicing strategies for thick elements, if needed
Strategy = xt.slicing.Strategy
Teapot = xt.slicing.Teapot

# Define a class for the septum aperture
class SeptumAperture:
    def __init__(self, first_wire_position: float = 68e-3) -> None:
        # Initialize the first wire position of the septum aperture
        self.first_wire_position = first_wire_position
        pass

    def interact(self, particles: xp.Particles) -> t.Optional[t.Dict]:
        # Method to define how particles interact with the septum aperture
        n_part = particles._num_active_particles  # Get the number of active particles

        # Update the state of each particle based on its x position relative to the first wire
        particles.state[:n_part] = np.where(
            particles.x[:n_part] >= self.first_wire_position, -1, 1
        )
        return None

# Define a rectangular limit for the septum aperture
septum = xt.LimitRect(min_x=-1., max_x=septum_aperture_size,
                      min_y=-1., max_y=1.)  # Rectangle defined by min/max values

# Insert the defined septum element into the line at a specified index
line.insert_element(name="zs_aperture", element=septum, index='ap.up.zs21633')
line.cycle("zs_aperture", inplace=True)  # Cycle the line to ensure proper setup
line.build_tracker()  # Build the tracker for the simulation

# Switch on extraction bump
line.vars["extr_bump_knob"] = 0.88

# Extraction tunes
opt = line.match(solve=False,
    vary=[
        xt.VaryList(['kqf', 'kqd'], step=1e-7),   # Varying quadrupole focal strengths
        xt.VaryList(['qph_setvalue', 'qpv_setvalue'], step=1e-4),   # Varying phase values
    ],
    targets=[
        xt.TargetSet(qx=26.666666666, qy=26.58, tol=1e-5),    # Desired target tunes
        xt.TargetSet(dqx=-1 * 26., dqy=0.47 * 26., tol=1e-3),    # Desired dispersion values
    ])

opt.solve()  # Solve the matching problem to find optimal parameters

# Switch on extraction sextupoles
line.vars["sps_on_extraction"] = 1.0

# Perform a twiss analysis of the line
tw = line.twiss(continue_on_closed_orbit_error=True)

# Specify number of particles for the simulation
n_part = 1000 

# Generate normalized 2D Gaussian distributions for particle coordinates and momenta
x_norm, px_norm = xp.generate_2D_gaussian(num_particles=n_part)
y_norm, py_norm = xp.generate_2D_gaussian(num_particles=n_part)

# Initialize longitudinal coordinate and momentum spread
zeta = 0.0
dpp = np.random.rand(n_part) * DPP

# Build the particle distribution based on the generated parameters
particles = line.build_particles(
    method="4d",
    zeta=zeta,
    delta=dpp,
    x_norm=x_norm,
    px_norm=px_norm,
    y_norm=y_norm,
    py_norm=py_norm,
    nemitt_x=N_EX,   # Normalized emittance in x direction
    nemitt_y=N_EY,   # Normalized emittance in y direction
)

# Discard previous tracker setup to prevent issues
line.discard_tracker()
# Build the tracker, explicitly using CPU context support
line.build_tracker(_context=xo.ContextCpu())  # Use CPU context for compatibility

# Set up for turn-by-turn tracking
num_turns = 100
positions_x = []
momenta_px = []

# Track the particles turn by turn and store their positions
for turn in range(num_turns):
    line.track(particles, num_turns=1)
    positions_x.append(particles.x.copy())
    momenta_px.append(particles.px.copy())

# Close any previous plots and prepare for new figures
plt.close("all")

# Prepare the animation figure
fig, ax = plt.subplots()
fig.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.15)  # Add margins
ax.set_xlim(0.0, 0.1)
ax.set_ylim(-0.002, 0.0)
scatter = ax.scatter([], [], s=1)

# Set up the animation function
def update(frame):
    ax.clear()
    ax.set_xlim(0.0, 0.1)
    ax.set_ylim(-0.002, 0.0)
    ax.axvline(septum_aperture_size, color="red", ls="--", label='Septum Aperture')
    ax.scatter(positions_x[frame], momenta_px[frame], s=1)
    ax.set_xlabel('X Position (m)')
    ax.set_ylabel('Px (rad)')
    ax.text(0.95, 0.95, f'Turn {frame+1}', ha='right', va='top', transform=ax.transAxes)  # Add turn number
    ax.legend(loc='lower left')  # Add legend
    return scatter,

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=num_turns, repeat=True)
#ani.save('animation.mp4', writer='ffmpeg', fps=10)
#ani.save('animation.gif', writer='pillow', fps=10)
plt.show()

