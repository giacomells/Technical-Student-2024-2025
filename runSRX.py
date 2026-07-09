from __future__ import annotations

import typing as t
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import xobjects as xo
import xpart as xp
import xtrack as xt
from Animations import phaseSpaceAnimation as phase_space_animation


line = phase_space_animation.configure_line()
phase_space_animation.match_extraction_tunes(line)
tw = line.twiss(continue_on_closed_orbit_error=True)

particles = phase_space_animation.generate_particles(line, tw, n_part=1000)
num_turns = 100
positions_x, momenta_px = phase_space_animation.track_and_collect(line, particles, num_turns)

plt.close("all")
ani = phase_space_animation.build_animation(positions_x, momenta_px, num_turns)

_gif_path = Path(__file__).resolve().parent / "phase_space_animation.gif"
ani.save(str(_gif_path), writer="pillow", fps=10)
print(f"Animation saved to {_gif_path}")

plt.show()