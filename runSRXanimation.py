from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from Animations import phaseSpaceAnimation as phase_space_animation
from Animations.save_sequence_SPS import save_sps_json


# SAVING THE SPS SEQUENCE JSON IF IT DOES NOT EXIST
_REPO_ROOT = Path(__file__).resolve().parent
_SPS_JSON = _REPO_ROOT / "database" / "sps_for_sx.json"

if not _SPS_JSON.exists():
    saved = save_sps_json(output_path=_SPS_JSON)
    print(f"Created SPS sequence JSON at {saved}")

# BUILDINg SEQUENCE, MATCHING EXTRACTION TUNES, AND TRACKING PARTICLES
line = phase_space_animation.configure_line(json_path=_SPS_JSON)
phase_space_animation.match_extraction_tunes(line)
tw = line.twiss(continue_on_closed_orbit_error=True)

# CREATING PARTICLE DISTRIBUTION AND TRACKING
particles = phase_space_animation.generate_particles(line, tw, n_part=1000)
num_turns = 100
positions_x, momenta_px = phase_space_animation.track_and_collect(line, particles, num_turns)

plt.close("all")
ani = phase_space_animation.build_animation(positions_x, momenta_px, num_turns)

# SAVING THE GIF ANIMATION
# _gif_path = Path(__file__).resolve().parent / "phase_space_animation.gif"
# ani.save(str(_gif_path), writer="pillow", fps=10)
# print(f"Animation saved to {_gif_path}")

plt.show()