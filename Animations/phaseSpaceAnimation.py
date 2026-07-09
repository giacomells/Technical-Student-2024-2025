"""
Turn-by-turn phase-space animation of SPS slow extraction.

Loads the SPS extraction line (``sps_for_sx.json``), configures the
extraction tunes and sextupoles, generates a Gaussian particle distribution,
and displays a matplotlib animation of the phase-space evolution over
successive turns.

Run ``save_sequence_SPS.py`` first if ``sps_for_sx.json`` is absent.

Usage
-----
    python phaseSpaceAnimation.py
"""

from __future__ import annotations

import typing as t
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import xobjects as xo
import xpart as xp
import xtrack as xt

# ---------------------------------------------------------------------------
# Physics constants
# ---------------------------------------------------------------------------

N_EX: float = 10e-6   # normalised horizontal emittance [m·rad]
N_EY: float = 5e-6    # normalised vertical emittance   [m·rad]
DPP: float = 1e-4     # relative momentum spread

septum_aperture_size: float = 68e-3  # ZS blade position [m]

# Default path to the serialised SPS line
_DEFAULT_LINE_JSON = Path(__file__).resolve().parent / "sps_for_sx.json"


# ---------------------------------------------------------------------------
# SeptumAperture
# ---------------------------------------------------------------------------


class SeptumAperture:
    """Simple aperture model that kills particles reaching the ZS wire.

    Particles whose horizontal coordinate exceeds *first_wire_position*
    are declared lost (``state = -1``).  No kick is applied.

    Parameters
    ----------
    first_wire_position : float
        Radial position of the inner ZS wire [m].  Default is 68 mm.
    """

    def __init__(self, first_wire_position: float = 68e-3) -> None:
        self.first_wire_position = first_wire_position

    def interact(self, particles: xp.Particles) -> t.Optional[t.Dict]:
        """Kill particles that have reached or crossed the wire."""
        n_part = particles._num_active_particles
        particles.state[:n_part] = np.where(
            particles.x[:n_part] >= self.first_wire_position, -1, 1
        )
        return None


# ---------------------------------------------------------------------------
# Line configuration
# ---------------------------------------------------------------------------


def configure_line(json_path: Path | None = None) -> xt.Line:
    """Load the SPS extraction line and insert the ZS aperture element.

    Parameters
    ----------
    json_path : Path or None
        Path to ``sps_for_sx.json``.  Defaults to the file in the same
        directory as this script.

    Returns
    -------
    xt.Line
        Line with the ZS aperture installed, cycled to start at the
        septum, and tracker built.
    """
    if json_path is None:
        json_path = _DEFAULT_LINE_JSON

    line = xt.Line.from_json(str(json_path))

    septum = xt.LimitRect(
        min_x=-1.0, max_x=septum_aperture_size, min_y=-1.0, max_y=1.0
    )
    line.insert_element(name="zs_aperture", element=septum, index="ap.up.zs21633")
    line.cycle("zs_aperture", inplace=True)
    line.build_tracker()
    return line


def match_extraction_tunes(line: xt.Line) -> None:
    """Match the extraction tunes (Qx ≈ 26.667, Qy ≈ 26.58) and enable sextupoles.

    Adjusts the quadrupole and chromaticity knobs so that the working
    point is near the third-integer resonance, then switches on the
    extraction sextupoles.

    Parameters
    ----------
    line : xt.Line
        Line returned by :func:`configure_line`.
    """
    line.vars["extr_bump_knob"] = 0.88

    opt = line.match(
        solve=False,
        vary=[
            xt.VaryList(["kqf", "kqd"], step=1e-7),
            xt.VaryList(["qph_setvalue", "qpv_setvalue"], step=1e-4),
        ],
        targets=[
            xt.TargetSet(qx=26.666666666, qy=26.58, tol=1e-5),
            xt.TargetSet(dqx=-1 * 26.0, dqy=0.47 * 26.0, tol=1e-3),
        ],
    )
    opt.solve()

    line.vars["sps_on_extraction"] = 1.0


def generate_particles(
    line: xt.Line,
    tw: xt.TwissTable,
    n_part: int = 1000,
) -> xp.Particles:
    """Build a Gaussian particle ensemble matched to the Twiss optics.

    Parameters
    ----------
    line : xt.Line
    tw : xt.TwissTable
        Twiss table obtained from ``line.twiss(...)``.
    n_part : int
        Number of macro-particles.  Default is 1000.

    Returns
    -------
    xp.Particles
        Particle object ready for tracking.
    """
    x_norm, px_norm = xp.generate_2D_gaussian(num_particles=n_part)
    y_norm, py_norm = xp.generate_2D_gaussian(num_particles=n_part)
    dpp = np.random.rand(n_part) * DPP

    particles = line.build_particles(
        method="4d",
        zeta=0.0,
        delta=dpp,
        x_norm=x_norm,
        px_norm=px_norm,
        y_norm=y_norm,
        py_norm=py_norm,
        nemitt_x=N_EX,
        nemitt_y=N_EY,
    )
    return particles


def track_and_collect(
    line: xt.Line,
    particles: xp.Particles,
    num_turns: int = 100,
) -> tuple[list, list]:
    """Track *particles* for *num_turns* and return per-turn coordinates.

    Parameters
    ----------
    line : xt.Line
    particles : xp.Particles
    num_turns : int
        Number of revolution turns.

    Returns
    -------
    positions_x : list of np.ndarray
    momenta_px  : list of np.ndarray
    """
    line.discard_tracker()
    line.build_tracker(_context=xo.ContextCpu())

    positions_x: list = []
    momenta_px: list = []
    for _ in range(num_turns):
        line.track(particles, num_turns=1)
        positions_x.append(particles.x.copy())
        momenta_px.append(particles.px.copy())

    return positions_x, momenta_px


def build_animation(
    positions_x: list,
    momenta_px: list,
    num_turns: int,
) -> animation.FuncAnimation:
    """Create a matplotlib FuncAnimation of the phase-space evolution.

    Parameters
    ----------
    positions_x, momenta_px : list of np.ndarray
        Per-turn horizontal coordinates and momenta from :func:`track_and_collect`.
    num_turns : int
        Total number of frames.

    Returns
    -------
    matplotlib.animation.FuncAnimation
    """
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.15)
    ax.set_xlim(0.0, 0.1)
    ax.set_ylim(-0.002, 0.0)
    scatter = ax.scatter([], [], s=1)

    def update(frame: int):
        ax.clear()
        ax.set_xlim(0.0, 0.1)
        ax.set_ylim(-0.002, 0.0)
        ax.axvline(septum_aperture_size, color="red", ls="--", label="Septum Aperture")
        ax.scatter(positions_x[frame], momenta_px[frame], s=1)
        ax.set_xlabel("X Position (m)")
        ax.set_ylabel("Px (rad)")
        ax.text(
            0.95, 0.95, f"Turn {frame + 1}",
            ha="right", va="top", transform=ax.transAxes,
        )
        ax.legend(loc="lower left")
        return (scatter,)

    return animation.FuncAnimation(fig, update, frames=num_turns, repeat=True)



