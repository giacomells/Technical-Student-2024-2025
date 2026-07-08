"""
run_extraction_simulation.py
============================
Simulate non-resonant slow extraction from the SPS driven by a crystal
(matrix-tracking method) and report extraction efficiency to the terminal.

Physics
-------
Particles start at the crystal location covering the Steinbach stopband.
The momentum offset (delta) is ramped turn-by-turn to drive particles into
the crystal angular acceptance.  A 6-D R-matrix propagates coordinates
between the crystal and the MST wire.  The crystal (EverestCrystal/xcoll)
imparts a channeling kick.  Particles landing inside MST_X0 +/- MST_DX/2
are counted as extracted; those that overshoot are flagged separately.

Usage
-----
    python run_extraction_simulation.py
    python run_extraction_simulation.py --n-particles 2000
    python run_extraction_simulation.py --output-dir results/

Requirements
------------
    xsuite + xcoll environment (xsuite_env)
    database/lhc_q20.json  (already in the repo)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_REPO_ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xtrack as xt

import crystal_extraction.plotters as plotters
import crystal_extraction.utils as utils
from crystal_extraction.steinbach import BeamArgs, SteinArgs, Steinbach, sps_crystal

SPS_T_REV = 23e-6    # SPS revolution period at 400 GeV/c [s]
MST_X0    = -7.5e-3  # MST blade centre [m]
MST_DX    =  5.0e-3  # MST blade full width [m]


def _section(title: str) -> None:
    print("\n" + "-" * 64)
    print(f"  {title}")
    print("-" * 64)


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

def build_simulation(line_path: Path, n_particles: int):
    _section("Loading optics")
    print(f"  Line JSON : {line_path}")
    xline = xt.Line.from_json(str(line_path))
    xline.particle_ref = xt.Particles(p0c=400e9, mass0=xt.PROTON_MASS_EV)
    print("  Twiss (4d) ...", end=" ", flush=True)
    tw0 = xline.twiss(method="4d")
    print("done")

    _section("Crystal and beam configuration")
    beam_args  = BeamArgs(nparticles=n_particles)
    stein_args = SteinArgs(approach_side=1, nsigma_emit=3)
    crystal    = sps_crystal.copy()
    crystal.bending_angle = 170e-6

    print(f"  Crystal location  : {stein_args.crystal_loc}")
    print(f"  Septum location   : {stein_args.septum_loc}")
    print(f"  Bending angle     : {crystal.bending_angle * 1e6:.1f} urad")
    print(f"  Emittance (geom.) : {beam_args.emit_rms * 1e9:.2f} nm.rad")
    print(f"  Momentum spread   : {beam_args.dpp_fullwidth * 1e3:.1f} per-mille (full width)")
    print(f"  Particles         : {n_particles:,}")

    ste = Steinbach(line=xline, crystal=crystal,
                    beam_args=beam_args, stein_args=stein_args)

    _section("Generating initial particle distribution")
    print("  create_xsuite_particles ...", end=" ", flush=True)
    ste.create_xsuite_particles()
    print("done")
    return ste, tw0


# ---------------------------------------------------------------------------
# R-matrices
# ---------------------------------------------------------------------------

def compute_rmatrices(ste: Steinbach, tw0):
    _section("Computing R-matrices  (crystal <-> MST)")
    new_line   = ste.line.cycle(name_first_element=ste.stein_args.crystal_loc)
    tw         = new_line.twiss(method="4d")
    tw_pandas  = ste.parse_xsuite_twiss(tw)
    tw_crystal = tw_pandas.loc[ste.stein_args.crystal_loc]
    ste.crystal.tilt = 0

    rmat_to_mst = tw.get_R_matrix(
        start=ste.stein_args.crystal_loc, end=ste.stein_args.septum_loc)
    line_from_mst = new_line.cycle(name_first_element=ste.stein_args.septum_loc)
    tw_from_mst   = line_from_mst.twiss(method="4d")
    rmat_to_cry   = tw_from_mst.get_R_matrix(
        start=ste.stein_args.septum_loc, end=ste.stein_args.crystal_loc)

    print(f"  dpp stopband rms : {tw_crystal['dpp_stopband_rms'] * 1e3:.4f} per-mille")
    return rmat_to_mst, rmat_to_cry, tw_crystal


# ---------------------------------------------------------------------------
# Tracking loop
# ---------------------------------------------------------------------------

def run_tracking(ste: Steinbach, rmat_to_mst, rmat_to_cry,
                 tw_crystal, extraction_time: float = 5.0):
    dpp_dturn = (ste.stein_args.approach_side
                 * ste.beam_args.dpp_fullwidth * SPS_T_REV / extraction_time)
    turns = int(2 * abs(tw_crystal["dpp_stopband_rms"]
                        * ste.stein_args.nstopbands_dpp / dpp_dturn))

    _section("Tracking")
    print(f"  Extraction time  : {extraction_time:.1f} s")
    print(f"  delta per turn   : {dpp_dturn:.3e}")
    print(f"  Number of turns  : {turns:,}\n")

    particles          = ste.xsuite_particles.copy()
    particles_init_mst = ste.xsuite_particles.copy()
    utils.track_with_mat(particles,          rmat_to_mst)
    utils.track_with_mat(particles_init_mst, rmat_to_mst)

    report_every = max(1, turns // 10)
    for turn in range(turns):
        utils.track_with_mat(particles, rmat_to_cry)
        ste.crystal.track(particles)
        utils.track_with_mat(particles, rmat_to_mst)
        particles.delta += dpp_dturn

        overshoot = particles.x < MST_X0 - MST_DX / 2
        extracted = np.abs(particles.x - MST_X0) < MST_DX / 2
        particles.state[overshoot] = -42
        particles.state[extracted] = -1

        if (turn + 1) % report_every == 0:
            alive  = int(np.sum(particles.state == 1))
            n_extr = int(np.sum(particles.state == -1))
            n_over = int(np.sum(particles.state == -42))
            pct    = n_extr / ste.beam_args.nparticles * 100
            print(f"  turn {turn+1:6d}/{turns}  |  alive {alive:5d}  "
                  f"extracted {n_extr:5d}  overshot {n_over:5d}  eff {pct:.1f}%")

    return particles, particles_init_mst


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def report_efficiency(particles, n_total: int) -> None:
    n_extr  = int(np.sum(particles.state == -1))
    n_over  = int(np.sum(particles.state == -42))
    n_alive = int(np.sum(particles.state == 1))
    n_other = n_total - n_extr - n_over - n_alive
    cap_den = n_extr + n_over

    _section("Extraction efficiency summary")
    print(f"  Initial particles         : {n_total:6d}  (100.0%)")
    print(f"  Extracted (MST window)    : {n_extr:6d}  ({n_extr/n_total*100:.1f}%)")
    print(f"  Overshot  (missed MST)    : {n_over:6d}  ({n_over/n_total*100:.1f}%)")
    print(f"  Still circulating         : {n_alive:6d}  ({n_alive/n_total*100:.1f}%)")
    print(f"  Other losses              : {n_other:6d}  ({n_other/n_total*100:.1f}%)")
    print(f"  " + "-" * 44)
    print(f"  Total extraction eff.     : {n_extr/n_total*100:.2f}%")
    if cap_den > 0:
        print(f"  MST capture fraction      : {n_extr/cap_den*100:.2f}%"
              "  (extracted / (extracted+overshot))")


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def save_plots(ste: Steinbach, particles, particles_init_mst,
               tw0, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Steinbach diagram
    print("\n  Steinbach diagram ...", end=" ", flush=True)
    figs, _ = ste.plot(tw0)
    for i, fig in enumerate(figs):
        p = output_dir / f"steinbach_{i}.png"
        fig.savefig(p, dpi=120, bbox_inches="tight")
        plt.close(fig)
    print(f"saved ({len(figs)} file(s))")

    # 2. Phase space at MST: initial vs end of spill
    print("  Phase-space at MST ...", end=" ", flush=True)
    f, ax = plotters.subplots(1, 2, sharex=True, sharey=True)
    factor    = 1e3
    isextract = (particles.state == -1).astype(int)

    ax[0].scatter(particles_init_mst.x * factor,
                  particles_init_mst.px * factor,
                  c=particles_init_mst.state, s=4, marker=".")
    ax[0].set_title("Initial")
    ax[1].scatter(particles.x * factor, particles.px * factor,
                  c=isextract, s=4, marker=".")
    ax[1].set_title("End of spill  (yellow = extracted)")

    for axis in ax:
        axis.axvline(factor * (MST_X0 - MST_DX / 2), color="fuchsia",
                     lw=1.2, label="MST blades")
        axis.axvline(factor * (MST_X0 + MST_DX / 2), color="fuchsia", lw=1.2)
        axis.set_xlabel("x [mm]")
        axis.set_ylabel("xp [mrad]")
        axis.legend(fontsize=7)

    n_extr = int(np.sum(particles.state == -1))
    n_tot  = len(particles.state)
    n_over = int(np.sum(particles.state == -42))
    f.suptitle(
        f"Phase space at MST  |  extracted {n_extr}/{n_tot} "
        f"({n_extr/n_tot*100:.1f}%)  overshot {n_over}",
        fontsize=9,
    )
    f.tight_layout()
    out = output_dir / "phase_space_mst.png"
    f.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(f)
    print(f"saved to {out}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="SPS non-resonant slow extraction with crystal — matrix tracking"
    )
    parser.add_argument("--json-path",
                        default=str(_REPO_ROOT / "database" / "lhc_q20.json"))
    parser.add_argument("--n-particles", type=int, default=1000)
    parser.add_argument("--extraction-time", type=float, default=5.0)
    parser.add_argument("--output-dir", default="outputs/extraction")
    args = parser.parse_args()

    line_path  = Path(args.json_path)
    output_dir = Path(args.output_dir)

    if not line_path.exists():
        print(f"\nERROR: line file not found: {line_path}")
        print("Generate it with crystal_extraction/xsuite_line_creation.py first.\n")
        sys.exit(1)

    print("\n" + "=" * 64)
    print("  SPS Non-Resonant Slow Extraction with Crystal")
    print("=" * 64)

    ste, tw0 = build_simulation(line_path, args.n_particles)
    rmat_to_mst, rmat_to_cry, tw_crystal = compute_rmatrices(ste, tw0)
    particles, p_init = run_tracking(ste, rmat_to_mst, rmat_to_cry,
                                     tw_crystal, args.extraction_time)
    report_efficiency(particles, args.n_particles)

    _section("Saving plots")
    save_plots(ste, particles, p_init, tw0, output_dir)

    print("\n" + "=" * 64)
    print(f"  Done.  Plots saved to: {output_dir.resolve()}")
    print("=" * 64 + "\n")


if __name__ == "__main__":
    main()
