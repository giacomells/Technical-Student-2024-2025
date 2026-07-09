from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xtrack as xt

import crystal_extraction.plotters as plotters
import crystal_extraction.utils as utils
from crystal_extraction.steinbach import BeamArgs, SteinArgs, Steinbach, sps_crystal

SPS_T_REV = 23e-6    # SPS revolution period at 400 GeV/c [s]
MST_X0 = -7.5e-3     # MST blade centre [m]
MST_DX = 5.0e-3      # MST blade full width [m]


def section(title: str) -> None:
    print("\n" + "-" * 64)
    print(f"  {title}")
    print("-" * 64)


def build_simulation(line_path: Path, n_particles: int):
    section("Loading optics")
    print(f"  Line JSON : {line_path}")
    xline = xt.Line.from_json(str(line_path))
    xline.particle_ref = xt.Particles(p0c=400e9, mass0=xt.PROTON_MASS_EV)
    print("  Twiss (4d) ...", end=" ", flush=True)
    tw0 = xline.twiss(method="4d")
    print("done")

    section("Crystal and beam configuration")
    beam_args = BeamArgs(nparticles=n_particles)
    stein_args = SteinArgs(approach_side=1, nsigma_emit=3)
    crystal = sps_crystal.copy()
    crystal.bending_angle = 170e-6

    print(f"  Crystal location  : {stein_args.crystal_loc}")
    print(f"  Septum location   : {stein_args.septum_loc}")
    print(f"  Bending angle     : {crystal.bending_angle * 1e6:.1f} urad")
    print(f"  Emittance (geom.) : {beam_args.emit_rms * 1e9:.2f} nm.rad")
    print(f"  Momentum spread   : {beam_args.dpp_fullwidth * 1e3:.1f} per-mille (full width)")
    print(f"  Particles         : {n_particles:,}")

    ste = Steinbach(line=xline, crystal=crystal,
                    beam_args=beam_args, stein_args=stein_args)

    section("Generating initial particle distribution")
    print("  create_xsuite_particles ...", end=" ", flush=True)
    ste.create_xsuite_particles()
    print("done")
    return ste, tw0


def compute_rmatrices(ste: Steinbach, tw0):
    del tw0  # kept for API compatibility with previous call-site

    section("Computing R-matrices  (crystal <-> MST)")
    new_line = ste.line.cycle(name_first_element=ste.stein_args.crystal_loc)
    tw = new_line.twiss(method="4d")
    tw_pandas = ste.parse_xsuite_twiss(tw)
    tw_crystal = tw_pandas.loc[ste.stein_args.crystal_loc]
    ste.crystal.tilt = 0

    rmat_to_mst = tw.get_R_matrix(
        start=ste.stein_args.crystal_loc, end=ste.stein_args.septum_loc)
    line_from_mst = new_line.cycle(name_first_element=ste.stein_args.septum_loc)
    tw_from_mst = line_from_mst.twiss(method="4d")
    rmat_to_cry = tw_from_mst.get_R_matrix(
        start=ste.stein_args.septum_loc, end=ste.stein_args.crystal_loc)

    print(f"  dpp stopband rms : {tw_crystal['dpp_stopband_rms'] * 1e3:.4f} per-mille")
    return rmat_to_mst, rmat_to_cry, tw_crystal


def run_tracking(ste: Steinbach, rmat_to_mst, rmat_to_cry,
                 tw_crystal, extraction_time: float = 5.0):
    dpp_dturn = (ste.stein_args.approach_side
                 * ste.beam_args.dpp_fullwidth * SPS_T_REV / extraction_time)
    turns = int(2 * abs(tw_crystal["dpp_stopband_rms"]
                        * ste.stein_args.nstopbands_dpp / dpp_dturn))

    section("Tracking")
    print(f"  Extraction time  : {extraction_time:.1f} s")
    print(f"  delta per turn   : {dpp_dturn:.3e}")
    print(f"  Number of turns  : {turns:,}\n")

    particles = ste.xsuite_particles.copy()
    particles_init_mst = ste.xsuite_particles.copy()
    utils.track_with_mat(particles, rmat_to_mst)
    utils.track_with_mat(particles_init_mst, rmat_to_mst)

    report_every = max(1, turns // 10)
    for turn in range(turns):
        utils.track_with_mat(particles, rmat_to_cry)
        ste.crystal.track(particles)
        utils.track_with_mat(particles, rmat_to_mst)
        particles.delta += dpp_dturn

        out = particles.x < MST_X0 - MST_DX / 2
        lost = np.abs(particles.x - MST_X0) < MST_DX / 2
        particles.state[out] = -42
        particles.state[lost] = -1

        if (turn + 1) % report_every == 0:
            alive = int(np.sum(particles.state == 1))
            nout = int(np.sum(particles.state == -42))
            nlost = int(np.sum(particles.state < 0))
            pct = nout / ste.beam_args.nparticles * 100
            print(f"  turn {turn + 1:6d}/{turns}  |  alive {alive:5d}  "
                  f"out {nlost:5d}  well_out {nout:5d}  portion {pct:.1f}%")

    return particles, particles_init_mst


def report_efficiency(particles, n_total: int) -> None:
    isout = particles.state == -42   # cleared the blade (notebook: isout)
    islost = particles.state < 0     # all removed    (notebook: islost)
    nout = int(np.sum(isout))
    nlost = int(np.sum(islost))
    n_alive = int(np.sum(particles.state == 1))
    n_other = n_total - nlost - n_alive

    section("Extraction efficiency summary")
    print(f"  Initial particles         : {n_total:6d}  (100.0%)")
    print(f"  Out   (all removed)       : {nlost:6d}  ({nlost / n_total * 100:.1f}%)")
    print(f"  Well out (past blade)     : {nout:6d}  ({nout / n_total * 100:.1f}%)")
    print(f"  Still circulating         : {n_alive:6d}  ({n_alive / n_total * 100:.1f}%)")
    print(f"  Other losses              : {n_other:6d}  ({n_other / n_total * 100:.1f}%)")
    print("  " + "-" * 44)
    if nlost > 0:
        print(f"  Portion (nout/nlost)      : {nout / nlost * 100:.2f}%"
              "  (well out / all removed)")


def save_plots(ste: Steinbach, particles, particles_init_mst,
               tw0, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n  Steinbach diagram ...", end=" ", flush=True)
    figs, _ = ste.plot(tw0)
    for i, fig in enumerate(figs):
        p = output_dir / f"steinbach_{i}.png"
        fig.savefig(str(p), dpi=120, bbox_inches="tight")
        plt.close(fig)
    print(f"saved ({len(figs)} file(s))")

    print("  Phase-space at MST ...", end=" ", flush=True)
    f, ax = plotters.subplots(1, 2, sharex=True, sharey=True)
    factor = 1e3
    islost = particles.state < 0
    isout = particles.state == -42
    nout = int(np.sum(isout))
    nlost = int(np.sum(islost))

    ax[0].scatter(particles_init_mst.x * factor,
                  particles_init_mst.px * factor,
                  c=particles_init_mst.state, s=4, marker=".", label="initial")
    ax[1].scatter(particles.x * factor, particles.px * factor,
                  c=islost.astype(int), s=4, marker=".", label="mid-spill")

    for axis in ax:
        axis.axvline(factor * (MST_X0 - MST_DX / 2), color="fuchsia",
                     lw=1.2, label="MST 5mm blade")
        axis.axvline(factor * (MST_X0 + MST_DX / 2), color="fuchsia", lw=1.2)
        axis.set_xlabel("x [mm]")
        axis.set_ylabel("xp [mrad]")
        axis.legend(fontsize=7)

    f.suptitle(
        f"Out: {nlost}, Well out: {nout}, Portion: {nout / nlost:.2f}" if nlost > 0
        else "No particles removed",
        fontsize=9,
    )
    f.tight_layout()
    out = output_dir / "phase_space_mst.png"
    f.savefig(str(out), dpi=120, bbox_inches="tight")
    plt.close(f)
    print(f"saved to {out}")
