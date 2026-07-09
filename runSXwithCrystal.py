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
imparts a channeling kick.  Particles that clear the MST blade
(x < MST_X0 - MST_DX/2, state=-42) are extracted; particles hitting
the blade (|x - MST_X0| < MST_DX/2, state=-1) are flagged as lost.

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

import crystal_extraction.extraction as extraction


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

    ste, tw0 = extraction.build_simulation(line_path, args.n_particles)
    rmat_to_mst, rmat_to_cry, tw_crystal = extraction.compute_rmatrices(ste, tw0)
    particles, p_init = extraction.run_tracking(
        ste, rmat_to_mst, rmat_to_cry, tw_crystal, args.extraction_time
    )
    extraction.report_efficiency(particles, args.n_particles)

    extraction.section("Saving plots")
    extraction.save_plots(ste, particles, p_init, tw0, output_dir)

    print("\n" + "=" * 64)
    print(f"  Done.  Plots saved to: {output_dir.resolve()}")
    print("=" * 64 + "\n")


if __name__ == "__main__":
    main()
