from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_REPO_ROOT))

import matplotlib
matplotlib.use("Agg")

import crystal_extraction.plotters as plotters

# PARSING INPUT ARGUMENTS
parser = argparse.ArgumentParser(
    description="Crystal single-pass characterisation: efficiency vs kick angle"
)
parser.add_argument(
    "--n-particles",
    type=int,
    default=50_000,
    help="Particles per configuration (default: 50000)",
)
parser.add_argument(
    "--output-dir",
    default="outputs/crystal",
    help="Directory for saved plots (default: outputs/crystal)",
)
args = parser.parse_args()

# RUNNING CRYSTAL CHARACTERISATION
plotters.run_crystal_characterisation(
    n_particles=args.n_particles,
    output_dir=Path(args.output_dir),
)
