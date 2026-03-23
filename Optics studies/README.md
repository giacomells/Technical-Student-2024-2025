# Optics studies

This directory contains optics-focused notebooks and the shared helper files they depend on.

## Shared files (must stay at directory root)

These files are used by multiple notebooks and helper modules and should remain in this folder:

- `elements.py` (shared optics/extraction helper definitions)
- `optimisers.py` (shared matching and optimization helpers)
- `sps_for_sx.json` (SPS optics sequence)
- `lhc_q22.json` (Q22 optics sequence)

The `database/` subdirectory stores copies/reference data, but the notebooks and helper functions are safer when the shared JSON files are also available at the directory root.

## Main notebooks

- `MaxDXatTECA.ipynb` (main tune-scan study)
- `Apertures.ipynb`
- `KnobsComparison.ipynb`
- `multiTurnSeparationTPST.ipynb`

## Usage notes

- Run notebooks with the working directory set to this folder.
- Notebook code uses local imports such as `from elements import ...` and `from optimisers import ...`.
- Some notebooks also refer to files in `database/`; keep those resources available.
- Avoid hardcoded absolute paths when exporting or reloading intermediate data.

