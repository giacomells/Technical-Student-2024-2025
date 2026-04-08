# Non-Resonant Extraction

This directory contains the non-resonant extraction studies and the shared files used by the related notebooks.

## Shared files (must stay at directory root)

These files are reused across notebooks and helper modules and should remain in this folder:

- `elements.py` (shared elements, crystal/septa helpers, line initialization helpers)
- `optimisers.py` (shared matching and optimization helpers)
- `sps_for_sx.json` (SPS sequence used by non-resonant extraction workflows)
- `lhc_q22.json` (Q22 optics sequence)

## Main notebooks

- `NRSXQ22.ipynb`
- `RecordingInteraction2.0.ipynb`

## Usage notes

- Run notebooks with the working directory set to this folder.
- Notebook code expects local imports such as `from elements import *` and `from optimisers import *`.
- Line initialization helpers in `elements.py` load local JSON sequences from this same directory.
- If the shared `.py` or `.json` files are moved or renamed, the notebooks will fail.

