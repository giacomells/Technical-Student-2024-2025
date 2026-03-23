# OrthogonalKnobsLSS4

This directory contains studies for orthogonal knob construction in LSS4 and the shared files used by the related notebooks.

## Shared files (must stay at directory root)

These files are used across notebooks and helper modules and should remain in this folder:

- `elements.py` (shared elements and line-initialization helpers)
- `optimisers.py` (shared knob construction and matching helpers)
- `lhc_q22.json` (Q22 optics sequence)
- `sps_for_sx.json` (SPS sequence used by helper functions)

## Main notebooks

- `OrthogonalKnobs.ipynb` (combined knob study)
- `Orthogonal_x_knob.ipynb` (horizontal knob example)
- `Orthogonal_px_knob.ipynb` (angular knob example)

## Usage notes

- Run notebooks with the working directory set to this folder.
- Notebook code expects local imports such as `from elements import *` and `from optimisers import *`.
- If shared `.py` or `.json` files are moved or renamed, notebook execution may fail.
- Keep these shared files at the directory root so the notebooks remain portable.