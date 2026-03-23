# MD studies

This directory contains notebooks and local data used for machine-development data inspection.

## Contents

- `view_dataSPS.ipynb`
- `view_dataLSS4.ipynb`
- `datablmALPS/` local pickle data used by the SPS view notebook
- `orbit.csv`

## Usage notes

- `view_dataSPS.ipynb` now loads its pickle files using repository-relative paths under `datablmALPS/`.
- Keep the `datablmALPS/` directory in place if you want the notebook to run without edits.
- Avoid reintroducing absolute local filesystem paths in these notebooks.
