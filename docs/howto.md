# How-to Guide

## How to run the resonant extraction study

1. Enter `Resonant Extraction example/`.
2. Open `SPSExtractionWithCrystalSeptum.ipynb` in Jupyter.
3. Run all cells top-to-bottom.  The notebook loads `sps_with_extraction_sliced_quads.json`, installs septa and tracks particles.

## How to run the non-resonant extraction study

1. Enter `Non-Resonant Extraction/`.
2. Open `NRSXQ22.ipynb`.
3. Run all cells.  The notebook calls `initialise_lineQ22()` which loads `lhc_q22.json` and installs the TECA crystal.

## How to explore optics

1. Enter `Optics studies/`.
2. Open `MaxDXatTECA.ipynb` for a scan of maximum dispersion at the TECA position versus tune.
3. For a pre-computed result load `twiss_df_15_05_2025.csv`; for a fresh scan run the simulation cells first.

> **Note:** `Apertures.ipynb` in `Optics studies/` references `lhc_q23.json` (Q22.99 optics) which is not tracked in this repository.  Skip or substitute that sub-study.

## How to run the phase-space animation

```bash
cd Animations
python phaseSpaceAnimation.py
```

The script reads `sps_for_sx.json` from the same directory.

## How to run tests

From the repository root:

```bash
pytest
```

Expected output: all tests pass in < 1 s.

## How to add a new test

1. Create `tests/test_<module>.py`.
2. Follow the pattern in `tests/test_septa_constants.py` — use `importlib.util` to load module files from study directories.
3. Run `pytest` to confirm before committing.
