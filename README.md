# Slow Extraction with Crystals (SPS)

This repository contains studies and simulation workflows for slow extraction from the SPS toward the North Area, including new [Non Resonant Slow Extraction Techniques with Crystals.](docs/TechnicalReport2025SlowExtractionWithCrystal.pdf).

The project combines notebooks for analysis/visualization and Python modules used in extraction studies.

## Documentation

Full documentation lives in `docs/`:

- [Theory and background](docs/theory.md) — physics of slow extraction and crystal channeling
- [API reference](docs/api.md) — functions and constants in `elements.py` / `optimisers.py`
- [How-to guide](docs/howto.md) — step-by-step instructions for running each study

  The full study is reported in this article:
  [Technical Report 2025 — Slow Extraction with Crystals](docs/TechnicalReport2025SlowExtractionWithCrystal.pdf)


## Repository structure

- `Non-Resonant Extraction/`: non-resonant extraction studies, elements and optimizers.
- `Resonant Extraction example/`: resonant extraction constants and notebooks.
- `Optics studies/`: optics optimization studies.
- `OrthogonalKnobsLSS4/`: LSS4 orthogonal knob studies.
- `Crystal Analysis/`: crystal-focused analysis notebooks.
- `Animations/`: animation scripts and sequence helpers.
- `MD studies/`: machine development related studies.
- `tests/`: automated tests.

## Requirements


Xsuite is the core dependency.  The recommended installation method is via the
official Xsuite guide (includes Miniforge/conda setup):

<https://xsuite.readthedocs.io/en/latest/installation.html>

Alternatively, install core dependencies directly:

```bash
pip install -r requirements-core.txt
```

To get an environment that perfectly matches the original study setup(done on MacBook with Intel chip), use the
complete pinned freeze:

```bash
pip install -r requirements.txt
```

`xsuite`/`xtrack`/`xcoll`/`xpart` are required by all study notebooks and modules.

Building a separated environment is strongly reccomended in order to avoid conflicts amongst dependencies. 

### Use of Jupyter notebooks

Jupyter notebooks are used throughout this project for analysis and visualization. If Jupyter is not already installed in your system, install it using:

```bash
pip install jupyter
```

Once installed, notebooks can be launched with:

```bash
jupyter notebook
```

## Quick start

1. Clone the repository.
2. Create and activate a virtual environment(preferably from <https://xsuite.readthedocs.io/en/latest/installation.html>).
3. Install dependencies.
4. Run tests:

```bash
pytest
```

5. Open notebooks for study-specific workflows.


## Where to start

For a practical entry point, start from:

- `Animations` folder to get along with XSuite and accelerator physics. Save the sequence and then run the phase space animation.

- To better understand how the accelerator is represented by the .json files, pull the gitlab repository of the SPS from the following link: <https://gitlab.cern.ch/acc-models/acc-models-sps.git>. Once the folder is downloaded, use the file 'BuildSequenceTemplate.py' to build a sequence with specific strengths. WARNING: if you did not create an environment with the requirements, the script won't work.

- `Resonant Extraction example` to get more deep knowledge and example about resonant extraction from the SPS.

## Testing

Tests live in `tests/` and are executed with `pytest`.

Current tests validate:
- Resonant extraction septa constants (`test_septa_constants.py`)
- Non-resonant extraction physics constants and septum interface (`test_elements.py`)

  ### Usage notes
#### Notebook code expects local imports such as `from elements import *` and `from optimisers import *`. In each directory the '.py' files have specific values inside of that could change the result of a simulation displayed in the notebook.
If the shared `.py` or `.json` files are moved or renamed, the notebooks could fail.




