# Slow Extraction with Crystals (SPS)

This repository contains studies and simulation workflows for slow extraction from the SPS toward the North Area, including new [Non Resonant Slow Extraction Techniques with Crystals.](docs/TechnicalReport2025SlowExtractionWithCrystal.pdf).

The project combines notebooks for analysis/visualization and Python modules used in extraction studies.


## Requirements
Xsuite is the core dependency.  The recommended installation method is via the
official Xsuite guide (includes Miniforge/conda setup):

<https://xsuite.readthedocs.io/en/latest/installation.html>

Building a separated environment is strongly reccomended in order to avoid conflicts amongst dependencies. 
### Automated setup (recommended)

Run the provided installer once. It detects your OS and CPU architecture,
downloads Miniforge if needed, creates a dedicated conda environment
(`xsuite_env`), and installs all dependencies:

```bash
chmod +x install.sh
./install.sh
```

Supported platforms: **macOS** (Intel & Apple Silicon) and **Linux**
(x86_64 & aarch64). On **Windows**, install
[WSL](https://learn.microsoft.com/en-us/windows/wsl/install) first and run
the script from the WSL terminal — the script prints step-by-step
instructions if it detects a Windows shell.

After installation, scripts can be run in three ways — no manual environment
selection is needed:

| Scenario | How |
|---|---|
| **VS Code** | `install.sh` writes `.vscode/settings.json` just press ▶ Run |
| **Terminal** (one-time per session) | `conda activate xsuite_env`, then `python script.py` as usual |
| **Terminal** (no activation) | `./run.sh Animations/save_sequence_SPS.py` or `./run.sh -m pytest` |

#### Manual setup

If you already have a Python environment, install the core dependencies
directly (range-pinned specifiers):

```bash
pip install -r requirements-core.txt
```

For a fully reproducible environment matching the original study setup
(macOS, Intel), use the complete pinned freeze:

```bash
pip install -r requirements.txt
```

## Documentation

Full documentation lives in `docs/`:

- [Theory and background](docs/theory.md) — physics of slow extraction and crystal channeling
- [API reference](docs/api.md) — functions and constants in `elements.py` / `optimisers.py`
- [How-to guide](docs/howto.md) — step-by-step instructions for running each study

The full study is reported in this article:
[Technical Report 2025 — Slow Extraction with Crystals](docs/TechnicalReport2025SlowExtractionWithCrystal.pdf)


## Quick start

1. Clone the repository.
2. Run `./install.sh` (see Requirements above).
3. Run the test suite:

```bash
pytest                         # A. after conda activate xsuite_env  
./run.sh -m pytest             # B. without activating
```

4. Run on of the .py files in the home directory:
```bash
./run.sh run_crystal_characterisation.py
./run.sh runSRXanimation.py
./run.sh runSXwithCrystal.py
```


### Testing

Tests live in `tests/` and are executed with `pytest`:

```bash
pytest                    # run all tests (fast unit + skippable integration)
pytest -m slow            # also run network / heavy integration tests
```

