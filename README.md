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

```
Animations/              # Runnable scripts for the phase-space animation workflow
    save_sequence_SPS.py    # Download SPS model from CERN GitLab → sps_for_sx.json
    phaseSpaceAnimation.py  # Load JSON and run turn-by-turn animation
    BuildSequenceTemplate.py  # Build LHC Q22 line from local MAD-X files → lhc_q22.json
    sps_for_sx.json         # Pre-built SPS Q26 extraction line (committed)

docs/                    # Project documentation
    theory.md
    api.md
    howto.md

tests/                   # Automated test suite (pytest)
    test_sps_sequence_creation.py  # Animation script unit + integration tests
```

## Requirements


Xsuite is the core dependency.  The recommended installation method is via the
official Xsuite guide (includes Miniforge/conda setup):

<https://xsuite.readthedocs.io/en/latest/installation.html>

Alternatively, install core dependencies with range-pinned specifiers:

```bash
pip install -r requirements-core.txt
```

For a fully reproducible environment matching the original study setup (macOS,
Intel), use the complete pinned freeze:

```bash
pip install -r requirements.txt
```

## Quick start

1. Clone the repository.
2. Create and activate a virtual environment (Xsuite guide recommended).
3. Install dependencies.
4. Run the test suite:

```bash
pytest
```

5. Generate the SPS extraction line (needed for the animation):

```bash
cd Animations && python save_sequence_SPS.py
```

6. Run the phase-space animation:

```bash
cd Animations && python phaseSpaceAnimation.py
```

### Optional: build the Q22 optics line

The optics integration tests require `lhc_q22.json`, generated from local
MAD-X model files placed in `Animations/acc-models-sps/`:

```bash
cd Animations && python BuildSequenceTemplate.py
```

## Testing

Tests live in `tests/` and are executed with `pytest`:

```bash
pytest                    # run all tests (fast unit + skippable integration)
pytest -m slow            # also run network / heavy integration tests
```

Current tests validate:
- Animation script module structure and constants (`test_sps_sequence_creation.py`)


