# Slow Extraction with Crystals (SPS)

This repository contains studies and simulation workflows for slow extraction from the SPS toward the North Area, including non-resonant extraction strategies with TECA crystals.

The project combines notebooks for analysis/visualization and Python modules used in extraction studies.

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

Python 3.10+ is recommended.

Install dependencies:

```bash
pip install -r requirements.txt
```

`xsuite`/`xtrack`/`xcoll` are required by several studies and notebooks.

## Quick start

1. Clone the repository.
2. Create and activate a virtual environment.
3. Install dependencies with `pip install -r requirements.txt`.
4. Run tests:

```bash
pytest
```

5. Open notebooks for study-specific workflows.

## Reproducibility notes

- Avoid hardcoded local absolute paths when running scripts/notebooks.
- Keep simulation, analysis and plotting as separate steps.
- Document machine/OS constraints if needed for specific workflows.
- If random processes are involved, expose and store random seeds.

## Where to start

For a practical entry point, start from:

- `Non-Resonant Extraction/RecordingInteraction2.0.ipynb`
- `Resonant Extraction example/SPSExtractionWithCrystalSeptum.ipynb`

## Testing

Tests live in `tests/` and are executed with `pytest`.

Current tests validate core configuration consistency from the resonant extraction setup module.

## Exam submission checklist

Before submission:

- [ ] Clean clone works on a second machine.
- [ ] `pip install -r requirements.txt` succeeds.
- [ ] Main study notebook/script run path is documented.
- [ ] `pytest` runs and passes.
- [ ] No hardcoded local absolute paths are required.
- [ ] Documentation is reachable from this README.
- [ ] Commit messages are clear and atomic.
- [ ] No secrets or credentials are committed.
