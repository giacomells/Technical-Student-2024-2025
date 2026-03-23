# Resonant Extraction example

This directory contains the resonant extraction study notebooks and the shared assets they depend on.

## Shared files (must stay at directory root)

These files are imported/loaded by multiple notebooks and should remain in this folder:

- `septa.py` (shared constants/dictionary for extraction elements)
- `sps_with_extraction.json` (ring sequence)
- `sps_with_extraction_sliced_quads.json` (ring sequence with sliced quadrupoles)

## Main notebooks

- `SPSExtractionWithCrystalSeptum.ipynb`
- `SPSExtractionPlots.ipynb`
- `ToyModel.ipynb`

## Usage notes

- Run the notebooks with the working directory set to this folder.
- Notebook code expects relative paths like `./sps_with_extraction_sliced_quads.json` and `import septa`.
- If files are moved or renamed, notebook imports/loading will fail.

## Quick sanity check

From repository root:

```bash
python - <<'PY'
import json
from pathlib import Path
base = Path("Resonant Extraction example")
json.load(open(base / "sps_with_extraction.json"))
json.load(open(base / "sps_with_extraction_sliced_quads.json"))
print("JSON files are readable")
PY
```



