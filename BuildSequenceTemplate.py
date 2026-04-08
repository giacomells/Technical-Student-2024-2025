from pathlib import Path

from cpymad.madx import Madx
import xtrack as xt


PROJECT_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PROJECT_ROOT / "acc-models-sps"
OUTPUT_JSON = PROJECT_ROOT / "lhc_q22.json"


def build_lhc_q22_json() -> Path:
    mad = Madx()
    mad.call(str(REPO_ROOT / "sps.seq"))
    mad.call(str(REPO_ROOT / "strengths" / "lhc_q22.str"))
    mad.call(str(REPO_ROOT / "beams" / "beam_lhc_injection.madx"))
    mad.use(sequence="sps")

    line = xt.Line.from_madx_sequence(
        mad.sequence.sps,
        deferred_expressions=False,
    )
    line.to_json(str(OUTPUT_JSON))
    return OUTPUT_JSON


if __name__ == "__main__":
    output_path = build_lhc_q22_json()
    print(output_path)