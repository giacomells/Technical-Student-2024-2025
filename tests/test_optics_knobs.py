import importlib.util
import math
import os
import sys
from contextlib import contextmanager
from pathlib import Path

import pytest

pytest.importorskip("xtrack")
pytest.importorskip("xpart")
pytest.importorskip("xobjects")
pytest.importorskip("xcoll")


@contextmanager
def optics_studies_context():
    repo_root = Path(__file__).resolve().parents[1]
    optics_dir = repo_root / "Optics studies"
    previous_cwd = Path.cwd()
    sys.path.insert(0, str(optics_dir))
    os.chdir(optics_dir)
    try:
        yield optics_dir
    finally:
        os.chdir(previous_cwd)
        sys.path.pop(0)


def load_module(module_name, module_path):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_optics_modules():
    with optics_studies_context() as optics_dir:
        elements = load_module("optics_elements", optics_dir / "elements.py")
        sys.modules.setdefault("elements", elements)
        optimisers = load_module("optics_optimisers", optics_dir / "optimisers.py")
    return elements, optimisers


def build_q22_line():
    elements, _ = load_optics_modules()
    with optics_studies_context():
        line = elements.initialise_lineQ22(change_aperture=False)
    line.twiss_default["method"] = "4d"
    line.build_tracker()
    return line


def test_lss4_x_knob_can_move_beam_at_teca_entry():
    _, optimisers = load_optics_modules()
    line = build_q22_line()

    knob = optimisers.set_x_knobLSS4(line)
    knob.step(20)

    twiss = line.twiss(method="4d")
    x_at_teca = twiss.rows["TECA.entry"].x[0]
    px_at_teca = twiss.rows["TECA.entry"].px[0]

    assert math.isclose(x_at_teca, 1e-3, abs_tol=2e-5)
    assert math.isclose(px_at_teca, 0.0, abs_tol=2e-6)


def test_lss4_px_knob_can_change_angle_without_offset():
    _, optimisers = load_optics_modules()
    line = build_q22_line()

    knob = optimisers.set_px_knobLSS4(line)
    knob.step(20)

    twiss = line.twiss(method="4d")
    x_at_teca = twiss.rows["TECA.entry"].x[0]
    px_at_teca = twiss.rows["TECA.entry"].px[0]

    assert math.isclose(px_at_teca, 1e-6, abs_tol=2e-7)
    assert math.isclose(x_at_teca, 0.0, abs_tol=2e-5)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "--color=yes", "-v"]))