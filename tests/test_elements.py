"""
Tests for the shared elements.py module (loaded from Non-Resonant Extraction/).

These tests validate physics constants and the SeptumInteraction interface
without requiring a full accelerator line to be built.
"""
import importlib.util
import math
from pathlib import Path


def load_elements_module():
    repo_root = Path(__file__).resolve().parents[1]
    elements_file = repo_root / "Non-Resonant Extraction" / "elements.py"
    spec = importlib.util.spec_from_file_location("elements_module", elements_file)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# Physics constants
# ---------------------------------------------------------------------------

def test_beam_momentum_is_400_GeV():
    el = load_elements_module()
    assert el.p == 400.0


def test_beam_rigidity_is_positive_and_matches_formula():
    el = load_elements_module()
    # Brho = p [GeV/c] * 3.3356  T·m/(GeV/c)
    expected = 400.0 * 3.3356
    assert math.isclose(el.Brho, expected, rel_tol=1e-9)


def test_emittances_are_positive():
    el = load_elements_module()
    assert el.N_EX > 0
    assert el.N_EY > 0
    assert el.N_EX < el.N_EY  # vertical emittance is larger in SPS


def test_momentum_spread_is_small_positive():
    el = load_elements_module()
    assert 0 < el.deltaP_P < 1e-2


def test_gamma_is_large_relativistic():
    el = load_elements_module()
    # At 400 GeV/c, gamma >> 1 (proton mass ~0.938 GeV)
    assert el.gamma > 100


def test_geometric_emittances_are_smaller_than_normalised():
    el = load_elements_module()
    # EX = N_EX / gamma  =>  EX << N_EX
    assert el.EX < el.N_EX
    assert el.EY < el.N_EY


# ---------------------------------------------------------------------------
# SeptumInteraction class interface
# ---------------------------------------------------------------------------

def test_septum_default_blade_position():
    el = load_elements_module()
    s = el.SeptumInteraction()
    assert math.isclose(s.blade_position, 68e-3, rel_tol=1e-9)


def test_septum_default_thickness_is_submillimetre():
    el = load_elements_module()
    s = el.SeptumInteraction()
    assert 0 < s.thickness < 1e-3


def test_septum_custom_parameters_stored_correctly():
    el = load_elements_module()
    s = el.SeptumInteraction(blade_position=50e-3, thickness=0.5e-3, kick=2e-3)
    assert math.isclose(s.blade_position, 50e-3)
    assert math.isclose(s.thickness, 0.5e-3)
    assert math.isclose(s.kick, 2e-3)


def test_septum_blade_position_is_positive():
    el = load_elements_module()
    s = el.SeptumInteraction()
    assert s.blade_position > 0


def test_septum_thickness_less_than_blade_position():
    """Blade thickness must be smaller than the blade position itself."""
    el = load_elements_module()
    s = el.SeptumInteraction()
    assert s.thickness < s.blade_position
