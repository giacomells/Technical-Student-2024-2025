"""
Tests for Animations/save_sequence_SPS.

The fast unit tests (no network, no heavy computation) verify the module
structure — constants, function signatures, return types, and path handling.

Integration tests that actually download the SPS model from CERN GitLab or
read large JSON files are marked ``@pytest.mark.slow`` and are skipped by
default.  Run them explicitly with ``pytest -m slow``.
"""

import inspect
from pathlib import Path

import pytest

requests = pytest.importorskip("requests")
pytest.importorskip("xtrack")
pytest.importorskip("xpart")
pytest.importorskip("xobjects")

import Animations.save_sequence_SPS as seq_mod  # noqa: E402

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SPS_JSON = _REPO_ROOT / "Animations" / "sps_for_sx.json"


# ---------------------------------------------------------------------------
# save_sequence_SPS — module-level constants
# ---------------------------------------------------------------------------


def test_sps_model_urls_is_a_list():
    """SPS_MODEL_URLS must be a list (not a tuple or other iterable)."""
    assert isinstance(seq_mod.SPS_MODEL_URLS, list)


def test_sps_model_urls_has_five_entries():
    """Five model files are needed to build the Q26 SPS sequence."""
    assert len(seq_mod.SPS_MODEL_URLS) == 5


def test_sps_model_urls_are_https():
    """All SPS model URLs must use HTTPS for secure transport."""
    assert all(url.startswith("https://") for url in seq_mod.SPS_MODEL_URLS)


def test_sps_model_urls_are_strings():
    """Every entry in SPS_MODEL_URLS must be a non-empty string."""
    assert all(isinstance(url, str) and url for url in seq_mod.SPS_MODEL_URLS)


def test_momentum_constant_is_400_GeV():
    """Beam momentum in save_sequence_SPS must be 400 GeV/c."""
    assert seq_mod.p == 400.0


def test_rigidity_constant_matches_formula():
    """Magnetic rigidity constant must equal p × 3.3356 T·m."""
    import math
    assert math.isclose(seq_mod.Brho, 400.0 * 3.3356, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# save_sequence_SPS — function signatures
# ---------------------------------------------------------------------------


def test_build_sps_madx_is_callable():
    """build_sps_madx must be a callable function."""
    assert callable(seq_mod.build_sps_madx)


def test_build_sps_madx_accepts_momentum_parameter():
    """build_sps_madx must accept a momentum_gev_c keyword argument."""
    sig = inspect.signature(seq_mod.build_sps_madx)
    assert "momentum_gev_c" in sig.parameters


def test_save_sps_json_is_callable():
    """save_sps_json must be a callable function."""
    assert callable(seq_mod.save_sps_json)


def test_save_sps_json_accepts_output_path_and_momentum():
    """save_sps_json must accept output_path and momentum_gev_c parameters."""
    sig = inspect.signature(seq_mod.save_sps_json)
    assert "output_path" in sig.parameters
    assert "momentum_gev_c" in sig.parameters


# ---------------------------------------------------------------------------
# Pre-built JSON files (if available)
# ---------------------------------------------------------------------------


def test_sps_json_is_valid_json_when_present():
    """When sps_for_sx.json exists it must be parseable JSON."""
    if not _SPS_JSON.exists():
        pytest.skip("sps_for_sx.json not present — run save_sequence_SPS.py first")
    import json
    with open(_SPS_JSON) as f:
        data = json.load(f)
    assert isinstance(data, dict), "Top-level JSON structure must be a dict"


# ---------------------------------------------------------------------------
# Integration tests (network required) — skipped by default
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_can_build_sps_sequence_from_cern_gitlab():
    """Smoke test: download SPS model, convert, verify line structure."""
    Madx = pytest.importorskip("cpymad.madx").Madx
    xt = pytest.importorskip("xtrack")
    xp = pytest.importorskip("xpart")

    mad = seq_mod.build_sps_madx(momentum_gev_c=400.0)
    line = xt.Line.from_madx_sequence(
        mad.sequence["sps"], deferred_expressions=True, allow_thick=True
    )
    line.particle_ref = xt.Particles(
        mass0=xp.PROTON_MASS_EV, gamma0=mad.sequence["sps"].beam.gamma
    )
    line.twiss_default["method"] = "4d"

    assert line.particle_ref is not None
    assert line.get_length() > 6000
    assert "actcse.31632" in line.element_names
    assert "ap.up.zs21633" in line.element_names


@pytest.mark.slow
def test_rf_setup_supports_6d_twiss_after_sequence_creation():
    """After enabling the RF cavity a valid 6-D Twiss solution must exist."""
    xt = pytest.importorskip("xtrack")
    xp = pytest.importorskip("xpart")

    mad = seq_mod.build_sps_madx(momentum_gev_c=400.0)
    line = xt.Line.from_madx_sequence(
        mad.sequence["sps"], deferred_expressions=True, allow_thick=True
    )
    line.particle_ref = xt.Particles(
        mass0=xp.PROTON_MASS_EV, gamma0=mad.sequence["sps"].beam.gamma
    )

    line.vv["v200"] = 0.0
    line.vv["freq200"] = 200e6
    line.vv["lag200"] = 180.0
    line.element_refs["actcse.31632"].voltage = line.vars["v200"]
    line.element_refs["actcse.31632"].frequency = line.vars["freq200"]
    line.element_refs["actcse.31632"].lag = line.vars["lag200"]

    line.vv["v200"] = 10e6
    twiss_6d = line.twiss(method="6d")
    assert twiss_6d.qs > 1e-3, "Synchrotron tune must be > 1e-3 for a valid 6-D solution"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "--color=yes", "-v"]))