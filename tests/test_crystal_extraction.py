"""
Tests for the crystal_extraction package.

Covers pure-Python / NumPy functions that require no network connection and no
pre-built lattice file.  Modules that depend on xcoll at import time are guarded
with ``pytest.importorskip`` so the test file is skipped cleanly when xcoll is
not installed.

Modules tested
--------------
crystal_extraction.utils          : helper maths and analysis functions
crystal_extraction.plotters       : covariance-ellipse utilities
crystal_extraction.dummy_crystal  : parametric crystal model (needs xcoll)
crystal_extraction.steinbach      : Steinbach diagram builder (needs xcoll)
crystal_extraction.xsuite_line_creation : lattice creation from CERN GitLab
"""

import inspect
import math

import numpy as np
import pandas as pd
import pytest

# xcoll is required by dummy_crystal and steinbach at module import time.
pytest.importorskip("xcoll")
pytest.importorskip("xtrack")

from crystal_extraction.dummy_crystal import DummyCrystal  # noqa: E402
from crystal_extraction.plotters import draw_ellipse, ellipse_frm_cov  # noqa: E402
from crystal_extraction.steinbach import BeamArgs, SteinArgs, Steinbach  # noqa: E402
from crystal_extraction.utils import mux_deg_to_ele, points_inside  # noqa: E402
import crystal_extraction.xsuite_line_creation as line_mod  # noqa: E402


# ===========================================================================
# crystal_extraction.utils
# ===========================================================================


class TestPointsInside:
    """Unit tests for the ``points_inside`` acceptance-fraction function."""

    def test_all_particles_inside_returns_100_percent(self):
        """When every particle is inside the crystal window the result is 100 %."""
        x = np.zeros(100)
        xp = np.zeros(100)
        result = points_inside(x, xp, d_crystal=1e-3, Delta_crystal=20e-6)
        assert math.isclose(result, 100.0)

    def test_no_particles_inside_returns_0_percent(self):
        """When no particle is inside the window the result is 0 %."""
        x = np.full(100, 10e-3)   # far outside ± 1 mm window
        xp = np.zeros(100)
        result = points_inside(x, xp, d_crystal=1e-3, Delta_crystal=20e-6)
        assert math.isclose(result, 0.0)

    def test_result_is_a_percentage_between_0_and_100(self):
        """Return value must always lie in [0, 100]."""
        rng = np.random.default_rng(42)
        x = rng.uniform(-5e-3, 5e-3, 500)
        xp = rng.uniform(-50e-6, 50e-6, 500)
        result = points_inside(x, xp, d_crystal=1e-3, Delta_crystal=20e-6)
        assert 0.0 <= result <= 100.0

    def test_angular_offset_shifts_accepted_region(self):
        """A non-zero Delta0 shifts the angular acceptance window."""
        x = np.zeros(100)
        # Particles at xp = +30e-6, outside centred window but inside shifted one.
        xp = np.full(100, 30e-6)
        without_offset = points_inside(x, xp, d_crystal=1e-3, Delta_crystal=20e-6, Delta0=0)
        with_offset = points_inside(x, xp, d_crystal=1e-3, Delta_crystal=20e-6, Delta0=20e-6)
        assert with_offset > without_offset

    def test_larger_window_accepts_more_particles(self):
        """Increasing the crystal window must not decrease the accepted fraction."""
        rng = np.random.default_rng(0)
        x = rng.uniform(-3e-3, 3e-3, 300)
        xp = rng.uniform(-60e-6, 60e-6, 300)
        small = points_inside(x, xp, d_crystal=0.5e-3, Delta_crystal=10e-6)
        large = points_inside(x, xp, d_crystal=2e-3, Delta_crystal=40e-6)
        assert large >= small


class TestMuxDegToEle:
    """Unit tests for ``mux_deg_to_ele`` phase-advance converter."""

    @pytest.fixture
    def twiss_df(self):
        """Minimal Twiss-like DataFrame with a monotonically increasing mux."""
        names = ["start", "mid", "target", "end"]
        return pd.DataFrame({"mux": [0.0, 0.1, 0.2, 0.35]}, index=names)

    def test_returns_pandas_series(self, twiss_df):
        """Output must be a pandas Series with the same index as the input."""
        result = mux_deg_to_ele(twiss_df, "target")
        assert isinstance(result, pd.Series)
        assert list(result.index) == list(twiss_df.index)

    def test_values_are_in_0_to_360_when_mod_true(self, twiss_df):
        """With mod=True all values must lie in [0°, 360°)."""
        result = mux_deg_to_ele(twiss_df, "target", mod=True)
        assert (result >= 0).all()
        assert (result < 360).all()

    def test_upstream_element_has_positive_phase_advance(self, twiss_df):
        """Elements upstream of the target must have a positive phase advance."""
        result = mux_deg_to_ele(twiss_df, "target", mod=True)
        # 'start' and 'mid' are upstream of 'target'
        assert result["start"] > 0
        assert result["mid"] > 0


# ===========================================================================
# crystal_extraction.plotters
# ===========================================================================


class TestEllipseFrmCov:
    """Unit tests for ``ellipse_frm_cov``."""

    def test_returns_three_floats(self):
        """Function must return a 3-tuple of floats (width, height, rotation)."""
        cov = np.eye(2)
        width, height, rotation = ellipse_frm_cov(cov)
        assert isinstance(width, float)
        assert isinstance(height, float)
        assert isinstance(rotation, float)

    def test_identity_covariance_gives_equal_width_and_height(self):
        """Unit covariance (circle) must yield equal width and height."""
        cov = np.eye(2)
        width, height, _ = ellipse_frm_cov(cov)
        assert math.isclose(width, height, rel_tol=1e-9)

    def test_nsig_scales_dimensions_linearly(self):
        """Doubling nsig must double both width and height."""
        cov = [[4.0, 1.0], [1.0, 2.0]]
        w1, h1, _ = ellipse_frm_cov(cov, nsig=1.0)
        w2, h2, _ = ellipse_frm_cov(cov, nsig=2.0)
        assert math.isclose(w2, 2 * w1, rel_tol=1e-9)
        assert math.isclose(h2, 2 * h1, rel_tol=1e-9)

    def test_width_and_height_are_positive(self):
        """Ellipse dimensions must always be positive."""
        cov = [[9.0, 3.0], [3.0, 4.0]]
        width, height, _ = ellipse_frm_cov(cov)
        assert width > 0
        assert height > 0

    def test_diagonal_covariance_has_zero_rotation(self):
        """An axis-aligned covariance matrix must produce a zero rotation angle."""
        cov = [[4.0, 0.0], [0.0, 1.0]]
        _, _, rotation = ellipse_frm_cov(cov)
        assert math.isclose(abs(rotation) % 90, 0.0, abs_tol=1e-9)


class TestDrawEllipse:
    """Unit tests for ``draw_ellipse``."""

    def test_returns_matplotlib_axes(self):
        """draw_ellipse must return a matplotlib Axes object."""
        import matplotlib
        matplotlib.use("Agg")  # non-interactive backend for CI
        ax = draw_ellipse(alpha=0.0, beta=50.0, eps=1e-9)
        import matplotlib.pyplot as plt
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_accepts_existing_axes(self):
        """When an Axes is passed it must be reused rather than a new one created."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        _, ax_in = plt.subplots()
        ax_out = draw_ellipse(alpha=0.5, beta=30.0, eps=1e-9, ax=ax_in)
        assert ax_out is ax_in
        plt.close("all")


# ===========================================================================
# crystal_extraction.dummy_crystal
# ===========================================================================


class TestDummyCrystalDefaults:
    """Tests for DummyCrystal default parameter values."""

    @pytest.fixture
    def crystal(self):
        return DummyCrystal()

    def test_default_width_is_1mm(self, crystal):
        """Default crystal half-width dx must be 1 mm."""
        assert math.isclose(crystal.dx, 1e-3)

    def test_default_channeling_probability_is_valid(self, crystal):
        """Channeling probability lambda_chan must be in (0, 1)."""
        assert 0 < crystal.lambda_chan < 1

    def test_default_tracker_is_none(self, crystal):
        """Tracker must be None until explicitly built."""
        assert crystal.tracker is None

    def test_default_bending_radius_is_positive(self, crystal):
        """Bending radius must be a positive number."""
        assert crystal.bending_radius > 0

    def test_default_active_length_is_positive(self, crystal):
        """Active length must be a positive number."""
        assert crystal.active_length > 0


class TestDummyCrystalKickProperty:
    """Tests for the mu_kick_chan computed property."""

    def test_mu_kick_chan_matches_formula(self):
        """mu_kick_chan must equal active_length / bending_radius."""
        crystal = DummyCrystal()
        expected = crystal.active_length / crystal.bending_radius
        assert math.isclose(crystal.mu_kick_chan, expected, rel_tol=1e-9)

    def test_mu_kick_chan_setter_updates_active_length(self):
        """Setting mu_kick_chan must update active_length consistently."""
        crystal = DummyCrystal()
        new_kick = 200e-6
        crystal.mu_kick_chan = new_kick
        assert math.isclose(crystal.active_length, new_kick * crystal.bending_radius, rel_tol=1e-9)

    def test_round_trip_mu_kick_chan(self):
        """Setting and reading back mu_kick_chan must be consistent."""
        crystal = DummyCrystal()
        crystal.mu_kick_chan = 150e-6
        assert math.isclose(crystal.mu_kick_chan, 150e-6, rel_tol=1e-9)

    def test_build_pandas_tracker_returns_callable(self):
        """build_pandas_tracker must return a callable and store it on self."""
        crystal = DummyCrystal()
        tracker = crystal.build_pandas_tracker()
        assert callable(tracker)
        assert crystal.tracker is tracker


# ===========================================================================
# crystal_extraction.steinbach
# ===========================================================================


class TestBeamArgs:
    """Tests for the BeamArgs dataclass defaults."""

    def test_default_emit_rms_is_positive(self):
        """Default normalised emittance must be positive."""
        args = BeamArgs()
        assert args.emit_rms > 0

    def test_default_dpp_fullwidth_is_positive(self):
        """Default full momentum spread must be positive."""
        args = BeamArgs()
        assert args.dpp_fullwidth > 0

    def test_default_nparticles_is_positive_integer(self):
        """Default particle count must be a positive integer."""
        args = BeamArgs()
        assert isinstance(args.nparticles, int)
        assert args.nparticles > 0

    def test_custom_values_stored_correctly(self):
        """Constructor must store custom values without modification."""
        args = BeamArgs(emit_rms=5e-9, dpp_fullwidth=1e-3, nparticles=500)
        assert math.isclose(args.emit_rms, 5e-9)
        assert math.isclose(args.dpp_fullwidth, 1e-3)
        assert args.nparticles == 500


class TestSteinArgs:
    """Tests for the SteinArgs dataclass defaults."""

    def test_default_npoints_is_positive(self):
        """Default number of Steinbach grid points must be positive."""
        args = SteinArgs()
        assert args.npoints > 0

    def test_default_crystal_angular_acceptance_is_positive(self):
        """Crystal angular acceptance must be a positive angle."""
        args = SteinArgs()
        assert args.crystal_angular_acceptance > 0

    def test_default_nsigma_emit_is_positive(self):
        """Default number of sigma for emittance cut must be positive."""
        args = SteinArgs()
        assert args.nsigma_emit > 0

    def test_approach_side_is_plus_or_minus_one(self):
        """approach_side must be +1 or -1."""
        args = SteinArgs()
        assert args.approach_side in (1, -1)


class TestSteinbachStaticMethods:
    """Tests for pure-computation static / internal methods of Steinbach."""

    def test_points_inside_all_inside(self):
        """_points_inside must return 100 when all particles are in window."""
        x = np.zeros(100)
        xp = np.zeros(100)
        result = Steinbach._points_inside(x, xp, d_crystal=1e-3, delta_crystal=20e-6)
        assert math.isclose(result, 100.0)

    def test_points_inside_none_inside(self):
        """_points_inside must return 0 when no particle is in window."""
        x = np.full(100, 10e-3)
        xp = np.zeros(100)
        result = Steinbach._points_inside(x, xp, d_crystal=1e-3, delta_crystal=20e-6)
        assert math.isclose(result, 0.0)

    def test_points_inside_result_in_range(self):
        """_points_inside result must be in [0, 100]."""
        rng = np.random.default_rng(7)
        x = rng.uniform(-5e-3, 5e-3, 200)
        xp = rng.uniform(-50e-6, 50e-6, 200)
        result = Steinbach._points_inside(x, xp, d_crystal=1e-3, delta_crystal=20e-6)
        assert 0.0 <= result <= 100.0

    def test_stein_boundaries_returns_four_arrays(self):
        """_stein_boundaries must return exactly 4 arrays of the same shape."""
        beam_args = BeamArgs()
        stein_args = SteinArgs()
        # Build a minimal Steinbach instance without a real line.
        sb = Steinbach.__new__(Steinbach)
        sb.beam_args = beam_args
        sb.stein_args = stein_args

        dpps = np.linspace(-2e-3, 2e-3, 20)
        tw_crystal = {"betx": 50.0, "dx": 2.0, "dpx": 0.05, "gammx": 0.02}

        result = sb._stein_boundaries(dpps, tw_crystal, d_crystal=1e-3, delta_crystal=20e-6)
        assert len(result) == 4
        for arr in result:
            assert arr.shape == dpps.shape

    def test_stein_boundaries_constant_when_zero_dispersion(self):
        """With zero dispersion (dx=dpx=0) boundaries must be independent of dpps."""
        sb = Steinbach.__new__(Steinbach)
        sb.beam_args = BeamArgs()
        sb.stein_args = SteinArgs()

        dpps = np.linspace(-2e-3, 2e-3, 10)
        tw_crystal = {"betx": 50.0, "dx": 0.0, "dpx": 0.0, "gammx": 0.02}

        pos_pos, neg_pos, pos_ang, neg_ang = sb._stein_boundaries(
            dpps, tw_crystal, d_crystal=1e-3, delta_crystal=20e-6
        )
        # All values must be the same since dpps has no effect without dispersion.
        assert np.allclose(pos_pos, pos_pos[0])
        assert np.allclose(neg_pos, neg_pos[0])
        assert np.allclose(pos_ang, pos_ang[0])
        assert np.allclose(neg_ang, neg_ang[0])


# ===========================================================================
# crystal_extraction.xsuite_line_creation
# ===========================================================================


class TestXsuiteLineCreationModule:
    """Tests for module-level constants and function signatures."""

    def test_default_urls_is_a_list(self):
        """DEFAULT_URLS must be a list."""
        assert isinstance(line_mod.DEFAULT_URLS, list)

    def test_default_urls_has_five_entries(self):
        """Five model files are needed to build the SPS Q20 line."""
        assert len(line_mod.DEFAULT_URLS) == 5

    def test_default_urls_are_https(self):
        """All default URLs must use HTTPS for secure transport."""
        assert all(url.startswith("https://") for url in line_mod.DEFAULT_URLS)

    def test_default_urls_are_non_empty_strings(self):
        """Every URL must be a non-empty string."""
        assert all(isinstance(url, str) and url for url in line_mod.DEFAULT_URLS)

    def test_create_xsuite_line_is_callable(self):
        """create_xsuite_line must be a callable function."""
        assert callable(line_mod.create_xsuite_line)

    def test_create_xsuite_line_accepts_sequence_parameter(self):
        """create_xsuite_line must accept a 'sequence' keyword argument."""
        sig = inspect.signature(line_mod.create_xsuite_line)
        assert "sequence" in sig.parameters

    def test_create_xsuite_line_accepts_file_path_parameter(self):
        """create_xsuite_line must accept a 'file_path' keyword argument."""
        sig = inspect.signature(line_mod.create_xsuite_line)
        assert "file_path" in sig.parameters

    def test_create_xsuite_line_accepts_urls_parameter(self):
        """create_xsuite_line must accept a 'urls' keyword argument."""
        sig = inspect.signature(line_mod.create_xsuite_line)
        assert "urls" in sig.parameters

    @pytest.mark.slow
    def test_create_xsuite_line_builds_sps_from_cern_gitlab(self):
        """Integration: download Q20 model and verify the resulting line."""
        import xtrack as xt
        line = line_mod.create_xsuite_line(sequence="sps")
        assert isinstance(line, xt.Line)
        assert line.get_length() > 6000


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "--color=yes", "-v"]))
