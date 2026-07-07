"""Tests for the Animations.phaseSpaceAnimation workflow.

These tests exercise the lightweight runtime contracts of the slow-extraction
animation helpers without building the full SPS lattice.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("xobjects")
pytest.importorskip("xpart")
pytest.importorskip("xtrack")

import Animations.phaseSpaceAnimation as phase_mod  # noqa: E402


class FakeParticles:
    """Minimal particle container for fast unit tests."""

    def __init__(self, x, px=None, state=None):
        self.x = np.array(x, dtype=float)
        self.px = np.array(px if px is not None else np.zeros_like(self.x), dtype=float)
        self.state = np.array(state if state is not None else np.ones(len(self.x), dtype=int), dtype=int)
        self._num_active_particles = len(self.x)


class TestSeptumAperture:
    """Behavioral tests for the septum aperture model."""

    def test_interact_marks_particles_at_or_beyond_wire_as_lost(self):
        """Particles reaching the wire must be marked lost."""
        particles = FakeParticles([0.01, 0.068, 0.08], state=[0, 0, 0])

        aperture = phase_mod.SeptumAperture(first_wire_position=0.068)
        result = aperture.interact(particles)

        assert result is None
        np.testing.assert_array_equal(particles.state, np.array([1, -1, -1]))


class TestConfigureLine:
    """Tests for SPS line loading and aperture installation."""

    def test_configure_line_inserts_aperture_cycles_and_builds_tracker(self, monkeypatch, tmp_path):
        """configure_line must install the ZS aperture and prepare the tracker."""
        captured = {}
        sentinel_limit = object()

        class FakeLine:
            def __init__(self):
                self.insert_calls = []
                self.cycled_to = None
                self.tracker_built = False

            def insert_element(self, name, element, index):
                self.insert_calls.append((name, element, index))

            def cycle(self, name, inplace):
                self.cycled_to = (name, inplace)

            def build_tracker(self):
                self.tracker_built = True

        fake_line = FakeLine()

        def fake_from_json(path):
            captured["json_path"] = path
            return fake_line

        def fake_limit_rect(**kwargs):
            captured["limit_rect"] = kwargs
            return sentinel_limit

        fake_xt = SimpleNamespace(
            Line=SimpleNamespace(from_json=fake_from_json),
            LimitRect=fake_limit_rect,
        )
        monkeypatch.setattr(phase_mod, "xt", fake_xt)

        custom_json = tmp_path / "line.json"
        result = phase_mod.configure_line(json_path=custom_json)

        assert result is fake_line
        assert captured["json_path"] == str(custom_json)
        assert captured["limit_rect"] == {
            "min_x": -1.0,
            "max_x": phase_mod.septum_aperture_size,
            "min_y": -1.0,
            "max_y": 1.0,
        }
        assert fake_line.insert_calls == [("zs_aperture", sentinel_limit, "ap.up.zs21633")]
        assert fake_line.cycled_to == ("zs_aperture", True)
        assert fake_line.tracker_built is True

    def test_configure_line_propagates_file_loading_errors(self, monkeypatch):
        """Missing JSON inputs must surface as file-loading failures."""

        def fake_from_json(path):
            raise FileNotFoundError(path)

        fake_xt = SimpleNamespace(
            Line=SimpleNamespace(from_json=fake_from_json),
            LimitRect=lambda **kwargs: None,
        )
        monkeypatch.setattr(phase_mod, "xt", fake_xt)

        with pytest.raises(FileNotFoundError):
            phase_mod.configure_line(json_path=Path("missing.json"))


class TestMatchExtractionTunes:
    """Tests for tune-matching side effects."""

    def test_match_extraction_tunes_sets_knobs_and_solves_match(self, monkeypatch):
        """match_extraction_tunes must solve the optics match and enable extraction."""
        captured = {}

        class FakeOptimizer:
            def __init__(self):
                self.solved = False

            def solve(self):
                self.solved = True

        optimizer = FakeOptimizer()

        class FakeLine:
            def __init__(self):
                self.vars = {}

            def match(self, solve, vary, targets):
                captured["solve"] = solve
                captured["vary"] = vary
                captured["targets"] = targets
                return optimizer

        def fake_vary_list(names, step):
            return {"names": names, "step": step}

        def fake_target_set(**kwargs):
            return kwargs

        fake_xt = SimpleNamespace(VaryList=fake_vary_list, TargetSet=fake_target_set)
        monkeypatch.setattr(phase_mod, "xt", fake_xt)

        line = FakeLine()
        phase_mod.match_extraction_tunes(line)

        assert line.vars["extr_bump_knob"] == pytest.approx(0.88)
        assert captured["solve"] is False
        assert captured["vary"] == [
            {"names": ["kqf", "kqd"], "step": 1e-7},
            {"names": ["qph_setvalue", "qpv_setvalue"], "step": 1e-4},
        ]
        assert captured["targets"][0] == {"qx": 26.666666666, "qy": 26.58, "tol": 1e-5}
        assert captured["targets"][1]["dqx"] == pytest.approx(-26.0)
        assert captured["targets"][1]["dqy"] == pytest.approx(0.47 * 26.0)
        assert captured["targets"][1]["tol"] == pytest.approx(1e-3)
        assert optimizer.solved is True
        assert line.vars["sps_on_extraction"] == pytest.approx(1.0)


class TestGenerateParticles:
    """Tests for matched particle generation."""

    def test_generate_particles_builds_4d_particles_with_expected_inputs(self, monkeypatch):
        """Particle generation must forward Gaussian samples and emittances to the line."""
        gaussian_samples = [
            (np.array([1.0, 2.0, 3.0]), np.array([-1.0, -2.0, -3.0])),
            (np.array([4.0, 5.0, 6.0]), np.array([-4.0, -5.0, -6.0])),
        ]
        captured = {}
        sentinel_particles = object()

        def fake_generate_2d_gaussian(num_particles):
            assert num_particles == 3
            return gaussian_samples.pop(0)

        class FakeLine:
            def build_particles(self, **kwargs):
                captured.update(kwargs)
                return sentinel_particles

        monkeypatch.setattr(
            phase_mod,
            "xp",
            SimpleNamespace(generate_2D_gaussian=fake_generate_2d_gaussian),
        )
        monkeypatch.setattr(phase_mod.np.random, "rand", lambda n: np.array([0.0, 0.25, 0.5]))

        result = phase_mod.generate_particles(FakeLine(), tw=object(), n_part=3)

        assert result is sentinel_particles
        assert captured["method"] == "4d"
        assert captured["zeta"] == pytest.approx(0.0)
        np.testing.assert_allclose(captured["delta"], np.array([0.0, 0.25, 0.5]) * phase_mod.DPP)
        np.testing.assert_array_equal(captured["x_norm"], np.array([1.0, 2.0, 3.0]))
        np.testing.assert_array_equal(captured["px_norm"], np.array([-1.0, -2.0, -3.0]))
        np.testing.assert_array_equal(captured["y_norm"], np.array([4.0, 5.0, 6.0]))
        np.testing.assert_array_equal(captured["py_norm"], np.array([-4.0, -5.0, -6.0]))
        assert captured["nemitt_x"] == pytest.approx(phase_mod.N_EX)
        assert captured["nemitt_y"] == pytest.approx(phase_mod.N_EY)


class TestTrackAndCollect:
    """Tests for per-turn coordinate collection."""

    def test_track_and_collect_rebuilds_cpu_tracker_and_stores_turnwise_copies(self, monkeypatch):
        """Tracking must rebuild the tracker on CPU and copy coordinates each turn."""
        captured = {"discarded": False, "contexts": [], "turns": []}
        fake_context = object()

        class FakeLine:
            def discard_tracker(self):
                captured["discarded"] = True

            def build_tracker(self, _context=None):
                captured["contexts"].append(_context)

            def track(self, particles, num_turns):
                captured["turns"].append(num_turns)
                particles.x += 1.0
                particles.px -= 1.0

        monkeypatch.setattr(phase_mod, "xo", SimpleNamespace(ContextCpu=lambda: fake_context))

        particles = FakeParticles([0.0], px=[0.0])
        positions_x, momenta_px = phase_mod.track_and_collect(FakeLine(), particles, num_turns=3)

        assert captured["discarded"] is True
        assert captured["contexts"] == [fake_context]
        assert captured["turns"] == [1, 1, 1]
        assert [coords.item() for coords in positions_x] == [1.0, 2.0, 3.0]
        assert [coords.item() for coords in momenta_px] == [-1.0, -2.0, -3.0]

        positions_x[0][0] = 99.0
        momenta_px[0][0] = 99.0
        assert particles.x[0] == pytest.approx(3.0)
        assert particles.px[0] == pytest.approx(-3.0)


class TestBuildAnimation:
    """Tests for matplotlib animation creation."""

    def test_build_animation_returns_funcanimation_and_update_runs(self):
        """Animation creation must succeed on a non-interactive backend."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.animation as mpl_animation
        import matplotlib.pyplot as plt

        anim = phase_mod.build_animation(
            positions_x=[np.array([0.02, 0.04])],
            momenta_px=[np.array([-1e-3, -5e-4])],
            num_turns=1,
        )

        assert isinstance(anim, mpl_animation.FuncAnimation)

        artists = anim._func(0)
        ax = anim._fig.axes[0]
        assert isinstance(artists, tuple)
        assert len(artists) == 1
        assert ax.get_xlabel() == "X Position (m)"
        assert ax.get_ylabel() == "Px (rad)"
        assert len(ax.lines) == 1

        anim._draw_was_started = True
        plt.close(anim._fig)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "--color=yes", "-v"]))