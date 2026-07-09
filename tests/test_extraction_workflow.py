"""Tests for crystal_extraction.extraction workflow helpers."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("xcoll")
pytest.importorskip("xtrack")

import crystal_extraction.extraction as extraction  # noqa: E402


class TestBuildSimulation:
    """Tests for build_simulation setup flow."""

    def test_build_simulation_loads_line_and_creates_particles(self, monkeypatch, tmp_path):
        captured = {}

        class FakeLine:
            def __init__(self):
                self.particle_ref = None

            def twiss(self, method):
                captured["twiss_method"] = method
                return "twiss_result"

        fake_line = FakeLine()

        class FakeBeamArgs:
            def __init__(self, nparticles):
                self.nparticles = nparticles
                self.emit_rms = 1.2e-9
                self.dpp_fullwidth = 2e-4

        class FakeSteinArgs:
            def __init__(self, approach_side, nsigma_emit):
                self.approach_side = approach_side
                self.nsigma_emit = nsigma_emit
                self.crystal_loc = "crystal"
                self.septum_loc = "septum"

        class FakeCrystal:
            def __init__(self):
                self.bending_angle = 0.0

            def copy(self):
                return FakeCrystal()

        class FakeSteinbach:
            def __init__(self, line, crystal, beam_args, stein_args):
                self.line = line
                self.crystal = crystal
                self.beam_args = beam_args
                self.stein_args = stein_args
                self.created = False

            def create_xsuite_particles(self):
                self.created = True

        def fake_from_json(path):
            captured["line_json"] = path
            return fake_line

        def fake_particles(**kwargs):
            captured["particle_ref_kwargs"] = kwargs
            return "particle_ref"

        monkeypatch.setattr(
            extraction,
            "xt",
            SimpleNamespace(
                Line=SimpleNamespace(from_json=fake_from_json),
                Particles=fake_particles,
                PROTON_MASS_EV=938e6,
            ),
        )
        monkeypatch.setattr(extraction, "BeamArgs", FakeBeamArgs)
        monkeypatch.setattr(extraction, "SteinArgs", FakeSteinArgs)
        monkeypatch.setattr(extraction, "Steinbach", FakeSteinbach)
        monkeypatch.setattr(extraction, "sps_crystal", FakeCrystal())

        line_path = tmp_path / "line.json"
        ste, tw0 = extraction.build_simulation(line_path=line_path, n_particles=321)

        assert captured["line_json"] == str(line_path)
        assert captured["particle_ref_kwargs"]["p0c"] == pytest.approx(400e9)
        assert captured["twiss_method"] == "4d"
        assert fake_line.particle_ref == "particle_ref"
        assert tw0 == "twiss_result"
        assert ste.created is True
        assert ste.beam_args.nparticles == 321


class TestComputeRmatrices:
    """Tests for R-matrix construction."""

    def test_compute_rmatrices_returns_expected_matrices(self):
        tw_crystal_df = pd.DataFrame(
            {"dpp_stopband_rms": [1.5e-3]}, index=["crystal"]
        )

        class FakeTwiss:
            def __init__(self, matrix_to_return):
                self._matrix_to_return = matrix_to_return

            def get_R_matrix(self, start, end):
                assert start in ("crystal", "septum")
                assert end in ("crystal", "septum")
                return self._matrix_to_return

        tw_first = FakeTwiss("to_mst")
        tw_second = FakeTwiss("to_cry")

        class FakeLine:
            def __init__(self):
                self.calls = 0

            def cycle(self, name_first_element):
                self.calls += 1
                if self.calls == 1:
                    return SimpleNamespace(
                        twiss=lambda method: tw_first,
                        cycle=lambda name_first_element: SimpleNamespace(
                            twiss=lambda method: tw_second
                        ),
                    )
                return SimpleNamespace(twiss=lambda method: tw_second)

        ste = SimpleNamespace(
            line=FakeLine(),
            stein_args=SimpleNamespace(crystal_loc="crystal", septum_loc="septum"),
            parse_xsuite_twiss=lambda tw: tw_crystal_df,
            crystal=SimpleNamespace(tilt=999.0),
        )

        r1, r2, twc = extraction.compute_rmatrices(ste, tw0="unused")

        assert r1 == "to_mst"
        assert r2 == "to_cry"
        assert twc["dpp_stopband_rms"] == pytest.approx(1.5e-3)
        assert ste.crystal.tilt == 0


class TestRunTracking:
    """Tests for matrix-tracking loop."""

    def test_run_tracking_marks_out_and_lost_and_increments_delta(self, monkeypatch):
        track_calls = []

        def fake_track_with_mat(particles, mat):
            track_calls.append(mat)

        monkeypatch.setattr(extraction.utils, "track_with_mat", fake_track_with_mat)

        class FakeParticles:
            def __init__(self):
                self.x = np.array([-0.011, -0.007, 0.001], dtype=float)
                self.state = np.array([1, 1, 1], dtype=int)
                self.delta = np.zeros(3, dtype=float)

            def copy(self):
                cp = FakeParticles()
                cp.x = self.x.copy()
                cp.state = self.state.copy()
                cp.delta = self.delta.copy()
                return cp

        crystal_calls = {"n": 0}

        def fake_crystal_track(particles):
            crystal_calls["n"] += 1

        ste = SimpleNamespace(
            stein_args=SimpleNamespace(approach_side=1, nstopbands_dpp=1),
            beam_args=SimpleNamespace(dpp_fullwidth=1.0, nparticles=3),
            xsuite_particles=FakeParticles(),
            crystal=SimpleNamespace(track=fake_crystal_track),
        )

        tw_crystal = {"dpp_stopband_rms": extraction.SPS_T_REV}
        particles, particles_init = extraction.run_tracking(
            ste,
            rmat_to_mst="to_mst",
            rmat_to_cry="to_cry",
            tw_crystal=tw_crystal,
            extraction_time=1.0,
        )

        # turns = int(2 * SPS_T_REV / SPS_T_REV) = 2
        assert crystal_calls["n"] == 2
        assert len(track_calls) == 6  # 2 pre-loop + 2 turns * 2 matrices
        assert np.allclose(particles.delta, np.full(3, 2 * extraction.SPS_T_REV))

        # x < -0.01 => out; -0.01 <= x <= -0.005 => lost at blade
        assert particles.state[0] == -42
        assert particles.state[1] == -1
        assert particles.state[2] == 1
        assert np.array_equal(particles_init.state, np.array([1, 1, 1]))


class TestReportEfficiency:
    """Tests for textual summary output."""

    def test_report_efficiency_prints_portion_when_particles_removed(self, capsys):
        particles = SimpleNamespace(state=np.array([-42, -1, 1, 1], dtype=int))
        extraction.report_efficiency(particles, n_total=4)
        out = capsys.readouterr().out
        assert "Extraction efficiency summary" in out
        assert "Portion (nout/nlost)" in out

    def test_report_efficiency_omits_portion_when_no_losses(self, capsys):
        particles = SimpleNamespace(state=np.array([1, 1, 1], dtype=int))
        extraction.report_efficiency(particles, n_total=3)
        out = capsys.readouterr().out
        assert "Extraction efficiency summary" in out
        assert "Portion (nout/nlost)" not in out


class TestSavePlots:
    """Tests for output plot generation."""

    def test_save_plots_writes_expected_files(self, tmp_path):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig1, _ = plt.subplots()
        fig2, _ = plt.subplots()

        ste = SimpleNamespace(plot=lambda tw0: ([fig1, fig2], None))

        particles = SimpleNamespace(
            x=np.array([-0.011, -0.007, 0.001]),
            px=np.array([0.0, 0.0, 0.0]),
            state=np.array([-42, -1, 1], dtype=int),
        )
        particles_init = SimpleNamespace(
            x=np.array([-0.002, -0.001, 0.0]),
            px=np.array([0.0, 0.0, 0.0]),
            state=np.array([1, 1, 1], dtype=int),
        )

        extraction.save_plots(
            ste=ste,
            particles=particles,
            particles_init_mst=particles_init,
            tw0="twiss",
            output_dir=tmp_path,
        )

        assert (tmp_path / "steinbach_0.png").exists()
        assert (tmp_path / "steinbach_1.png").exists()
        assert (tmp_path / "phase_space_mst.png").exists()


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "--color=yes", "-v"]))
