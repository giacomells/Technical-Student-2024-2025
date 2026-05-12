import os

import pytest

# When this file is executed directly, enable the integration checks automatically.
# Plain pytest collection still keeps them opt-in unless the variable is exported.
if __name__ == "__main__":
    os.environ.setdefault("RUN_CERN_GITLAB_TESTS", "1")

RUN_CERN_GITLAB_TESTS = os.environ.get("RUN_CERN_GITLAB_TESTS") == "1"

requests = pytest.importorskip("requests")
xp = pytest.importorskip("xpart")
xt = pytest.importorskip("xtrack")
pytest.importorskip("xobjects")
Madx = pytest.importorskip("cpymad.madx").Madx


# CERN GitLab sources needed to recreate the SPS sequence used by the example.
SEQUENCE_URLS = [
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/SPS_LS2_2020-05-26.seq",
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/strengths/ft_q26_extr.str",
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/toolkit/macro.madx",
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/aperture/aperturedb_classes.madx",
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/aperture/aperturedb_elements.madx",
]

# Extra MAD-X definitions copied from the SPS example so the imported sequence
# includes the extraction bump, installed marker and chromaticity knobs.
EXTRA_MADX_INPUT = """
extr_bump_knob = 0;
kMPLH21431 := 4.9e-4 * extr_bump_knob;
kMPLH21995 := 2.503e-4 * extr_bump_knob;
kMPLH22195 := -3.5585e-4 * extr_bump_knob;
kMPNH21732 := 3.3309e-4 * extr_bump_knob;
kMPSH21202 := -7.6765e-5 * extr_bump_knob;

seqedit, sequence = sps;
install, element=ap.up.zs21633, class = marker, at=-1 * zs.21633->L/2, from = zs.21633;
endedit;

LSDA0 = -0.149628261;
LSDB0 = -0.145613183;
LSFA0 = 0.063256459;
LSFB0 = 0.121416689;
LSFC0 = 0.063256459;

logical.LSDAQPH = .011283;
logical.LSDBQPH = -.040346;
logical.LSFAQPH = .04135;
logical.LSFBQPH = .079565;
logical.LSFCQPH = .04135;

logical.LSDAQPV = -.11422;
logical.LSDBQPV = -.08606;
logical.LSFAQPV = .0097365;
logical.LSFBQPV = .016931;
logical.LSFCQPV = .0097365;

kLSDA := logical.LSDAQPH*QPH_setvalue + logical.LSDAQPV*QPV_setvalue + LSDA0;
kLSDB := logical.LSDBQPH*QPH_setvalue + logical.LSDBQPV*QPV_setvalue + LSDB0;
kLSFA := logical.LSFAQPH*QPH_setvalue + logical.LSFAQPV*QPV_setvalue + LSFA0;
kLSFB := logical.LSFBQPH*QPH_setvalue + logical.LSFBQPV*QPV_setvalue + LSFB0;
kLSFC := logical.LSFCQPH*QPH_setvalue + logical.LSFCQPV*QPV_setvalue + LSFC0;
"""


def _load_madx_from_cern_gitlab(momentum_gev_c=400.0):
    # Build a MAD-X instance from the same remote sequence and strengths files
    # used in the animation example.
    mad = Madx(stdout=False)

    for url in SEQUENCE_URLS:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        mad.input(response.text)

    mad.command.beam(particle="PROTON", pc=momentum_gev_c, charge=1)
    mad.input("BRHO = BEAM->PC * 3.3356;")
    mad.use(sequence="sps")
    mad.input(EXTRA_MADX_INPUT)
    mad.use(sequence="sps")
    return mad


def _build_sps_line_from_gitlab():
    # Convert the MAD-X sequence into an Xsuite line so Twiss and RF checks can
    # be performed directly in Python.
    mad = _load_madx_from_cern_gitlab()
    line = xt.Line.from_madx_sequence(
        mad.sequence["sps"],
        deferred_expressions=True,
        allow_thick=True,
    )
    line.particle_ref = xt.Particles(
        mass0=xp.PROTON_MASS_EV,
        gamma0=mad.sequence["sps"].beam.gamma,
    )
    line.twiss_default["method"] = "4d"
    return line


@pytest.mark.skipif(
    not RUN_CERN_GITLAB_TESTS,
    reason="Set RUN_CERN_GITLAB_TESTS=1 to run CERN GitLab sequence creation checks.",
)
def test_can_build_sps_sequence_from_cern_gitlab():
    # Smoke test: the sequence should load, convert, and expose the expected SPS
    # line length and reference elements.
    line = _build_sps_line_from_gitlab()

    assert line.particle_ref is not None
    assert line.get_length() > 6000
    assert "actcse.31632" in line.element_names
    assert "ap.up.zs21633" in line.element_names


@pytest.mark.skipif(
    not RUN_CERN_GITLAB_TESTS,
    reason="Set RUN_CERN_GITLAB_TESTS=1 to run CERN GitLab sequence creation checks.",
)
def test_rf_setup_supports_6d_twiss_after_sequence_creation():
    # Recreate the RF configuration from the example and verify that enabling the
    # cavity produces a valid 6D Twiss solution.
    line = _build_sps_line_from_gitlab()

    line.vv["v200"] = 0.0
    line.vv["freq200"] = 200e6
    line.vv["lag200"] = 180.0
    line.element_refs["actcse.31632"].voltage = line.vars["v200"]
    line.element_refs["actcse.31632"].frequency = line.vars["freq200"]
    line.element_refs["actcse.31632"].lag = line.vars["lag200"]

    line.vv["v200"] = 10e6
    twiss_6d = line.twiss(method="6d")

    assert twiss_6d.qs > 1e-3


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "--color=yes", "-v"]))