"""
Build and serialise the SPS Q26 extraction line from CERN GitLab model files.

Use save_sps_json() from an application script to generate
database/sps_for_sx.json, which can then be loaded without a network
connection.
"""

from __future__ import annotations

from pathlib import Path

import requests
import xpart as xp
import xtrack as xt
from cpymad.madx import Madx

# ---------------------------------------------------------------------------
# Physics constants
# ---------------------------------------------------------------------------

p: float = 400.0       # beam momentum [GeV/c]
Brho: float = p * 3.3356  # magnetic rigidity [T·m]
N_EX: float = 10e-6   # normalised horizontal emittance [m·rad]
N_EY: float = 5e-6    # normalised vertical emittance   [m·rad]
DPP: float = 1e-4     # relative momentum spread

# ---------------------------------------------------------------------------
# CERN GitLab SPS model URLs (2021 tag)
# ---------------------------------------------------------------------------

SPS_MODEL_URLS: list = [
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/SPS_LS2_2020-05-26.seq",
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/strengths/ft_q26_extr.str",
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/toolkit/macro.madx",
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/aperture/aperturedb_classes.madx",
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/aperture/aperturedb_elements.madx",
]

# MAD-X input that installs the extraction bump knob, the ZS marker and the
# chromaticity knob variables used in the Xsuite deferred-expression system.
_EXTRA_MADX_INPUT: str = """
extr_bump_knob = 0;
kMPLH21431 := 4.9e-4 * extr_bump_knob;
kMPLH21995 := 2.503e-4 * extr_bump_knob;
kMPLH22195 := -3.5585e-4 * extr_bump_knob;
kMPNH21732 := 3.3309e-4 * extr_bump_knob;
kMPSH21202 := -7.6765e-5 * extr_bump_knob;

! Install marker at ZS
seqedit, sequence = sps;
install, element=ap.up.zs21633, class = marker, at=-1 * zs.21633->L/2, from = zs.21633;
endedit;

! Build chromaticity knobs
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


def build_sps_madx(momentum_gev_c: float = 400.0) -> Madx:
    """Download SPS model files from CERN GitLab and return a MAD-X instance.

    Each URL in :data:`SPS_MODEL_URLS` is fetched over HTTPS and fed to
    MAD-X.  The beam is set up at *momentum_gev_c* GeV/c and the extraction
    knob / chromaticity variables are defined.

    Parameters
    ----------
    momentum_gev_c : float
        Proton beam momentum [GeV/c].  Default is 400 GeV/c.

    Returns
    -------
    cpymad.madx.Madx
        Fully configured MAD-X instance with ``sps`` as the active sequence.
    """
    mad = Madx(stdout=False)
    for url in SPS_MODEL_URLS:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        mad.input(response.text)

    mad.command.beam(particle="PROTON", pc=momentum_gev_c, charge=1)
    mad.input("BRHO = BEAM->PC * 3.3356;")
    mad.use(sequence="sps")
    mad.input(_EXTRA_MADX_INPUT)
    mad.use(sequence="sps")
    return mad


def save_sps_json(
    output_path: Path | None = None,
    momentum_gev_c: float = 400.0,
) -> Path:
    """Build the SPS extraction line and write it to a JSON file.

    Calls :func:`build_sps_madx`, converts the MAD-X sequence to an Xsuite
    line, configures the RF cavity so that a 6-D Twiss can be verified, and
    serialises the result.

    Parameters
    ----------
    output_path : Path or None
        Destination file.  Defaults to ``database/sps_for_sx.json`` at the
        repository root.
    momentum_gev_c : float
        Proton beam momentum [GeV/c].

    Returns
    -------
    Path
        Absolute path of the written JSON file.
    """
    if output_path is None:
        output_path = Path(__file__).resolve().parent.parent / "database" / "sps_for_sx.json"
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    mad = build_sps_madx(momentum_gev_c)

    line = xt.Line.from_madx_sequence(
        mad.sequence["sps"], deferred_expressions=True, allow_thick=True
    )
    line.particle_ref = xt.Particles(
        mass0=xp.PROTON_MASS_EV,
        gamma0=mad.sequence["sps"].beam.gamma,
    )
    line.twiss_default["method"] = "4d"

    # Verify RF closes a 6-D solution, then switch cavity off before saving.
    line.vv["v200"] = 0.0
    line.vv["freq200"] = 200e6
    line.vv["lag200"] = 180.0
    line.element_refs["actcse.31632"].voltage = line.vars["v200"]
    line.element_refs["actcse.31632"].frequency = line.vars["freq200"]
    line.element_refs["actcse.31632"].lag = line.vars["lag200"]

    line.vv["v200"] = 10e6
    tw6d = line.twiss(method="6d")
    assert tw6d.qs > 1e-3, "6-D Twiss failed: synchrotron tune is too small."
    line.vv["v200"] = 0.0

    line.to_json(str(output_path))
    return output_path