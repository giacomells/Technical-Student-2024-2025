import pickle
import typing as t

import matplotlib.pyplot as plt
import numpy as np
import requests
import xobjects as xo
import xpart as xp
import xtrack as xt
from cpymad.madx import Madx

p = 400.0  # beam momentum (GeV/c)
momentum = p  # beam momentum (GeV/c)
Brho = p * 3.3356  # beam rigidity ???

N_EX = 10e-6
N_EY = 5e-6
DPP = 1e-4

mad = Madx(stdout=True)

mad.input(
    requests.get(
        "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/SPS_LS2_2020-05-26.seq"
    ).text
)
mad.input(
    requests.get(
        "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/strengths/ft_q26_extr.str"
    ).text
)
mad.input(
    requests.get(
        "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/toolkit/macro.madx"
    ).text
)
mad.input(
    requests.get(
        "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/aperture/aperturedb_classes.madx"
    ).text
)

mad.input(
    requests.get(
        "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/aperture/aperturedb_elements.madx"
    ).text
)

mad.command.beam(particle="PROTON", pc=p, charge=1)
mad.input("BRHO = BEAM->PC * 3.3356;")
mad.use(sequence="sps")

mad.input("""
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

""")

mad.use(sequence="sps")

line = xt.Line.from_madx_sequence(mad.sequence["sps"], deferred_expressions=True,
                                  allow_thick=True)

line.particle_ref = xt.Particles(mass0=xp.PROTON_MASS_EV,
                                gamma0=mad.sequence["sps"].beam.gamma)
line.twiss_default['method'] = '4d'

# Prepare rf
line.vv['v200'] = 0.
line.vv['freq200'] = 200e6
line.vv['lag200'] = 180.
line.element_refs['actcse.31632'].voltage = line.vars['v200']
line.element_refs['actcse.31632'].frequency = line.vars['freq200']
line.element_refs['actcse.31632'].lag = line.vars['lag200']

line.vv['v200'] = 10e6
tw6d = line.twiss(method='6d')
assert tw6d.qs > 1e-3

line.vv['v200'] = 0.

line.to_json("sps_for_sx.json")