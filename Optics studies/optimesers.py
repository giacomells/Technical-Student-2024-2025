import numpy as np
import pandas as pd
import scipy.stats as stats
import typing as t

import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from shapely.geometry import Polygon

import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc

from elements import TECA

from cpymad.madx import Madx

def match_tunes(line, qx, qy):
    # Extraction tunes
    opt = line.match(solve=False,
                     method='4d',
        vary=[
            xt.VaryList(['kqf', 'kqd'], step=1e-7),   # Varying quadrupole focal strengths
            xt.VaryList(['qph_setvalue', 'qpv_setvalue'], step=1e-4),   # Varying phase values
        ],
        targets=[
            xt.TargetSet(qx=qx, qy=qy, tol=1e-5),# Desired target tunes
            #xt.TargetSet(dqx=-1 * qx, dqy=0.47 * qy, tol=1e-3),   # Desired target chromaticities
        ])
    return opt

def match_chromaticity(line, qx, qy):
    # Extraction tunes
    opt = line.match(solve=False,
                     method='4d',
        vary=[
            xt.VaryList(['klsfa', 'klsda', 'klsdb', 'klsfb'], step=1e-7),   # Varying setupoles strengths ,  klsfa, klsda, klsfb, klsdb????
        ],
        targets=[
            xt.TargetSet(dqx=-1 * int(qx), dqy = 0.47 * int(qy), tol=1e-3),   # Desired target chromaticities
        ])
    return opt



def adjust_bumps(line):
    opt =line.match(
        start='mq.30l8.b1', end='mq.23l8.b1',
        betx=1, bety=1, y=0, py=0, # <-- conditions at start
        vary=xt.VaryList(['acbv30.l8b1', 'acbv28.l8b1', 'acbv26.l8b1', 'acbv24.l8b1'],
                    step=1e-10, limits=[-1e-3, 1e-3]),
        targets = [
            xt.TargetSet(y=3e-3, py=0, at='mb.b28l8.b1'),
            xt.TargetSet(y=0, py=0, at=xt.END)
        ])
