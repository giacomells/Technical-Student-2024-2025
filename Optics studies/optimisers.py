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


def horizontal_bumpLSS2(line, x_target, px_target):
    """
    Adjusts the horizontal bump at tpst.21760_entry using four elements.
    
    Parameters:
        line (xt.Line): The accelerator beamline.
        x_target (float): Desired horizontal position at tpst.21760_entry.
        px_target (float): Desired horizontal angle at tpst.21760_entry.
    """
    opt = line.match(
        start='mpsh.21202', end='mplh.22195',  # Bump region
        betx=1, bety=1, x=0, px=0,  # Keep initial conditions unchanged
        vary=xt.VaryList(
            ['kmpsh21202', 'kmplh21431', 'kmplh21995', 'kmplh22195'],  # Selected correctors
            step=1e-10,
            limits=[-1e-3, 1e-3]  # Define kick limits to avoid excessive changes
        ),
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START),  # Ensure bump is closed
            xt.TargetSet(x=x_target, px=px_target, at='tpst.21760_entry'),  # Apply bump at target location
            #xt.TargetSet(x=x_target, px=px_target, at='ap.up.mst21774'),  # Apply bump at target location
            xt.TargetSet(x=0, px=0, at=xt.END)  # Ensure bump is closed
        ]
    )
    return opt




def horizontal_bumpLSS4(line):
    """
    Adjusts the horizontal bump at tpst.21760_entry using four elements.
    
    Parameters:
        line (xt.Line): The accelerator beamline.
        x_target (float): Desired horizontal position at tpst.21760_entry.
        px_target (float): Desired horizontal angle at tpst.21760_entry.
    """
    opt = line.match(
        start='mpsh.41402', end='mpsh.42198',  # Bump region
        betx=1, bety=1, x=0, px=0,  # Keep initial conditions unchanged
        vary=xt.VaryList(
            ['kmpsh41402', 'kmplh41658', 'kmplh41994', 'kmpsh42198'],  # Selected correctors
            step=1e-10,
            limits=[-1e-3, 1e-3]  # Define kick limits to avoid excessive changes
        ),
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START), 
            xt.TargetSet(x = TECA.jaw + TECA.width, px=TECA.tilt, at='TECA.entry'),  # Apply bump at target location
            xt.TargetSet(x=0, px=0, at=xt.END)  # Ensure bump is closed
        ]
    )
    return opt