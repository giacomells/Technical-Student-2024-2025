import numpy as np
import pandas as pd
import scipy.stats as stats
import typing as t

import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import minimize

import xobjects as xo
import xpart as xp
import xtrack as xt
import xcoll as xc

from elements import TECA
from elements import open_blocking_apertures, remove_ZS_apertures, remove_inner_sideLimits_closeTECA, save_df_Limit_elements_features, deltaP_P

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

def match_tunesQ22(line):
    # Extraction tunes
    opt = line.match(solve=False,
                     method='4d',
        vary=[
            xt.VaryList(['kqf', 'kqd'], step=1e-7),   # Varying quadrupole focal strengths
            xt.VaryList(['qph_setvalue', 'qpv_setvalue'], step=1e-4),   # Varying phase values
        ],
        targets=[
            xt.TargetSet(qx=22.13, qy=22.18, tol=1e-5),# Desired target tunes
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

def match_chromaticityQ22(line, qx, qy):
    # Extraction tunes
    opt = line.match(solve=False,
                     method='4d',
        vary=[
            xt.VaryList(['klsfa', 'klsda', 'klsdb', 'klsfb'], step=1e-7),   # Varying setupoles strengths ,  klsfa, klsda, klsfb, klsdb????
        ],
        targets=[
            xt.TargetSet(dqx= 0.288 * int(qx), dqy = 0.1944 * int(qy), tol=1e-3),   # Desired target chromaticities
        ])
    return opt


def horizontal_bumpLSS2(line, x_target = 0.039): # I REMOVED THE PX_TARGET
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
            ['kmpsh21202', 'kmplh21431', 'kmpnh21732', 'kmplh21995', 'kmplh22195'],  # Selected correctors. added kmpnh21732 in date 04.04.2025 from Aleksandr example
            step=1e-10,
            limits=[-4e-4, 1e-3]  # Define kick limits to avoid excessive changes
        ),
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START),  # Ensure bump is closed
            #xt.TargetSet(x=x_target, px=px_target, at='zs21633.entry.p1mm'),  # Apply bump at target location
            #xt.TargetSet(x=x_target, px=px_target, at='tpst.21760_entry'),  # Apply bump at target location
            # xt.TargetSet(x=x_target, at='tpst.21760_entry'),  # Apply bump at target location
            # xt.TargetSet(px=px_target, at='tpst.21760_entry'),  # Apply bump at target location
            xt.Target(x=x_target, at='tpst.21760_entry'),  # Apply bump at target location
            xt.Target(px= xt.GreaterThan(0), at='tpst.21760_entry'),  # Apply bump at target location
            # #xt.TargetSet(x=x_target, px=px_target, at='ap.up.mst21774'),  # Apply bump at target location
            xt.TargetSet(x=0, px=0, at=xt.END)  # Ensure bump is closed
            #xt.TargetSet(x=0, px=0, at="qf.22010")  # Ensure bump is closed
            
        ]
    )
    return opt

# def horizontal_bumpLSS2(line, x_target, px_target):
#     """
#     Adjusts the horizontal bump at tpst.21760_entry using four elements.
    
#     Parameters:
#         line (xt.Line): The accelerator beamline.
#         x_target (float): Desired horizontal position at tpst.21760_entry.
#         px_target (float): Desired horizontal angle at tpst.21760_entry.
#     """
#     opt = line.match(
#         start='mpsh.21202', end='mplh.22195',  # Bump region
#         betx=1, bety=1, x=0, px=0,  # Keep initial conditions unchanged
#         vary=xt.VaryList(
#             ['kmpsh21202', 'kmplh21431', 'kmpnh21732', 'kmplh21995', 'kmplh22195'],  # Selected correctors. added kmpnh21732 in date 04.04.2025 from Aleksandr example
#             step=1e-10,
#             limits=[-1e-3, 1e-3]  # Define kick limits to avoid excessive changes
#         ),
#         targets=[
#             xt.TargetSet(x=0, px=0, at=xt.START),  # Ensure bump is closed
#             #xt.TargetSet(x=x_target, px=px_target, at='zs21633.entry.p1mm'),  # Apply bump at target location
#             xt.TargetSet(x=x_target, px=px_target, at='tpst.21760_entry'),  # Apply bump at target location
#             #xt.TargetSet(x=x_target, px=px_target, at='ap.up.mst21774'),  # Apply bump at target location
#             xt.TargetSet(x=0, px=0, at=xt.END)  # Ensure bump is closed
#             #xt.TargetSet(x=0, px=0, at="qf.22010")  # Ensure bump is closed
            
#         ]
#     )
#     return opt



def horizontal_bumpLSS2atZS(line, x_target, px_target):
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
            ['kmpsh21202', 'kmplh21431', 'kmpnh21732', 'kmplh21995', 'kmplh22195'],  # Selected correctors. added kmpnh21732 in date 04.04.2025 from Aleksandr example
            step=1e-10,
            limits=[-1e-3, 1e-3]  # Define kick limits to avoid excessive changes
        ),
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START),  # Ensure bump is closed
            xt.TargetSet(x=x_target, px=px_target, at='zs21633.entry.p1mm'),  # Apply bump at target location
            #xt.TargetSet(x=x_target, px=px_target, at='tpst.21760_entry'),  # Apply bump at target location
            #xt.TargetSet(x=x_target, px=px_target, at='ap.up.mst21774'),  # Apply bump at target location
            xt.TargetSet(x=0, px=0, at=xt.END)  # Ensure bump is closed
            #xt.TargetSet(x=0, px=0, at="qf.22010")  # Ensure bump is closed
            
        ]
    )
    return opt

def ensure_closed_orbit(line):
    """
    Ensures that the bump is closed at the start and end of the line.
    
    Parameters:
        line (xt.Line): The accelerator beamline.
    """
    opt = line.match(
        start='begi.10010', end='end.10010',  # Bump region
        betx=1, bety=1, x=0, px=0,  # Keep initial conditions unchanged
        vary=xt.VaryList(
            #['kmpsh21202', 'kmplh21431', 'kmpnh21732', 'kmplh21995', 'kmplh22195', 'kmpsh41402', 'kmplh41658', 'kmplh41994', 'kmpsh42198'],  # Selected correctors. added kmpnh21732 in date 04.04.2025 from Aleksandr example
            ['kmpsh41402', 'kmplh41658', 'kmplh41994', 'kmpsh42198'],  # Selected correctors. added kmpnh21732 in date 04.04.2025 from Aleksandr example
            
            step=1e-10,
            limits=[-1e-3, 1e-3]  # Define kick limits to avoid excessive changes
        ),
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START),  # Ensure bump is closed
            xt.TargetSet(x=0, px=0, at= 'drift_779'),  # Ensure bump is closed
            xt.TargetSet(x=0, px=0, at=xt.END)  # Ensure bump is closed
        ]
    )
    return opt


TECA39MM = xc.EverestCrystal(
    length=2e-3, 
    material=xc.materials.SiliconCrystal, 
    bending_angle = - 174e-6 ,              # THIS CRYSTAL IS CHANNELING TOWARDS THE INSIDE OF THE RING!!!
    side="-",
    lattice="strip",
    jaw = - 39e-3,  # CLOSEST POSITION REACHABLE TO THE CENTER OF THE BEAM PIPE
    #jaw = - 52e-3,  #original setting
    tilt = - 1.23e-3,
    width = 0.8e-3,
    height = 50e-3
    )

def horizontal_bumpLSS4_39MM(line, x_target = (TECA39MM.jaw - TECA39MM.width), px_target = TECA39MM.tilt):
    """
    Adjusts the horizontal bump at tpst.21760_entry using four elements.
    
    Parameters:
        line (xt.Line): The accelerator beamline.
        x_target (float): Desired horizontal position at tpst.21760_entry.
        px_target (float): Desired horizontal angle at tpst.21760_entry.
    """
    opt = line.match(
        start='mpsh.41402', end='mpsh.42198',  # Bump region
        #start='begi.10010', end='end.10010',
        betx=1, bety=1, x=0, px=0,  # Keep initial conditions unchanged
        vary=xt.VaryList(
            ['kmpsh41402', 'kmplh41658', 'kmplh41994', 'kmpsh42198'],  # Selected correctors
            #['kmplh41658', 'kmplh41994', 'kmpsh42198'],  # Selected correctors with one corrector removed (kmpsh41402)
            step=1e-10,
            limits=[-1e-3, 1e-3]  # Define kick limits to avoid excessive changes
        ),
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START), 
            xt.TargetSet(x = x_target, px = px_target, at='TECA.entry'),  # Apply bump at target location
            #xt.TargetSet(x=0, px=0, at= 'drift_779'),  # Ensure bump is closed
            xt.TargetSet(x=0, px=0, at=xt.END)  # Ensure bump is closed
        ]
    )
    return opt

def horizontal_bumpLSS4NEW(line, x_target = (TECA.jaw - TECA.width), px_target = TECA.tilt + 0 * 10e-6):
    """
    Adjusts the horizontal bump at tpst.21760_entry using four elements.
    
    Parameters:
        line (xt.Line): The accelerator beamline.
        x_target (float): Desired horizontal position at tpst.21760_entry.
        px_target (float): Desired horizontal angle at tpst.21760_entry.
    """
    opt = line.match(
        start='mpsh.41402', end='mpsh.42198',  # Bump region
        #start='begi.10010', end='end.10010',
        betx=1, bety=1, x=0, px=0,  # Keep initial conditions unchanged
        vary=xt.VaryList(
            ['kmpsh41402', 'kmplh41658', 'kmplh41994', 'kmpsh42198'],  # Selected correctors
            #['kmplh41658', 'kmplh41994', 'kmpsh42198'],  # Selected correctors with one corrector removed (kmpsh41402)
            step=1e-10,
            limits=[-1e-3, 1e-3]  # Define kick limits to avoid excessive changes
        ),
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START), 
            xt.TargetSet(x = x_target, px = px_target, at='TECA.entry'),  # Apply bump at target location
            #xt.TargetSet(x=0, px=0, at= 'drift_779'),  # Ensure bump is closed
            xt.TargetSet(x=0, px=0, at=xt.END)  # Ensure bump is closed
        ]
    )
    return opt



def horizontal_bumpLSS4_Francesco(line):
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
            xt.TargetSet(x = TECA.jaw + TECA.width * 0.7, px=TECA.tilt, at='teca'),  # Apply bump at target location
            xt.TargetSet(x=0, px=0, at=xt.END)  # Ensure bump is closed
        ]
    )
    return opt


def optimize_bumps(line, x_target, px_target):
    opt = horizontal_bumpLSS4(line)
    opt = horizontal_bumpLSS2(line, x_target, px_target)
    
    return opt




def rematch_optics(line, tune = 24.39, change_aperture = True):
    
    optTune = match_tunes(line, tune, tune - 0.02)
    optChromaticity = match_chromaticity(line, tune, tune - 0.02)

    optTune.step(10)
    optTune.target_status()
    optTune.vary_status()

    # CHanging the chromaticity afterwards
    optChromaticity.step(10)

    optChromaticity.target_status()
    optChromaticity.vary_status() 

    line.discard_tracker()
    
        
    return optTune, optChromaticity

def rematch_opticsQ22(line, tune = 22.13, change_aperture = True):
    
    optTune = match_tunes(line, tune, tune + 0.05)
    optChromaticity = match_chromaticity(line, tune, tune + 0.05)

    optTune.step(10)
    optTune.target_status()
    optTune.vary_status()

    # CHanging the chromaticity afterwards
    optChromaticity.step(10)

    optChromaticity.target_status()
    optChromaticity.vary_status() 

    line.discard_tracker()
    
        
    return optTune, optChromaticity








#### ORTHOGONAL KNOBS FOR THE BUMPS

def set_x_knobLSS4(line):
    """
    Create a knob that moves the beam by +1 mm in x at TECA.entry.
    """
    opt = line.match(
        #knob_name="x_teca_knob",
        #run=True,  # Run the matching now
        method="4d",
        vary=[
            xt.VaryList(
                [
                    "kmplh41658",
                    "kmplh41994",
                    "kmpsh41402",
                    "kmpsh42198",
                ],
                step=1e-6,
            )
        ], 
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START),  # Ensure bump is closed
            xt.TargetSet(["x", "px"], 0.0, at="qf.42210"),
            xt.Target("x",  1e-3, at="TECA.entry"),
            xt.Target(px=0, at="TECA.entry"),            #ENSURING THIS KNOB IS ORTHOGONAL TO THE px TECA KNOB
            xt.TargetSet(["x", "px"], 0.0, at="qd.41310"),
            xt.TargetSet(x=0, px=0, at=xt.END),  # Ensure bump is closed

    
        ],
    )
    return opt

def set_x_knobLSS2(line):
    """
    Create a knob that moves the beam by +1 mm in x at TECA.entry.
    """
    opt = line.match(
        #knob_name="x_teca_knob",
        #run=True,  # Run the matching now
        method="4d",
        vary=[
            xt.VaryList(
                [
                'kmpsh21202', 'kmplh21431', 'kmpnh21732', 'kmplh21995', 'kmplh22195',  # Selected correctors. added kmpnh21732 in date 04.04.2025 from Aleksandr example
                ],
                step=1e-6,
            )
        ], 
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START),  # Ensure bump is closed
            #xt.TargetSet(["x", "px"], 0.0, at="qf.42210"),
            xt.Target("x",  1e-3, at='zs21633.entry.p1mm'),
            xt.Target(px=0, at='zs21633.entry.p1mm'),           
            xt.Target('x', xt.LessThan(2e-3), at='mplh.21995'), # <-- inequality
            #xt.TargetSet(["x", "px"], 0.0, at="qd.41310"),
            xt.TargetSet(x=0, px=0, at=xt.END),  # Ensure bump is closed

    
        ],
    )
    return opt

def normalized_x_knobLSS4(line):
    """
    Create a knob that moves the beam by +1 mm in x at TECA.entry.
    """
    opt = line.match(
        #knob_name="x_teca_knob",
        #run=True,  # Run the matching now
        method="4d",
        vary=[
            xt.VaryList(
                [
                    "kmplh41658",
                    "kmplh41994",
                    "kmpsh41402",
                    "kmpsh42198",
                ],
                step=1e-6,
            )
        ], 
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START),  # Ensure bump is closed
            xt.TargetSet(["x", "px"], 0.0, at="qf.42210"),
            xt.Target("x", TECA.jaw - TECA.width, at="TECA.entry"),
            xt.Target(px=0, at="TECA.entry"),            #ENSURING THIS KNOB IS ORTHOGONAL TO THE px TECA KNOB
            xt.TargetSet(["x", "px"], 0.0, at="qd.41310"),
            xt.TargetSet(x=0, px=0, at=xt.END),  # Ensure bump is closed

    
        ],
    )
    return opt



def set_px_knobLSS4(line):
    """
    Create a knob that changes beam angle by +1 µrad (1e-6 rad) at TECA.entry.
    """
    opt = line.match(
        #knob_name="px_teca_knob",
        #run=True,
        method="4d",
        vary=[
            xt.VaryList(
                [
                    "kmplh41658",
                    "kmplh41994",
                    "kmpsh41402",
                    "kmpsh42198",
                ],
                step=1e-9,
            )
        ],
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START),  # Ensure bump is closed
            xt.TargetSet(["x", "px"], 0.0, at="qf.42210"),
            xt.Target("px",  1e-6, at="zs21633.entry.p1mm"),
            xt.Target(x = 0, at = "zs21633.entry.p1mm"),            #ENSURING THIS KNOB IS ORTHOGONAL TO THE px TECA KNOB
            xt.TargetSet(["x", "px"], 0.0, at="qd.41310"),
            xt.TargetSet(x=0, px=0, at=xt.END),  # Ensure bump is closed
        ],
    )
    return opt

def set_px_knobLSS4(line):
    """
    Create a knob that changes beam angle by +1 µrad (1e-6 rad) at TECA.entry.
    """
    opt = line.match(
        #knob_name="px_teca_knob",
        #run=True,
        method="4d",
        vary=[
            xt.VaryList(
                [
                    "kmplh41658",
                    "kmplh41994",
                    "kmpsh41402",
                    "kmpsh42198",
                ],
                step=1e-9,
            )
        ],
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START),  # Ensure bump is closed
            xt.TargetSet(["x", "px"], 0.0, at="qf.42210"),
            xt.Target("px",  1e-6, at="TECA.entry"),
            xt.Target(x = 0, at = "TECA.entry"),            #ENSURING THIS KNOB IS ORTHOGONAL TO THE px TECA KNOB
            xt.TargetSet(["x", "px"], 0.0, at="qd.41310"),
            xt.TargetSet(x=0, px=0, at=xt.END),  # Ensure bump is closed
        ],
    )
    return opt

def set_px_knobLSS2(line):
    """
    Create a knob that changes beam angle by +1 µrad (1e-6 rad) at zs21633.entry.p1mm.
    """
    opt = line.match(
        #knob_name="px_teca_knob",
        #run=True,
        method="4d",
        vary=[
            xt.VaryList(
                [
                    'kmpsh21202', 'kmplh21431', 'kmpnh21732', 'kmplh21995', 'kmplh22195',  # Selected correctors. added kmpnh21732 in date 04.04.2025 from Aleksandr example
                ],
                step=1e-9,
            )
        ],
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START),  # Ensure bump is closed
            #xt.TargetSet(["x", "px"], 0.0, at="qf.42210"),
            xt.Target("px", 1e-6, at='tpst.21760_entry'),
            xt.Target(x=0, at='tpst.21760_entry'),  # ENSURING THIS KNOB IS ORTHOGONAL TO THE px TECA KNOB
            #xt.TargetSet(["x", "px"], 0.0, at="qd.41310"),
            xt.TargetSet(x=0, px=0, at=xt.END),  # Ensure bump is closed
        ],
    )
    return opt




def normalized_px_knobLSS4(line):
    """
    Create a knob that changes beam angle by +1 µrad (1e-6 rad) at TECA.entry.
    """
    opt = line.match(
        #knob_name="px_teca_knob",
        #run=True,
        method="4d",
        vary=[
            xt.VaryList(
                [
                    "kmplh41658",
                    "kmplh41994",
                    "kmpsh41402",
                    "kmpsh42198",
                ],
                step=1e-9,
            )
        ],
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START),  # Ensure bump is closed
            xt.TargetSet(["x", "px"], 0.0, at="qf.42210"),
            xt.Target("px",  TECA.tilt, at="TECA.entry"),
            xt.Target(x = 0, at = "TECA.entry"),            #ENSURING THIS KNOB IS ORTHOGONAL TO THE px TECA KNOB
            xt.TargetSet(["x", "px"], 0.0, at="qd.41310"),
            xt.TargetSet(x=0, px=0, at=xt.END),  # Ensure bump is closed
        ],
    )
    return opt


# class ActionMeasAmplDet(xt.Action):

#     def __init__(self, line, num_turns, nemitt_x, nemitt_y):

#         self.line = line
#         self.num_turns = num_turns
#         self.nemitt_x = nemitt_x
#         self.nemitt_y = nemitt_y

#     def run(self):

#         det_coefficients = self.line.get_amplitude_detuning_coefficients(
#                                 nemitt_x=self.nemitt_x, nemitt_y=self.nemitt_y,
#                                 num_turns=self.num_turns)

#         out = {'d_xx': det_coefficients['det_xx'],
#                'd_yy': det_coefficients['det_yy']}

#         return out

# action = ActionMeasAmplDet(line=line, nemitt_x=2.5e-6, nemitt_y=2.5e-6,
#                            num_turns=128)

# opt = line.match(vary=xt.VaryList(['kof.a23b1', 'kod.a23b1'], step=1.),
#                  targets=[action.target('d_xx', 1000., tol=0.1),
#                           action.target('d_yy', 2000., tol=0.1)])

# opt.target_status()
# # prints:
# #
# # Target status:
# # id state tag tol_met     residue current_val target_val description
# #  0 ON           True   0.0844456     1000.08       1000 'd_xx', val=1000, ...
# #  1 ON           True -0.00209987        2000       2000 'd_yy', val=2000, ...

# # Complete source: xtrack/examples/match/006_match_action.py




def horizontal_bumpLSS4_MD(line, x_target = (TECA.jaw - TECA.width), px_target = TECA.tilt):
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
        vary=[
            xt.VaryList(
                [
                    "kmpsh41402",
                    "kmpsh42198",
                ],
                step=1e-9, limits=[-5e-4, 5e-4]
            ),
            xt.VaryList(
                [
                    "kmplh41658",
                    "kmplh41994",
                ],
                step=1e-9, limits=[-7e-4, 7e-4]
            ),
        ],  # ADD A COSTRAINT THAT LIMITS THE SQUARE SUMM OF THE STRENGTHS
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START), 
            # xt.TargetSet(x=x_target, px=px_target, at='TECA.entry'),  # Apply bump at target location, JUST THE X POSITION
            xt.TargetSet(x = - 40e-3, px_target = px_target, at='TECA.entry'),  # Apply bump at target location, JUST THE X POSITION
            
            xt.TargetSet(x=0, px=0, at=xt.END)  # Ensure bump is closed
        ]
    )
    return opt


def horizontal_bumpLSS4(line, x_target = (TECA.jaw - TECA.width), px_target = TECA.tilt):
    """
    Adjusts the horizontal bump at tpst.21760_entry using four elements.
    
    Parameters:
        line (xt.Line): The accelerator beamline.
        x_target (float): Desired horizontal position at tpst.21760_entry.
        px_target (float): Desired horizontal angle at tpst.21760_entry.
    """
    opt = line.match(
        start='mpsh.41402', end='mpsh.42198',  # Bump region
        #start='begi.10010', end='end.10010',
        betx=1, bety=1, x=0, px=0,  # Keep initial conditions unchanged
        vary=xt.VaryList(
            ['kmpsh41402', 'kmplh41658', 'kmplh41994', 'kmpsh42198'],  # Selected correctors
            #['kmplh41658', 'kmplh41994', 'kmpsh42198'],  # Selected correctors with one corrector removed (kmpsh41402)
            step=1e-10,
            limits=[-5e-4, 5e-4]  # Define kick limits to avoi excessive changes
        ),  # ADD A COSTRAINT THAT LIMITS THE SQUARE SUMM OF THE STRENGTHS
        targets=[
            xt.TargetSet(x=0, px=0, at=xt.START), 
            #xt.TargetSet(x = x_target, px = px_target, at='TECA.entry'),  # Apply bump at target location
            xt.Target(x = x_target, at='TECA.entry'),  # Apply bump at target location
            #xt.TargetSet(x = x_target, px = 0, at='TECA.entry'),  # Apply bump at target location, JUST THE X POSITION!!! AS YANN SUGGESTED INTHE CCC THE 14.07.2025
            #xt.TargetSet(x=0, px=0, at= 'drift_779'),  # Ensure bump is closed
            xt.TargetSet(x=0, px=0, at=xt.END)  # Ensure bump is closed
        ]
    )
    return opt