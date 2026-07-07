"""@author: Y. Dutheil"""

import typing as t
from math import ceil

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def ellipse_frm_cov(cov: np.ndarray, nsig: float=1.) -> t.Tuple[float, float, float]:
    r"""
    Retrurns width, hight and angle of ellipse, from the covariance matrix

    taken from `https://stackoverflow.com/questions/12301071/multidimensional-confidence-intervals`

    Parameters
    ----------
    cov : 2-D array or list
        Covariance matrix

    msig : scalar
        To be checked and explained, and illustrated what this nsig is

    Returns
    -------
    list
        list of width, heigh and rotation of the ellipse

    """
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    
    vals, vecs = vals[order], vecs[:,order]
    rotation = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * nsig * np.sqrt(vals)

    return width, height, rotation


def draw_ellipse(alpha: float, beta: float, eps: float, x0: t.Tuple[float, float]=(0, 0), ax: t.Optional[plt.Axes]=None, plot_kwargs: t.Optional[dict]=None) -> plt.Axes:
    """
    Draws ellipse of emittance 'eps' given a center and Courant Snyder parameters

    Parameters
    ----------

    alpha: float
        Courant Snyder alpha paramter
    beta: float
        Courant Snyder beta parameter
    eps: float
        Beam emittance
    ax: matplotlib Axes
        Axes on which to draw the ellipse. If None, axes is created
    x0: 2-tuple
        Center of ellipse

    Returns
    -------

    ax: matplotlib Axes
        Axes containing drawing of the ellipse

    """

    default_plot_kwargs = {
        "facecolor": None,
        "fill": None,
        "color": "black",
        "linewidth": 2,
        "alpha": 0.8,
        "zorder": 10,
    }

    if plot_kwargs is None:
        plot_kwargs = default_plot_kwargs
    else:
        default_plot_kwargs.update(plot_kwargs)
        plot_kwargs = default_plot_kwargs

    if ax is None:
        f, ax = plt.subplots()

    cov = eps*np.array([[beta, -alpha],[-alpha, (1+alpha**2)/beta]])
    width, height, angle = ellipse_frm_cov(cov, nsig=1)
    ax.add_patch(
        mpl.patches.Ellipse(x0, width=width, height=height, angle=angle, **plot_kwargs)
    )

    return ax


def my_mpl_style(smooth=False):
    '''Sets my preferred style options for matplotlib.'''

    import matplotlib as mpl
    from cycler import cycler

    cmap = 'cool'

    if smooth:
        colors = ['skyblue', 'dodgerblue', 'b',
                  'indigo', 'darkmagenta', 'fuchsia', 'deeppink', 'black']
        # Color choices
        mpl.rcParams['axes.prop_cycle'] = cycler(color=colors)

    mpl.rcParams['image.cmap'] = cmap

    # Font
    mpl.rcParams["mathtext.fontset"] = 'cm'
    mpl.rcParams["figure.facecolor"] = 'white'
    mpl.rcParams['text.latex.preamble'] = r'\boldmath'
    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['legend.fontsize'] = 12
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    mpl.rcParams['axes.formatter.limits'] = (-3, 4)


def subplots(nrows=1, ncols=1, width_single=4, height_single=3, nsubplots=None, col_numbers=None, **kwargs):
    if col_numbers is None:
        col_numbers = [3, 4, 5]
    if nsubplots is not None:
        min_subplots = [ceil(nsubplots/col_number)*col_number for col_number in col_numbers]
        chosen_col_numbers = [col_number for col_number, min_subplot
                              in zip(col_numbers, min_subplots)
                              if min_subplot == min(min_subplots)]
        ncols = max(chosen_col_numbers)
        nrows = ceil(nsubplots/ncols)
    f, ax = plt.subplots(nrows, ncols,
                         figsize=(ncols*width_single, nrows*height_single), **kwargs)
    if type(ax) != np.ndarray:
        ax = np.array([ax])
    ax = ax.flatten()
    return f, ax
