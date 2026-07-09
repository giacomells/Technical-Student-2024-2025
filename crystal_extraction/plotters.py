from __future__ import annotations

"""@author: Y. Dutheil"""

import typing as t
from math import ceil
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
import xpart as xp

from crystal_extraction.dummy_crystal import DummyCrystal


# ---------------------------------------------------------------------------
# Crystal characterisation defaults
# ---------------------------------------------------------------------------

BENDING_RADII = np.array([5, 10, 15, 20, 25, 30])
KICKS = np.array([100, 150, 200, 300, 400, 500, 700, 1000]) * 1e-6
SIGMA_CUTOFF = 5
NBINS = 200
P0C = 400e9

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


def section(title: str) -> None:
    print("\n" + "-" * 64)
    print(f"  {title}")
    print("-" * 64)


def _gauss(x, std, amp, mean=0.0):
    return amp * np.exp(-0.5 * ((x - mean) / std) ** 2) / (std * np.sqrt(2 * np.pi))


def _double_gauss(x, m1, s1, a1, m2, s2, a2):
    return _gauss(x, s1, a1, mean=m1) + _gauss(x, s2, a2, mean=m2)


def fit_px_distribution(rest_pxs: np.ndarray, kick: float, cry: DummyCrystal):
    hist, bins = np.histogram(rest_pxs, bins=NBINS)
    centres = 0.5 * (bins[1:] + bins[:-1])
    p0 = [0.0, cry.sigma_kick_ac, float(np.max(hist)),
          kick, cry.sigma_kick_chan, float(np.max(hist))]
    try:
        params, _ = curve_fit(_double_gauss, centres, hist, p0=p0, maxfev=10_000)
    except RuntimeError:
        params = np.array(p0)
    return params, hist, centres


def characterise_one(cry: DummyCrystal, kick: float, n_particles: int) -> dict:
    cry.mu_kick_chan = kick
    particles = xp.Particles(
        p0c=P0C,
        x=np.zeros(n_particles), px=np.zeros(n_particles),
        y=np.zeros(n_particles), py=np.zeros(n_particles),
        zeta=np.zeros(n_particles), delta=np.zeros(n_particles),
    )
    particles._init_random_number_generator()
    xsuite_cry = cry.build_xsuite_crystal()
    xsuite_cry.track(particles)

    all_pxs = particles.px
    lost = particles.state <= 0
    rest_pxs = all_pxs[~lost]

    large = (rest_pxs > kick + 1e-4) | (rest_pxs < -1e-4)
    rest_pxs_clean = rest_pxs[~large]

    params, hist, centres = fit_px_distribution(rest_pxs_clean, kick, cry)
    ac_mu, ac_sig, _, ch_mu, ch_sig, _ = params

    lost_pct = np.sum(lost) / n_particles
    large_pct = np.sum(large) / n_particles
    chan_pct = np.sum(np.abs(rest_pxs_clean - kick) <= SIGMA_CUTOFF * abs(ch_sig)) / n_particles
    ac_pct = np.sum(np.abs(rest_pxs_clean - ac_mu) <= SIGMA_CUTOFF * abs(ac_sig)) / n_particles
    dechan_pct = np.sum(
        (rest_pxs_clean - ac_mu > SIGMA_CUTOFF * abs(ac_sig))
        & (rest_pxs_clean - kick < -SIGMA_CUTOFF * abs(ch_sig))
    ) / n_particles

    return {
        "chan_pct": chan_pct,
        "ac_pct": ac_pct,
        "lost_pct": lost_pct,
        "large_pct": large_pct,
        "dechan_pct": dechan_pct,
        "ac_mu": ac_mu,
        "ac_sig": ac_sig,
        "ch_mu": ch_mu,
        "ch_sig": ch_sig,
        "params": params,
        "hist": hist,
        "centres": centres,
    }


def save_fit_plots(bending_radius: float, results: list[dict], output_dir: Path) -> Path:
    f, ax = subplots(nsubplots=len(KICKS))
    f.suptitle(f"Angular-kick distributions  -  R = {bending_radius} m", fontsize=10)

    for i, (kick, d) in enumerate(zip(KICKS, results)):
        y_fit = _double_gauss(d["centres"], *d["params"])
        ax[i].bar(
            d["centres"],
            d["hist"],
            width=d["centres"][1] - d["centres"][0],
            alpha=0.7,
            label="data",
        )
        ax[i].plot(d["centres"], y_fit, color="red", lw=1.5, label="fit")
        ax[i].set_title(f"kick = {kick * 1e6:.0f} urad", fontsize=8)
        ax[i].set_xlabel("px  [rad]")
        ax[i].set_ylabel("counts")
        ax[i].legend(fontsize=7)

    f.tight_layout()
    fname = output_dir / f"fits_R{int(bending_radius):02d}m.png"
    f.savefig(str(fname), dpi=110, bbox_inches="tight")
    plt.close(f)
    return fname


def save_summary_plots(results_by_radius: dict, output_dir: Path) -> None:
    kicks_urad = KICKS * 1e6

    f_pct, ax_pct = subplots(nsubplots=5, col_numbers=[2])
    f_stats, ax_stats = subplots(nsubplots=4, col_numbers=[2])

    pct_keys = ["chan_pct", "ac_pct", "lost_pct", "large_pct", "dechan_pct"]
    pct_labels = ["Channeled %", "Amorphous %", "Lost %", "Large-angle %", "Dechanneled %"]
    stat_keys = ["ac_mu", "ch_mu", "ac_sig", "ch_sig"]
    stat_labels = ["AC mean [urad]", "Ch mean [urad]", "AC sigma [urad]", "Ch sigma [urad]"]

    for r, data_list in results_by_radius.items():
        for ax, key, label in zip(ax_pct, pct_keys, pct_labels):
            ax.plot(kicks_urad, [d[key] * 100 for d in data_list], marker="o", label=f"R={r} m")
            ax.set_ylabel(label)
        for ax, key, label in zip(ax_stats, stat_keys, stat_labels):
            scale = 1e6 if "mu" in key or "sig" in key else 1.0
            ax.plot(kicks_urad, [d[key] * scale for d in data_list], marker="o", label=f"R={r} m")
            ax.set_ylabel(label)

    for ax in ax_pct:
        ax.set_xlabel("Kick angle [urad]")
    ax_pct[0].legend(title="Bending radius", fontsize=7)

    for ax in ax_stats:
        ax.set_xlabel("Kick angle [urad]")
    ax_stats[0].legend(title="Bending radius", fontsize=7)

    f_pct.suptitle("Crystal characterisation - particle fractions", fontsize=10)
    f_stats.suptitle("Crystal characterisation - Gaussian fit statistics", fontsize=10)
    f_pct.tight_layout()
    f_stats.tight_layout()

    f_pct.savefig(str(output_dir / "crystal_fractions.png"), dpi=120, bbox_inches="tight")
    f_stats.savefig(str(output_dir / "crystal_fit_stats.png"), dpi=120, bbox_inches="tight")
    plt.close(f_pct)
    plt.close(f_stats)


def print_table(results_by_radius: dict) -> None:
    section("Summary table")
    header = (f"  {'R[m]':>5}  {'kick[urad]':>10}  "
              f"{'chan%':>6}  {'ac%':>6}  {'lost%':>6}  {'dechan%':>7}")
    print(header)
    print("  " + "-" * (len(header) - 2))
    for r, data_list in results_by_radius.items():
        for kick, d in zip(KICKS, data_list):
            print(f"  {r:>5}  {kick * 1e6:>10.0f}  "
                  f"{d['chan_pct'] * 100:>6.1f}  "
                  f"{d['ac_pct'] * 100:>6.1f}  "
                  f"{d['lost_pct'] * 100:>6.1f}  "
                  f"{d['dechan_pct'] * 100:>7.1f}")
        print()


def run_crystal_characterisation(n_particles: int, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 64)
    print("  Crystal single-pass characterisation")
    print(f"  Bending radii : {BENDING_RADII} m")
    print(f"  Kick angles   : {KICKS * 1e6} urad")
    print(f"  Particles     : {n_particles:,} per configuration")
    print(f"  Configurations: {len(BENDING_RADII) * len(KICKS)}")
    print("=" * 64)

    cry = DummyCrystal()
    cry.active_length = 2.5e-3

    results_by_radius: dict = {}

    for r in BENDING_RADII:
        section(f"Bending radius R = {r} m")
        cry.bending_radius = r
        data_list: list[dict] = []

        for kick in KICKS:
            print(f"  kick {kick * 1e6:6.0f} urad  "
                  f"length {cry.active_length * 1e3:.2f} mm ...", end=" ", flush=True)
            d = characterise_one(cry, kick, n_particles)
            data_list.append(d)
            print(f"chan {d['chan_pct'] * 100:5.1f}%  "
                  f"ac {d['ac_pct'] * 100:5.1f}%  "
                  f"lost {d['lost_pct'] * 100:5.1f}%")

        results_by_radius[r] = data_list

        fname = save_fit_plots(r, data_list, output_dir)
        print(f"  -> fit plots: {fname}")

    print_table(results_by_radius)

    section("Saving summary plots")
    save_summary_plots(results_by_radius, output_dir)
    print("  crystal_fractions.png and crystal_fit_stats.png saved")

    print("\n" + "=" * 64)
    print(f"  Done.  All plots in: {output_dir.resolve()}")
    print("=" * 64 + "\n")
