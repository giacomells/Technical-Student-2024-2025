"""
run_crystal_characterisation.py
================================
Sweep crystal bending radius and kick angle, track single-pass particles,
fit the angular-kick distribution with a double Gaussian, and report
channeling/loss fractions both on the terminal and as saved plots.

Method (from crystal_parametrisation.ipynb)
--------------------------------------------
1. Build a DummyCrystal -> EverestCrystal (xcoll) for each configuration.
2. Fire n_particles perfectly-aligned particles through it.
3. Separate lost vs surviving particles; discard unphysical large-angle tail.
4. Fit the px histogram with gauss_amorphous + gauss_channeled.
5. Classify surviving particles into: channeled / amorphous / dechanneled.
6. Print a compact table per configuration; save all plots.

Usage
-----
    python run_crystal_characterisation.py
    python run_crystal_characterisation.py --n-particles 20000
    python run_crystal_characterisation.py --output-dir results/crystal

Requirements
------------
    xsuite + xcoll environment (xsuite_env)
    No network access or JSON file needed — DummyCrystal is standalone.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_REPO_ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import xpart as xp

import crystal_extraction.plotters as plotters
from crystal_extraction.dummy_crystal import DummyCrystal

# ---------------------------------------------------------------------------
# Parameter grid
# ---------------------------------------------------------------------------

BENDING_RADII = np.array([5, 10, 15, 20, 25, 30])                        # m
KICKS         = np.array([100, 150, 200, 300, 400, 500, 700, 1000]) * 1e-6  # rad
SIGMA_CUTOFF  = 5    # n-sigma window for classifying particles
NBINS         = 200
P0C           = 400e9  # GeV/c


def _section(title: str) -> None:
    print("\n" + "-" * 64)
    print(f"  {title}")
    print("-" * 64)


# ---------------------------------------------------------------------------
# Double-Gaussian fit helpers
# ---------------------------------------------------------------------------

def _gauss(x, std, amp, mean=0.0):
    return amp * np.exp(-0.5 * ((x - mean) / std) ** 2) / (std * np.sqrt(2 * np.pi))


def _double_gauss(x, m1, s1, a1, m2, s2, a2):
    return _gauss(x, s1, a1, mean=m1) + _gauss(x, s2, a2, mean=m2)


def fit_px_distribution(rest_pxs: np.ndarray, kick: float, cry: DummyCrystal):
    hist, bins = np.histogram(rest_pxs, bins=NBINS)
    centres    = 0.5 * (bins[1:] + bins[:-1])
    p0 = [0.0, cry.sigma_kick_ac,   float(np.max(hist)),
          kick, cry.sigma_kick_chan, float(np.max(hist))]
    try:
        params, _ = curve_fit(_double_gauss, centres, hist,
                               p0=p0, maxfev=10_000)
    except RuntimeError:
        params = np.array(p0)   # fall back to initial guess
    return params, hist, centres


# ---------------------------------------------------------------------------
# Single-configuration tracking and classification
# ---------------------------------------------------------------------------

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

    all_pxs  = particles.px
    lost     = particles.state <= 0
    rest_pxs = all_pxs[~lost]

    # Remove unphysical large-angle tail before fitting
    large = (rest_pxs > kick + 1e-4) | (rest_pxs < -1e-4)
    rest_pxs_clean = rest_pxs[~large]

    params, hist, centres = fit_px_distribution(rest_pxs_clean, kick, cry)
    ac_mu, ac_sig, _, ch_mu, ch_sig, _ = params

    lost_pct   = np.sum(lost)  / n_particles
    large_pct  = np.sum(large) / n_particles
    chan_pct   = np.sum(np.abs(rest_pxs_clean - kick)  <= SIGMA_CUTOFF * abs(ch_sig)) / n_particles
    ac_pct     = np.sum(np.abs(rest_pxs_clean - ac_mu) <= SIGMA_CUTOFF * abs(ac_sig)) / n_particles
    dechan_pct = np.sum(
        (rest_pxs_clean - ac_mu > SIGMA_CUTOFF * abs(ac_sig)) &
        (rest_pxs_clean - kick  < -SIGMA_CUTOFF * abs(ch_sig))
    ) / n_particles

    return {
        "chan_pct":    chan_pct,
        "ac_pct":     ac_pct,
        "lost_pct":   lost_pct,
        "large_pct":  large_pct,
        "dechan_pct": dechan_pct,
        "ac_mu":      ac_mu,   "ac_sig": ac_sig,
        "ch_mu":      ch_mu,   "ch_sig": ch_sig,
        "params":     params,
        "hist":       hist,
        "centres":    centres,
    }


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def save_fit_plots(bending_radius: float, results: list[dict],
                   output_dir: Path) -> Path:
    """Save per-kick angular-kick fit plots for one bending radius."""
    f, ax = plotters.subplots(nsubplots=len(KICKS))
    f.suptitle(f"Angular-kick distributions  —  R = {bending_radius} m", fontsize=10)

    for i, (kick, d) in enumerate(zip(KICKS, results)):
        y_fit = _double_gauss(d["centres"], *d["params"])
        ax[i].bar(d["centres"], d["hist"],
                  width=d["centres"][1] - d["centres"][0],
                  alpha=0.7, label="data")
        ax[i].plot(d["centres"], y_fit, color="red", lw=1.5, label="fit")
        ax[i].set_title(f"kick = {kick*1e6:.0f} urad", fontsize=8)
        ax[i].set_xlabel("px  [rad]")
        ax[i].set_ylabel("counts")
        ax[i].legend(fontsize=7)

    f.tight_layout()
    fname = output_dir / f"fits_R{int(bending_radius):02d}m.png"
    f.savefig(fname, dpi=110, bbox_inches="tight")
    plt.close(f)
    return fname


def save_summary_plots(results_by_radius: dict, output_dir: Path) -> None:
    """Save fraction and fit-statistic summary plots across all radii."""
    kicks_urad = KICKS * 1e6

    f_pct,   ax_pct   = plotters.subplots(nsubplots=5, col_numbers=[2])
    f_stats, ax_stats = plotters.subplots(nsubplots=4, col_numbers=[2])

    pct_keys   = ["chan_pct", "ac_pct", "lost_pct", "large_pct", "dechan_pct"]
    pct_labels = ["Channeled %", "Amorphous %", "Lost %",
                  "Large-angle %", "Dechanneled %"]
    stat_keys   = ["ac_mu", "ch_mu", "ac_sig", "ch_sig"]
    stat_labels = ["AC mean [urad]", "Ch mean [urad]",
                   "AC sigma [urad]", "Ch sigma [urad]"]

    for r, data_list in results_by_radius.items():
        for ax, key, label in zip(ax_pct, pct_keys, pct_labels):
            ax.plot(kicks_urad, [d[key] * 100 for d in data_list],
                    marker="o", label=f"R={r} m")
            ax.set_ylabel(label)
        for ax, key, label in zip(ax_stats, stat_keys, stat_labels):
            scale = 1e6 if "mu" in key or "sig" in key else 1.0
            ax.plot(kicks_urad, [d[key] * scale for d in data_list],
                    marker="o", label=f"R={r} m")
            ax.set_ylabel(label)

    for ax in ax_pct:
        ax.set_xlabel("Kick angle [urad]")
    ax_pct[0].legend(title="Bending radius", fontsize=7)

    for ax in ax_stats:
        ax.set_xlabel("Kick angle [urad]")
    ax_stats[0].legend(title="Bending radius", fontsize=7)

    f_pct.suptitle("Crystal characterisation — particle fractions", fontsize=10)
    f_stats.suptitle("Crystal characterisation — Gaussian fit statistics", fontsize=10)
    f_pct.tight_layout();   f_stats.tight_layout()

    f_pct.savefig(output_dir / "crystal_fractions.png",  dpi=120, bbox_inches="tight")
    f_stats.savefig(output_dir / "crystal_fit_stats.png", dpi=120, bbox_inches="tight")
    plt.close(f_pct);  plt.close(f_stats)


# ---------------------------------------------------------------------------
# Terminal summary table
# ---------------------------------------------------------------------------

def print_table(results_by_radius: dict) -> None:
    _section("Summary table")
    header = (f"  {'R[m]':>5}  {'kick[urad]':>10}  "
              f"{'chan%':>6}  {'ac%':>6}  {'lost%':>6}  {'dechan%':>7}")
    print(header)
    print("  " + "-" * (len(header) - 2))
    for r, data_list in results_by_radius.items():
        for kick, d in zip(KICKS, data_list):
            print(f"  {r:>5}  {kick*1e6:>10.0f}  "
                  f"{d['chan_pct']*100:>6.1f}  "
                  f"{d['ac_pct']*100:>6.1f}  "
                  f"{d['lost_pct']*100:>6.1f}  "
                  f"{d['dechan_pct']*100:>7.1f}")
        print()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Crystal single-pass characterisation: efficiency vs kick angle"
    )
    parser.add_argument("--n-particles", type=int, default=50_000,
                        help="Particles per configuration (default: 50000)")
    parser.add_argument("--output-dir", default="outputs/crystal",
                        help="Directory for saved plots (default: outputs/crystal)")
    args = parser.parse_args()

    output_dir  = Path(args.output_dir)
    n_particles = args.n_particles
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 64)
    print("  Crystal single-pass characterisation")
    print(f"  Bending radii : {BENDING_RADII} m")
    print(f"  Kick angles   : {KICKS*1e6} urad")
    print(f"  Particles     : {n_particles:,} per configuration")
    print(f"  Configurations: {len(BENDING_RADII) * len(KICKS)}")
    print("=" * 64)

    cry = DummyCrystal()
    cry.active_length = 2.5e-3

    results_by_radius: dict = {}

    for r in BENDING_RADII:
        _section(f"Bending radius R = {r} m")
        cry.bending_radius = r
        data_list: list[dict] = []

        for kick in KICKS:
            print(f"  kick {kick*1e6:6.0f} urad  "
                  f"length {cry.active_length*1e3:.2f} mm ...",
                  end=" ", flush=True)
            d = characterise_one(cry, kick, n_particles)
            data_list.append(d)
            print(f"chan {d['chan_pct']*100:5.1f}%  "
                  f"ac {d['ac_pct']*100:5.1f}%  "
                  f"lost {d['lost_pct']*100:5.1f}%")

        results_by_radius[r] = data_list

        fname = save_fit_plots(r, data_list, output_dir)
        print(f"  -> fit plots: {fname}")

    print_table(results_by_radius)

    _section("Saving summary plots")
    save_summary_plots(results_by_radius, output_dir)
    print(f"  crystal_fractions.png and crystal_fit_stats.png saved")

    print("\n" + "=" * 64)
    print(f"  Done.  All plots in: {output_dir.resolve()}")
    print("=" * 64 + "\n")


if __name__ == "__main__":
    main()
