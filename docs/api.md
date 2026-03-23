# API Reference

The shared Python modules `elements.py` and `optimisers.py` appear across the
study directories (`Non-Resonant Extraction/`, `OrthogonalKnobsLSS4/`,
`Optics studies/`).  Each is a self-contained copy tailored to that study.

---

## `elements.py`

### Module-level constants

| Name | Value | Description |
|------|-------|-------------|
| `p` | 400.0 GeV/c | Beam momentum |
| `Brho` | `p * 3.3356` T·m | Beam rigidity |
| `N_EX` | 1e-6 m·rad | Normalised horizontal emittance |
| `N_EY` | 5e-6 m·rad | Normalised vertical emittance |
| `BANDWIDTH` | 1e-4 | Fractional momentum bandwidth |
| `deltaP_P` | 1.5e-3 | RMS momentum spread |
| `gamma` | `400 / proton_mass_GeV` | Relativistic Lorentz factor |

---

### `class SeptumInteraction`

Thin-septum beam-interaction element compatible with `xtrack.BeamInteraction`.

```python
SeptumInteraction(
    blade_position: float = 68e-3,   # [m] transverse position of blade inner edge
    thickness: float = 0.3e-3,       # [m] blade thickness
    kick: float = 1e-3,              # [rad] angular kick applied to extracted particles
)
```

**Methods**

- `interact(particles: xp.Particles) -> None`  
  Loses particles hitting the blade; applies `kick` to particles beyond blade.

---

### `install_septa(line, install_zs=True, septum_aperture_size=68e-3)`

Installs the full septum chain (ZS, TPST, MST, MSE) into an `xtrack.Line`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `line` | — | `xt.Line` object to modify in-place |
| `install_zs` | `True` | Whether to install the five ZS septa |
| `septum_aperture_size` | 68e-3 m | Blade position for ZS septa |

Returns `list[str]` — names of all inserted septum + aperture elements.

---

### `initialise_line(json_path="sps_for_sx.json")` / `initialise_lineQ22(...)`

Loads the SPS lattice from a JSON file, sets the reference particle, installs
septa and TECA crystal, cycles the line to `TECA.entry`, and returns a
configured `xt.Line`.

---

### `draw_synoptic(ax, line, line_df)`

Draws a machine synoptic (quadrupoles, bends, sextupoles) on `ax`.

---

### `plot_twiss(fig, twiss, line)`

Plots $\beta_x$, $\beta_y$, $D_x$ and the synoptic on a three-panel figure.

---

### `print_optics_features(twiss)`

Prints tune, chromaticity, and dispersion at the TECA position.

---

## `optimisers.py`

### `match_tunes(line, qx, qy)`

Matches horizontal and vertical tunes of `line` to `(qx, qy)` using
`xtrack`'s built-in optimizer.

### `match_chromaticity(line, dqx, dqy)`

Matches first-order chromaticity; adjusts sextupole families.

### `horizontal_bumpLSS4(line)` / `horizontal_bumpLSS2(line)`

Sets up a local horizontal closed-orbit bump in LSS4 (or LSS2) and returns
the configured `xtrack.Optimize` object.

### `ensure_closed_orbit(line)`

Verifies a closed orbit exists; raises if the line has not been matched.

---

## `Animations/phaseSpaceAnimation.py`

Standalone script.  Run from the `Animations/` directory:

```bash
python phaseSpaceAnimation.py
```

Reads `sps_for_sx.json` from the same directory.  Produces a phase-space
animation of tracked particles.

---

## Notes on code duplication

`elements.py` and `optimisers.py` are intentionally copied with minor
study-specific modifications into each subdirectory.  A future refactor could
extract a shared package; the current approach keeps each study self-contained
and runnable without a package install step.
