# Theory and Background

## Slow Extraction from the SPS

The SPS (Super Proton Synchrotron) at CERN deflects a stored proton beam
toward the North Area (NA) fixed-target experiments via a process called
**slow extraction**. Instead of emptying the ring in a single turn, the beam
is pushed gradually across a thin electrostatic septum over thousands of turns,
producing a quasi-continuous spill.

## Resonant Extraction (third-integer resonance)

The standard SPS slow extraction exploits the horizontal third-integer resonance
($Q_x \approx 26/3$). Sextupole magnets create a stable triangle in phase space
(the **separatrix**). Particles outside the separatrix spiral outward and are
captured by the Zero-degree Electrostatic Septum (ZS), then bent by the magnetic
septa MST and MSE and transported to TT20.

Key parameters:
- **Beam momentum**: 400 GeV/c
- **Beam rigidity**: $B\rho = p / (qc) \approx 1334~\text{T·m}$
- **Normalised horizontal emittance**: $\varepsilon_x^N = 1\!\times\!10^{-6}~\text{m·rad}$
- **Momentum spread**: $\Delta p/p = 1.5\!\times\!10^{-3}$

The septum chain typically consists of:
| Element | Type | Function |
|---------|------|----------|
| ZS (×5) | Electrostatic | First deflection, thin wire ~0.3 mm |
| TPST | Magnetic | Pre-MST bump correction |
| MST (×3) | Magnetic | Main deflection |
| MSE (×5) | Magnetic | Final deflection to TT20 |

## Non-Resonant Extraction with TECA Crystal

An alternative extraction scheme studied here replaces the resonant excitation
with a bent silicon crystal (**TECA** — Test Crystal in the SPS) positioned
inside the circulating beam.  Channeled protons are deflected by the crystal
lattice planes at a bending angle sufficient to steer them into the extraction
septum aperture.

The crystal is modelled as an `EverestCrystal` element from `xcoll`, using the
Geant4-based Everest physics model.

## Optics

Studies are performed at the nominal SPS flat-top optics ($Q_x \approx 26.13$,
$Q_y \approx 26.18$, labelled `lhc_q22`) and, for comparison, at the Q22 optics
($Q_x \approx 20.13$, $Q_y \approx 20.18$).  Machine sequences are stored as
`xtrack` JSON files (`sps_for_sx.json`, `lhc_q22.json`).

## References

- Xsuite documentation: <https://xsuite.readthedocs.io>
- xcoll documentation: <https://xcoll.readthedocs.io>
- SPS slow extraction review: CERN-ACC-2018-0009
