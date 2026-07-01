# Notation - ICFF (Katritch, Totrov, Abagyan 2003)

Units are **kcal/mol** for energy, **angstrom (Å)** for distance, **degrees** (grid) or radians (trig argument) for angles throughout. No reduced units. Force field is MMFF94s (Halgren).

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $E_{\text{tor}}$ | ICFF composite torsion energy for one bond | kcal/mol | eq. 1; absorbs stretch/bend/OOP/1-4 vdW implicitly |
| $\theta$ | torsion angle about a rotatable bond | deg (grid) / rad (trig) | scanned 0°..360° |
| $C_0$ | constant offset in torsion Fourier series | kcal/mol | eq. 1 |
| $A_k$ | Fourier cosine coefficient, harmonic $k$ | kcal/mol | $k=0..6$ (sixfold) |
| $B_k$ | Fourier sine coefficient, harmonic $k$ | kcal/mol | $k=0..6$ |
| $k$ | Fourier harmonic index | int, $0 \le k \le 6$ | small coeffs (<0.01) nullified |
| $E_{\text{vw}(1\text{-}5,1\text{-}6)}$ | soft empirical repulsion for 1-5 / 1-6 pairs | kcal/mol | eq. 2; only for $R_{IJ}<R_{IJ}^*$ |
| $\varepsilon_{IJ}$ | vdW well depth for atom types $I,J$ | kcal/mol | minimum-energy magnitude of source vdW surface |
| $R_{IJ}$ | actual distance between atoms $I,J$ | Å | - |
| $R_{IJ}^*$ | distance at vdW energy minimum | Å | from source MMFF94 vdW params |
| $R_{IJ}^0$ | distance where vdW surface crosses zero | Å | from source MMFF94 vdW params |
| $C$ | harmonic/cubic blend factor in eq. 2 | dimensionless, $0\le C\le1$ | only adjustable ICFF param; chosen $C=0.55$ |
| $E_{\text{MMFF}}$ | full MMFF94(s) Cartesian energy | kcal/mol | eq. 3; source force field |
| $EB, EA, EBA, EOOP, ET, EvdW, EQ$ | MMFF term components | kcal/mol | bond/angle/stretch-bend/OOP/torsion/vdW/electrostatic |
| $E_r$ | torsion restraint potential (profile generation) | kcal/mol | harmonic about $\theta^0$ |
| $\theta^0$ | target/constrained torsion angle | deg | grid value |
| $C_r$ | restraint force constant | 10000 kcal (per angle^2) | keeps deviation < 0.1° |
| $W(\theta)$ | least-squares weight for Fourier fit | dimensionless | $[E-\min E+1]^{-1}$, inverse relative energy |
| $E_{\text{cut}}$ | energy cutoff for RMSD reporting | kcal/mol | e.g. 1,3,5,7,10,20 kcal |
| "1-4" | atoms separated by 3 covalent bonds | - | absorbed into ICFF torsion term |
| "1-5","1-6" | atoms separated by 4 / 5 covalent bonds | - | soft repulsion eq. 2 |
| "1-7+" (6+ bonds) | atoms separated by >= 6 bonds | - | original MMFF94 hard vdW term |

## Conventions / regime notes
- ICFF is a **rigid-covalent-geometry** internal-coordinate force field: bond lengths, bond angles, phase angles are FIXED; only torsion angles are free variables. This reduces free variables ~10-fold vs Cartesian.
- Torsion parameters are derived per **bond type** (defined by atom types + connection topology of the torsion fragment), generated on-the-fly at ~0.1 s per drug-like molecule.
- Fixed covalent geometry is obtained by Cartesian minimization of the "local" MMFF94s energy (bonded terms + 1-4 vdW only; no long-range vdW, no electrostatics), starting from best-energy torsion values.
- Tri-coordinate nitrogens are forced planar in ICFF (source of some geometry deviation).
- Hydrogen-bond atom pairs keep the original "hard" vdW term (not the soft eq. 2), to preserve donor-acceptor attraction.
- The soft eq. 2 mimics only *average* bond flexibility over 4-5 bonds.
