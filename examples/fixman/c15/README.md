# C15 idealized bead chain

Linear 15-bead chain, planar zigzag starting geometry, per
`docs/specs/fixman-idealized-chains-validation.md` Sec. 3.1. The heaviest of
the four idealized chains -- the Tier-1 capstone (T1.3).

| Spec Sec. 3.1 parameter | Value | Where |
|---|---|---|
| Beads | 15 | `c15.mol2` `@<TRIPOS>ATOM` |
| Bonds | 14 | `c15.mol2` `@<TRIPOS>BOND` |
| Free torsions | 12 | implied: 12 non-terminal bonds |
| Bond angle | 109 deg | `c15.frcmod` `ANGLE` |
| Bead mass | 14 amu | `c15.frcmod` `MASS` |
| Bond length | 1.54 A | `c15.mol2` coordinates, `c15.frcmod` `BOND` |
| Bond spring | 83.66 kcal/mol/A^2 | `c15.frcmod` `BOND` |
| Angle spring | 43.46 kcal/mol | `c15.frcmod` `ANGLE` |
| Torsion term | k=0 (identically zero) | `c15.frcmod` `DIHE` |
| Nonbonded | sigma/epsilon = 0, no charges | `c15.frcmod` `NONBON`, `c15.mol2` `NO_CHARGES` |

Coordinates are a planar zigzag (all torsions == 0 at the reference geometry),
bond length 1.54 A, interior bond angle 109 deg at every vertex.

Regenerate with:

```
$CONDA_PREFIX/bin/tleap -f examples/fixman/c15/tleap.in
```
