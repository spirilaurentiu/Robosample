# C5 idealized bead chain

Linear 5-bead chain, planar zigzag starting geometry, per
`docs/specs/fixman-idealized-chains-validation.md` Sec. 3.1.

| Spec Sec. 3.1 parameter | Value | Where |
|---|---|---|
| Beads | 5 | `c5.mol2` `@<TRIPOS>ATOM` |
| Bonds | 4 | `c5.mol2` `@<TRIPOS>BOND` |
| Free torsions | 2 | implied: 2 non-terminal bonds |
| Bond angle | 109 deg | `c5.frcmod` `ANGLE` |
| Bead mass | 14 amu | `c5.frcmod` `MASS` |
| Bond length | 1.54 A | `c5.mol2` coordinates, `c5.frcmod` `BOND` |
| Bond spring | 83.66 kcal/mol/A^2 | `c5.frcmod` `BOND` |
| Angle spring | 43.46 kcal/mol | `c5.frcmod` `ANGLE` |
| Torsion term | k=0 (identically zero) | `c5.frcmod` `DIHE` |
| Nonbonded | sigma/epsilon = 0, no charges | `c5.frcmod` `NONBON`, `c5.mol2` `NO_CHARGES` |

Coordinates are a planar zigzag (all torsions == 0 at the reference geometry
-- a valid, non-overlapping starting point; the sampler explores from there),
generated with bond length 1.54 A and interior bond angle 109 deg at every
vertex.

Regenerate with:

```
$CONDA_PREFIX/bin/tleap -f examples/fixman/c5/tleap.in
```
