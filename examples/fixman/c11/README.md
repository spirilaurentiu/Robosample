# C11 idealized bead chain

Linear 11-bead chain, planar zigzag starting geometry, per
`docs/specs/fixman-idealized-chains-validation.md` Sec. 3.1.

| Spec Sec. 3.1 parameter | Value | Where |
|---|---|---|
| Beads | 11 | `c11.mol2` `@<TRIPOS>ATOM` |
| Bonds | 10 | `c11.mol2` `@<TRIPOS>BOND` |
| Free torsions | 8 | implied: 8 non-terminal bonds |
| Bond angle | 109 deg | `c11.frcmod` `ANGLE` |
| Bead mass | 14 amu | `c11.frcmod` `MASS` |
| Bond length | 1.54 A | `c11.mol2` coordinates, `c11.frcmod` `BOND` |
| Bond spring | 83.66 kcal/mol/A^2 | `c11.frcmod` `BOND` |
| Angle spring | 43.46 kcal/mol | `c11.frcmod` `ANGLE` |
| Torsion term | k=0 (identically zero) | `c11.frcmod` `DIHE` |
| Nonbonded | sigma/epsilon = 0, no charges | `c11.frcmod` `NONBON`, `c11.mol2` `NO_CHARGES` |

Coordinates are a planar zigzag (all torsions == 0 at the reference geometry),
bond length 1.54 A, interior bond angle 109 deg at every vertex.

Regenerate with:

```
$CONDA_PREFIX/bin/tleap -f examples/fixman/c11/tleap.in
```

Per spec Sec. 5 (Tier 1, T1.3), C11 MAY be dropped from the statistical tier
since C4+C5+C15 already exercise the scaling; the prmtop/rst7 are still
shipped here for completeness and for anyone extending the Tier-1 suite.
