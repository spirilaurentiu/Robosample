# C4 idealized bead chain

Mass-corrected copy of `examples/c4bond3` (bead mass 14 amu instead of
12.011), per `docs/specs/fixman-idealized-chains-validation.md` Sec. 3.1.

| Spec Sec. 3.1 parameter | Value | Where |
|---|---|---|
| Beads | 4 | `c4.mol2` `@<TRIPOS>ATOM` |
| Bonds | 3 | `c4.mol2` `@<TRIPOS>BOND` |
| Free torsions | 1 (dihedral A-B-C-D) | implied: 1 non-terminal bond (B-C) |
| Bond angle | 90 deg | `c4.frcmod` `ANGLE` |
| Bead mass | 14 amu | `c4.frcmod` `MASS` |
| Bond length | 1.54 A | `c4.mol2` coordinates, `c4.frcmod` `BOND` |
| Bond spring | 83.66 kcal/mol/A^2 | `c4.frcmod` `BOND` |
| Angle spring | 43.46 kcal/mol | `c4.frcmod` `ANGLE` |
| Torsion term | k=0 (identically zero) | `c4.frcmod` `DIHE` |
| Nonbonded | sigma/epsilon = 0, no charges | `c4.frcmod` `NONBON`, `c4.mol2` `NO_CHARGES` |

Regenerate with:

```
$CONDA_PREFIX/bin/tleap -f examples/fixman/c4/tleap.in
```

This is the Tier-1 (Python full-stack) fixture. Tier 0 (C++, always-on) does
NOT load this file -- it hand-builds an equivalent `RobotModel` directly (see
`tests/TestFixmanIdealizedChains.cpp`); Tier 1's C4 TORSIONAL histogram is
compared against the SAME external oracle (Jain 2013 eq:24) as Tier 0, which
is what ties the shipped prmtop back to the paper-exact geometry (spec Sec.
3.3 NOTE).
