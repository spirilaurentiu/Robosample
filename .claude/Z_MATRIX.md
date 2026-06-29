# Z-matrix

A molecule's internal coordinates are defined by a **Z-matrix**: each atom is placed by a bond length `r`, a bond angle `theta`, and a dihedral (torsion) `tau` relative to three previously placed atoms. The Z-matrix therefore induces the **bond-angle-torsion (BAT)** decomposition `q_int = (b, theta, tau)` and quivalently, the kinematic tree -- each joint of the tree corresponds to one internal coordinate. Torsional dynamics is the case in which `b` and `theta` are frozen (their joints welded) and only the `tau` (a subset of the Z-matrix dihedrals) remain mobile.

The BAT coordinates connect the program's reduced sampling to the Cartesian partition function. The configurational partition function is `Z = integral exp(-beta*U) dx` over Cartesian `x`. The Cartesian -> BAT change of variables has the Jacobian

```text
dx = |J_{BAT}| * db * d\theta * d\tau * de
|J_BAT| = ( product_i r_i^2 ) * ( product_j sin theta_j ) * const
```

(the standard Z-matrix volume element; `de` are external DOF). This Jacobian is **independent of the torsions**, which is why torsional dynamics can target the correct conditional without an extra positional Jacobian. The remaining configuration-dependence that *does* matter -- the mass-metric determinant `det M(q)` -- is handled by the Fixman term.

The builder returns four equal-length lists `(z_i, z_j, z_k, z_l)`, one row per atom, in local atom indices, with the sentinel -1 where a reference does not yet exist. Row `r` describes the atom placed at step r and the already-built atoms it is measured against:

- `z_i`: atom placed at step `r`
- `z_j`: bond reference (bonded to `z_i`; `-1` for `r=0`)
- `z_k`: angle reference (bonded to `z_j`; `-1` for `r<2`)
- `z_l`: dihedral reference (bonded to `z_k`; `-1` for `r<3`)

Internal coordinate(s) of `z_i`:

- Bond length `(z_i, z_j)`
- Angle `(z_i, z_j, z_k)`
- Torsion `(z_i, z_j, z_k, z_l)`

The four atoms of a full row are distinct and bonded in sequence `z_i-z_j-z_k-z_l`, and every reference points **toward atoms already placed** (toward the root). The Z-matrix is therefore a strict build order: an atom is positioned only after the three atoms it is defined against. This is the *references toward the root* convention of the BAT construction (Chang, Potter & Gilson 2003; Hikiri, Yoshidome & Ikeguchi 2016) and of the AlGDock/MDAnalysis BAT implementation it derives from.

The first three atoms form a root triplet that fixes the molecular frame rather than torsions: they define the translation, the first two bond lengths, the first bond angle. In the reference BAT these three atoms also carry the six external DOF -- the first atom's Cartesian position plus the axis-angle rotation (polar, azimuthal, and a third angle). In Robosample those six are the root body's external joint DOF, with `U_Jacobian(q)` restoring their Haar measure.

Atoms are placed by a mass-prioritized walk outward from the root -- heaviest-first, ties broken by ascending atom index for reproducibility. The angle reference `z_k` must be non-terminal (degree > 1 in the full bond graph): a leaf cannot anchor a stable angle/torsion chain, so a terminal atom is never used as `z_k`. `z_l` must differ from `z_j`. When a row's references are not yet placed it is deferred and retried, which is what lets a single outward pass yield a consistent order without back-tracking.

The walk runs on the molecular bond graph with ring-closing bonds deleted -- the same acyclic tree used for the rigid-body forest and the frame
build. No Z-matrix reference is ever taken across a ring-closing bond; those bonds re-enter only as explicit holonomic constraints. Non-terminality, however, is judged on the *full* bond graph, so a ring atom is correctly treated as non-terminal even after its closing bond is removed.

Multiple atoms can share the same central bond `(z_j, z_k)`, meaning they rotate together about that bond: the first dihedral is the proper torsion, which describes the overall bond rotation, while the remaining improper torsions measure only the relative distortions between attached atoms, separating collective motions such as methyl rotation from local deformations, much like a spinning propeller whose blades can also bend slightly independently.

The volume element is the factor that tells you how infinitesimal volumes in Cartesian space transform when you change coordinates; it is the Jacobian determinant of the transformation and ensures that probability densities remain correctly normalized under a change of variables. In molecular internal coordinates, this is crucial because naive sampling of bond lengths, angles, and torsions would otherwise distort the equilibrium distribution.

For a full system of N atoms, the BAT representation provides exactly 3N−6 internal coordinates plus 6 external rigid-body degrees of freedom, forming a complete non-redundant reparameterization of Cartesian space (3N Cartesian <-> 3N BAT); its volume element factorizes into `(product_i r_i^2) * (product_j sin theta_j)` times a constant, independent of torsions, with any remaining Jacobian correction carried by the Fixman term.
