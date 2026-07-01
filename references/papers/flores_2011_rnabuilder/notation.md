# Notation - Flores 2011, RNABuilder

Unit system (inherited from Molmodel/Simbody, reduced/knowledge-based - time,
energy, temperature not physically meaningful here):

| quantity | unit |
|---|---|
| length | nanometers (nm) |
| time | picoseconds (ps) |
| energy | kJ/mol |
| mass | Daltons |

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| q | generalized (internal) coordinates | nq-vector | joint coords relating each body to its parent; q = {q_i} |
| u | generalized speeds (internal velocities) | n-vector | u = {u_i} |
| q_i, u_i | coords/speeds of the i-th mobilizer | 1-6 each | number = DOFs granted by that mobilizer |
| nq | number of internal coordinates | int | nq >= n (can exceed n) |
| n | number of internal velocities / mobilities (DOFs) | int | O(n) time scaling in n |
| M(q) | composite system mass matrix | n×n | symmetric, configuration-dependent |
| N(q) | kinematic coupling matrix | nq×n | block-diagonal |
| f(t,q,u) | generalized forces | n-vector | applied + Coriolis |
| g(t,q) | constraint equations | m-vector | index-3 DAE; g=0 on manifold |
| G | constraint Jacobian ∂g/∂q | m×nq | - |
| λ | Lagrange multipliers | m-vector | unknown constraint forces |
| ${}^{X}R^{Y}$ | rotation of frame Y expressed in frame X | 3×3 ∈ SO(3) | left superscript = reference frame |
| G (frame) | ground / global reference frame | frame | fixed origin |
| A1 | attachment frame on residue 1's base | frame | constant orientation ${}^{B1}R^{A1}$ in body 1 |
| B1, B2 | body frames of residues 1, 2 | frame | centered on glycosidic nitrogen |
| O1, O2 | body origins of bases 1, 2 | position (nm) | force actually applied here, not at A1/B2 |
| r | translational distance between A1 and B2 | nm | r = |x_B2 - x_A1| |
| $\hat{r}$ | radial unit vector | 3-vector | (x_B2 - x_A1)/r |
| θ | Euler rotation angle of ${}^{A1}R^{B2}$ | rad | domain (-π, π); measures frame misalignment |
| $\hat{\theta}$ | angular (rotation-axis) unit direction | 3-vector | Euler axis direction |
| κ (kappa) | angular stiffness constant | per interaction type | typically positive |
| k | radial depth constant | kJ/mol scale, per interaction type | typically NEGATIVE for this potential; U=k·(...) at inflection |
| c | radial range / cutoff | nm, global | inflection point of g at r=c; harmonic for r<c, 1/r for r>=c |
| m | global force/energy scaling factor | dimensionless | "forceMultiplier" input; e.g. 10 or 20 |
| U(r,θ) | base-pairing potential energy | kJ/mol | eq:5 |
| g(r,k,c) | radial shape function | energy scale | eq:6; g(c)=k |
| g'(r,k,c) | dg/dr | eq:9 | - |
| $\vec{F}$ | force from potential (gradient) | force | eq:8 |
| $\vec{f}_{A1}, \vec{f}_{B2}$ | translational forces on the two bases | force | equal and opposite, eq:10 |
| $\vec{\tau}^{*}$ | adjusted torque (with moment-arm correction) | torque | eq:11 |

## Mobilizer (bond) types (Fig. 3 default RNA residue)

| mobilizer | meaning | RNA default |
|---|---|---|
| Torsion (pin) | fixed bond length & angle, free to rotate about axis; 1 DOF | backbone bonds, most ribose ring bonds, C2'-O2' bond, glycosidic N bond |
| Free | no restriction; 6 DOF | O4'-C1' ribose ring-closing bond (allows puckering) |
| Rigid | no freedom; 0 DOF | base bonds, bonds of single-coordinated atoms (H, some phosphate O) |

## Conventions / gotchas

- Rigid body has ZERO inherent DOF in Simbody; DOFs exist only where a mobilizer grants them (opposite of Cartesian MD where everything moves unless constrained).
- Sign convention: κ typically > 0, k typically < 0 for the base-pairing potential.
- Radial shape function g is continuous at r=c with value k.
- Force is applied at the body origin O, not the attachment/body frame A1/B2; eq:11 adds the (x - x_O) × f moment correction.
