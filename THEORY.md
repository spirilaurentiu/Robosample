# Robosample -- Theory of Operation

Robosample performs **Gibbs sampling coupled with Hamiltonian Monte Carlo (HMC)**. A configuration
is updated by a sweep of *Gibbs blocks* (executed in a chosen order -- systematic or random; the
order is **not** required to be fixed); each block proposes a new configuration by integrating
molecular dynamics over a chosen subset of coordinates and accepts or rejects it with a
Metropolis-Hastings (MH) test on a Hamiltonian function. The coordinates of a block may be
**Cartesian** (ordinary MD) or, more powerfully, a **reduced set of generalized (internal)
coordinates**. Torsional dynamics -- where bonds and angles are frozen and only torsions move -- is
the most efficient and currently the only fully implemented internal-coordinate case, but it is a
**special case** of the generalized-coordinate formalism used throughout; nothing below is specific
to torsions unless stated.

By default the target is the canonical (constant-temperature) Boltzmann distribution, but the
framework is **not restricted to NVT** (Section 2). Energies are evaluated on the full atomistic
model through **OpenMM** using an **AMBER or CHARMM** force field; the solvent treatment is
configurable -- **vacuum, implicit (e.g. GBSA-OBC2, the current default), or explicit** -- and is
not mandated by the theory (Section 8).

Notation convention: "~" means "distributed as"; "proportional to" is written out; superscripts use
"^" (e.g. M^-1, p^T, sin^2), and "1/2" denotes one half.

---

## 1. Notation

| Symbol | Meaning |
|---|---|
| q | generalized coordinates **mobile in the current block** -- the block's *restricted space*. In torsional dynamics: a subset of torsions, plus the 6 external DOF of any free root body. In general: any subset of internal/Cartesian DOF |
| p | generalized momenta conjugate to the mobile q |
| \mathcal{M}(q) | configuration-dependent generalized mass metric tensor (mass matrix) of the block's reduced multibody tree; the Riemannian metric on the restricted space (Section 6) |
| D_b | per-body articulated factor (the dof x dof block produced by the articulated-body algorithm at body b); ln det M(q) = sum_b ln det(D_b), computed in O(n) without forming M |
| M_{3N} | constant 3N x 3N Cartesian mass matrix (diagonal atomic masses); the metric used in Cartesian blocks. The *effective* reference in the Fixman term for torsional dynamics -- see M_B for the exact one |
| M_B | full 3N x 3N BAT mass-metric (M_B = J_{BAT}^T M_{3N} J_{BAT}); the **exact** Fixman reference (Jain et al. 2013). det M_B = \|J_{BAT}\|^2 * det M_3N, so M_B = M_3N up to the configuration-dependent BAT volume element \|J_{BAT}\|^2 (Section 6) |
| W_BAT(q) | internal BAT-volume compensating potential for any world that mobilizes a *hard* internal coordinate (Ball: angles; Cyl: bonds). W_BAT = -R*T* [ sum_{mobile angles k} ln sin theta_k + sum_{mobile bonds i} ln r_i^2 ]; identically zero (constant, cancels) under torsional dynamics (Section 6) |
| V(q) | potential energy (force field + chosen solvent model), evaluated on the full atomistic model through OpenMM |
| \mathcal{H}(q,p) | the **Hamiltonian function** of the block, H = V + K + F + J (+ W_BAT for Ball/Cyl worlds) (Section 6) |
| K(p,q) | kinetic energy, (1/2) p^T M(q)^-1 p (the Riemannian cometric form) |
| F(q) | Fixman compensating potential (Section 6) |
| J(q) | external-rotation Jacobian for free root bodies (Section 6) |
| Z-matrix | the internal-coordinate construction that places each atom by a bond length, bond angle, and dihedral relative to three previously placed atoms; it defines the BAT coordinates and the kinematic tree (Section 3) |
| BAT | bond-angle-torsion coordinates; the (b, theta, tau) decomposition induced by the Z-matrix |
| T, R | temperature; molar gas constant R = 8.3144626e-3 kJ/(mol*K). Energies are per mole (kJ/mol); beta = 1/(R*T) |
| sigma(q), d0 | ring-closure constraint sigma(q) = \|r_AB\|^2 - d0^2 = 0, with target distance d0 (the carried-over closure distance) |
| G(q) | constraint Jacobian d(sigma)/dq of the holonomic ring-closure constraints sigma(q) = 0. **Distinct** from the articulated gain Ga of Section 3.5 |
| gamma2_b | polar ("pitch") Euler angle of free root body b, extracted from its orientation quaternion |
| n, nu | n = number of bodies in the block's tree (the O(n) cost); nu = total mobilities (length of u). Distinct from M_3N's N = atom count |
| Phi, H, P, D, Ga, Mk_G | spatial-operator quantities of the multibody engine: rigid shift Phi, hinge map H, articulated body inertia P, hinge inertia D = H^T P H, articulated gain Ga = P H D^-1, body spatial inertia Mk_G (all defined in Section 3.5) |
| \pi | target distribution; Z is the configurational partition function |

---

## 2. Target distribution

By default the target is the canonical Boltzmann distribution

```
\pi(q) = Z^-1 * exp( -beta * V(q) )
```

(constant N, T), with pi defined on the flat Cartesian configuration measure dx; the
generalized-coordinate machinery of Sections 3 and 6 reproduces this same target, the
coordinate-change Jacobians absorbed into F and J. With implicit solvent there is no simulation box,
so there is no volume or PV term;
"NVT" here means constant-temperature canonical sampling. This is the **implemented default, not a
limitation**: the machinery samples any target of the form \pi(q) proportional to exp(-beta*phi(q))
for a configurational weight phi, so other ensembles or biased/tempered targets are admissible by
substituting the corresponding weight into the acceptance Hamiltonian (Section 6). Explicit-solvent
constant-pressure sampling, for instance, would add the appropriate PV/box terms; nothing in the
Gibbs/HMC construction (Sections 4-5, 13) assumes the canonical form.

Generalized-coordinate (e.g. torsional) sampling is designed so that the marginal distribution of
configurations matches the marginal of the **fully-flexible Cartesian** Boltzmann distribution
restricted to the sampled subspace; the Fixman term (Section 6) removes the mass-metric artifact
that would otherwise distort this marginal, and the mixed Gibbs scan (Section 13) supplies the
relaxation of the frozen coordinates.

---

## 3. Internal-coordinate representation (robot model)

Each molecule is represented as a **robot** (kinematic tree of rigid bodies connected by joints, or
*mobilizers*). An **O(n) articulated-body recursion** -- Featherstone's articulated-body algorithm,
in spatial-operator form -- propagates positions, velocities, accelerations, and forces up and down
the tree without forming or inverting the mass matrix M (Section 3.5). Per-atom Cartesian forces
(-grad V) computed by OpenMM are reduced to per-body spatial forces (net force + torque about the
body origin) for the articulated-body solver.

*Implementation note.* The recursion is realized in a custom, data-oriented (structure-of-arrays)
multibody engine; the global M, M^-1, sqrt(M), and det M are never assembled, only applied as O(n)
operators (Section 3.5).

### 3.1 Z-matrix, BAT coordinates, and the partition function

A molecule's internal coordinates are defined by a **Z-matrix**: each atom is placed by a bond
length r, a bond angle theta, and a dihedral (torsion) tau relative to three previously placed
atoms. The Z-matrix therefore induces the **bond-angle-torsion (BAT)** decomposition q_int =
(b, theta, tau) and, equivalently, the kinematic tree -- each joint of the tree corresponds to one
internal coordinate. Torsional dynamics is the case in which b and theta are frozen (their joints
welded) and only the tau (a subset of the Z-matrix dihedrals) remain mobile.

The BAT coordinates connect the program's reduced sampling to the Cartesian partition function. The
configurational partition function is Z = integral exp(-beta*V) dx over Cartesian x. The Cartesian
-> BAT change of variables has the Jacobian

```
dx = |J_{BAT}| * db * d\theta * d\tau * de,     |J_BAT| = ( product_i r_i^2 ) * ( product_j sin theta_j ) * const
```

(the standard Z-matrix volume element; de are external DOF). This Jacobian is **independent of the
torsions**, which is why torsional dynamics can target the correct conditional without an extra
positional Jacobian (Section 13). The remaining configuration-dependence that *does* matter -- the
mass-metric determinant det M(q) -- is handled by the Fixman term (Section 6). Sections 12-13 make
the link between this volume element, the Fixman correction, and the sampled marginal precise.

**The computed Z-matrix rows.** The builder returns four equal-length lists (z_i, z_j, z_k, z_l), one
row per atom, in local atom indices, with the sentinel -1 where a reference does not yet exist. Row r
describes the atom placed at step r and the already-built atoms it is measured against:

```
row r:  z_i = atom placed at step r
        z_j = bond reference      (bonded to z_i;  -1 for r = 0)
        z_k = angle reference     (bonded to z_j;  -1 for r < 2)
        z_l = dihedral reference  (bonded to z_k;  -1 for r < 3)
internal coordinate(s) of z_i:  bond |z_i z_j|,  angle (z_i,z_j,z_k),  torsion (z_i,z_j,z_k,z_l)
```

The four atoms of a full row are distinct and bonded in sequence z_i-z_j-z_k-z_l, and every reference
points **toward atoms already placed** (toward the root). The Z-matrix is therefore a strict build
order: an atom is positioned only after the three atoms it is defined against. This is the
"references toward the root" convention of the BAT construction (Chang, Potter & Gilson 2003; Hikiri,
Yoshidome & Ikeguchi 2016) and of the AlGDock/MDAnalysis BAT implementation it derives from.

**Root triplet (the first three rows).** The first three atoms carry sentinels and fix the molecular
frame rather than torsions: row 0 (z_i only) is the translation atom; row 1 adds one bond (z_j), the
r01 distance; row 2 adds a second bond and the first angle (z_j, z_k), the r12 distance and the
(z_i,z_j,z_k) angle. In the reference BAT these three atoms also carry the six external DOF -- the
first atom's Cartesian position plus the axis-angle rotation (polar, azimuthal, and a third angle);
in Robosample those six are the root body's external joint DOF (Section 3.3), with J(q) restoring
their Haar measure (Section 6). The builder picks the triplet by mass: a root atom, then its heaviest
non-terminal neighbour, then that neighbour's heaviest non-terminal neighbour.

**Traversal, tie-breaking, and the non-terminal rule.** Atoms are placed by a mass-prioritized walk
outward from the root -- heaviest-first, ties broken by ascending atom index for reproducibility. The
angle reference z_k must be **non-terminal** (degree > 1 in the full bond graph): a leaf cannot anchor
a stable angle/torsion chain, so a terminal atom is never used as z_k. z_l must differ from z_j. When
a row's references are not yet placed it is deferred and retried, which is what lets a single outward
pass yield a consistent order without back-tracking.

**Acyclic tree (ring closures removed).** The walk runs on the molecular bond graph with ring-closing
bonds deleted -- the same acyclic tree used for the rigid-body forest (Section 3.2) and the frame
build (Section 3.6). No Z-matrix reference is ever taken across a ring-closing bond; those bonds
re-enter only as explicit holonomic constraints (Section 7). Non-terminality, however, is judged on
the *full* bond graph, so a ring atom is correctly treated as non-terminal even after its closing
bond is removed.

**Atoms sharing a central bond; proper vs improper torsions.** Several atoms may carry the same
(z_j, z_k) pair -- the same central bond -- and hence reference the same triple. Chang's example is
two substituents on one bond: they share their three reference atoms and rotate as a unit, so their
torsions are strongly correlated. The first torsion about a given central bond is the *proper*
torsion; the others may be re-expressed as *improper* torsions (each one's difference from the
proper), which separates the shared rotation from the in-place distortion. Hikiri et al. show this
removes the three-body correlation of, for example, a methyl group's three hydrogens (one proper
torsion for the spin, two impropers for the distortion), and that the proper torsions are the
markedly anharmonic, multimodal coordinates (trans/gauche), whereas bond lengths, bond angles, and
improper torsions stay close to harmonic. The Z-matrix rows expose this grouping through the shared
(z_j, z_k); the improper *differencing* is a value operation applied downstream and is not part of the
row indices.

**Completeness and the volume element.** For N atoms the rows define exactly 3N - 6 internal
coordinates, plus the 6 external DOF carried by the root triplet -- a complete, non-redundant BAT set
(3N Cartesian <-> 3N BAT). The associated volume element is ( product_i r_i^2 ) * ( product_j sin
theta_j ) times a constant (Chang eq. 5), independent of the torsions, which is what lets a torsional
world target the correct conditional without a positional Jacobian, det M(q) being carried by the
Fixman term (Section 6).

**Robosample's builder vs the reference BAT, and vs the C++ frame kernel.** Robosample's Z-matrix
generalizes the AlGDock/MDAnalysis BAT in two indexing-relevant ways: the root is an arbitrary
caller-chosen atom (the reference forces the heaviest *terminal* atom), and the walk runs on the
acyclic graph with ring bonds removed (the reference walks the *full* bond graph, so on a cycle it can
reference across a ring bond). Robosample also orders the walk heaviest-first throughout, whereas the
reference grows its torsion list lightest-first; this changes the row order and, where an atom admits
several references, which neighbour becomes z_k or z_l -- but not the row schema, the four-atom
bonded-in-sequence rule, or the toward-root direction. Distinct again is the C++ frame kernel of
Section 3.6: it shares this acyclic tree (same z_j / z_k backbone) but takes its *local dihedral
reference* from the atom's first child (toward the leaves), so its fourth atom need not equal z_l. The
three trees agree on topology and on the bond and angle references; they differ only in the dihedral
reference atom and in row order.

### 3.2 Freezing by rigid-body lumping (no zero-velocity holds, no Schur complement)

A block that freezes a coordinate does **not** hold that DOF at zero velocity. Instead, rigidly
connected atoms are lumped into composite rigid bodies (welds), so the frozen coordinates simply do
not exist in that block's reduced tree. Consequently each block has a **single reduced articulated
mass matrix M(q)**, and the same M(q) drives both momentum sampling and the Fixman determinant.
Because (i) momenta are fully resampled at the start of every block -- nothing is carried across
blocks to condition on -- and (ii) each block's M(q) is the genuine articulated mass of its own
reduced tree rather than a sub-block of a larger matrix, **no Schur complement arises**. (Frozen
DOF are removed structurally here; this is distinct from the holonomic ring-closure constraints of
Section 7, which are imposed as Lagrange-type constraints on a manifold.)

### 3.3 The ground frame and how robots are rooted

The robot model has a single **Ground** (inertial) frame. Every mobilized body has exactly one
parent and one inboard mobilizer, and the tree is rooted at Ground. Each molecule's root body
attaches to Ground through a **root mobilizer**:

- **Weld** -- the root is rigidly fixed to Ground; no external DOF, J = 0 (single-molecule
  torsional sampling in a fixed frame).
- **Free** -- the root floats with 6 external DOF (3 translation + 3 rotation); the orientation
  contributes the Jacobian J(q) of Section 6.
- other mobilizers (e.g. Cartesian/translational roots) are possible.

For **multiple molecules**, each molecular robot is rooted **independently to Ground** (a forest of
trees rooted at the common Ground), *not* chained to one another: robot 0's root attaches to Ground
via its root mobilizer, and robot 1's root likewise attaches to Ground -- not to the last body of
robot 0. Chaining molecules would impose an artificial inter-molecular kinematic link and is not
done unless a joint is deliberately defined between them. This independent-grounding topology is why
the per-root mobility (Weld/Free/...) is specified separately for each robot, and why J(q) is summed
over free root bodies in Section 6.

### 3.4 The mobilizer frame chain (Simbody convention)

Each mobilized body B is positioned by a chain of frames between Ground and B. Following the
Simbody convention (Simbody Theory Manual, Sec. 5), the frames are:

- **G** -- Ground (inertial) frame.
- **P** -- parent body frame (the inboard body's frame).
- **F** -- the inboard "fixed" mobilizer frame, rigidly attached to the parent P.
- **M** -- the outboard "moving" mobilizer frame, rigidly attached to the child body B.
- **B** -- the child body frame, anchored at the body's root (inboard) atom; every atom of the body
  is stored as a constant *station* r_a in B.

Two transforms are constant (built once at setup from the Z-matrix geometry) and one carries the
joint motion:

```
X_PF      : P -> F, the inboard frame's placement on the parent      (constant)
X_BM      : B -> M, the outboard frame's placement on the child body (constant)
X_FM(q)   : F -> M, the mobilizer transform -- the ONLY configuration-dependent factor (joint DOF)
```

The body's pose in Ground is the ordered composition

```
X_GB = X_GP * X_PF * X_FM(q) * X_MB,        with   X_MB = (X_BM)^-1
```

equivalently `X_GB = X_GF . X_FM . X_MB` with `X_GF = X_GP . X_PF` (F's pose in Ground) and X_MB the
constant pose of the body frame in M. Atom positions follow as `r_a(Ground) = X_GB . r_a`. Only
X_FM(q) changes as q moves; X_PF and X_BM are fixed. For a Pin (torsion) joint the rotation is about
the mobilizer Z axis, so the setup bakes a fixed bond->Z alignment (a -90 deg rotation about Y in
this engine) into X_PF/X_BM, placing Z along the rotated bond.

**Common pitfall (X_MB vs X_BM).** The last factor is X_MB = (X_BM)^-1, **not** X_BM. The stored
constant is X_BM (B->M; the Simbody manual: "moving frame M attached to B with constant transform
B X M"); the composition needs its inverse. Concretely, X_BM.p() is Mo's position expressed in B
(the vector Bo->Mo in B), whereas X_MB.p() is Bo's position expressed in M (the vector Mo->Bo in M).
Substituting X_BM (or X_BM.p()) into the composition -- or into the r_MB shift term of the velocity
hinge map H -- has the wrong frame **and** the wrong sign; it silently corrupts every body whose M
frame is offset from Bo (i.e. every non-root torsion body, offset = bond length). The position path
can still look plausible while velocities and the kinetic energy are wrong, so an OpenMM energy
check on a welded test will not catch it -- golden-test the hinge map H directly. A guide that
writes "X_GB = X_GP X_PF X_FM X_BM" has exactly this error: the final factor must be X_MB.

**Orientation kinematics (exact specification).** For joints with a quaternion-parameterized
rotation (Free, and Ball when ported), the four orientation coordinates are the quaternion
q = (qw, qx, qy, qz) representing R_FM, and the three rotational generalized speeds are the angular
velocity w_FM = (wx, wy, wz) of M relative to F **expressed in the parent F frame** (Simbody
convention). The rotation map, the quaternion derivative, its inverse, and the second derivative are
exactly Simbody's (Rotation.h: calcUnnormalizedNForQuaternion / NInv):

```
R_FM(q) = [ 1-2(qy^2+qz^2)    2(qx*qy - qw*qz)   2(qx*qz + qw*qy)
            2(qx*qy + qw*qz)  1-2(qx^2+qz^2)     2(qy*qz - qw*qx)
            2(qx*qz - qw*qy)  2(qy*qz + qw*qx)   1-2(qx^2+qy^2) ]      (q normalized)

qdot = N(q) * w_FM ,    N(q) = (1/2) * [ -qx  -qy  -qz       (4x3)
                                          qw   qz  -qy
                                         -qz   qw   qx
                                          qy  -qx   qw ]

w_FM = NInv(q) * qdot , NInv(q) = 2 * [ -qx   qw  -qz   qy    (3x4)
                                        -qy   qz   qw  -qx
                                        -qz  -qy   qx   qw ]

qddot = N(q) * wdot_FM + N(qdot) * w_FM      (N is linear and homogeneous in its first argument)
```

Component form of the derivative -- the form the integrator and the SHAKE q-update (Section 7) must
use:

```
qdot_w = (1/2) * ( -qx*wx - qy*wy - qz*wz )
qdot_x = (1/2) * (  qw*wx + qz*wy - qy*wz )
qdot_y = (1/2) * ( -qz*wx + qw*wy + qx*wz )
qdot_z = (1/2) * (  qy*wx - qx*wy + qw*wz )
```

This is the **parent-frame** map. The body-frame map `qdot = (1/2) q (x) (0, w)` negates the
off-diagonal (cross-coupling) terms of N -- giving `qdot_x = (1/2)(qw*wx - qz*wy + qy*wz)`, and
likewise for qdot_y, qdot_z -- and expects w expressed in M. It **must not** be used here, because
w_FM is expressed in F and R_FM(q) above is the standard map: pairing the standard R_FM with the
body-frame N integrates the orientation with the wrong-handed angular velocity. The error is
invisible on welded tests (Pin/Translation/Weld have qdot = u and never touch N) and on the position
path (which never touches the qdot map), so it must be checked by golden-testing N against Simbody
on a free body. The quaternion is renormalized after each q-update.

### 3.5 The multibody dynamics engine (articulated-body recursions)

Every kinematic and mass-matrix quantity is computed by O(n) recursions over the n bodies of the
tree -- **Featherstone's articulated-body algorithm**, in spatial-operator form. The global mass
matrix M, its inverse M^-1, its square root, and its determinant are **never assembled as dense
matrices**; each is *applied* as an operator by one or two sweeps. This linear cost is the reason
the internal-coordinate representation exists, and it is what makes per-block sampling scale to large
systems.

**Spatial quantities and operators.** Each body b carries, in the Ground frame, a spatial velocity
V_GB = (w ; v) (angular **over** linear), a spatial force F = (torque ; force) (same ordering), and
a spatial inertia Mk_G[b] (the body's rigid-body spatial inertia about its own origin Bo,
re-expressed in Ground; mass, first moment mass*c, and the inertia tensor). The spatial inner
product is the straight component pairing < V, F > = w . torque + v . force, so a generalized force
is tau = H^T F under that pairing. Four operators relate a body to its parent p:

- **Phi (rigid shift / transmission).** With r = p_PB the vector from the parent origin Po to the
  body origin Bo (in Ground), Phi is the 6x6 built from 3x3 blocks [[1, r_x],[0, 1]], r_x the
  cross-product matrix of r. ~Phi (its transpose) shifts velocities and accelerations **outward**
  (parent -> child); Phi shifts forces and inertias **inward** (child -> parent). It carries no
  mass -- pure geometry.
- **H (hinge / joint map).** Maps a body's dof generalized speeds to its cross-mobilizer spatial
  velocity, V_PB_G = H u. H is built in the F frame (Section 3.4) and rotated to Ground.
- **P (articulated body inertia).** The apparent spatial inertia a body exhibits as the free base
  of the subtree hanging below it. Unlike a rigid-body inertia, it is a full symmetric 6x6 (the
  apparent mass depends on the direction of push, and there is no single center of mass).
- **D, Ga.** The mobility-space (hinge) inertia D = H^T P H (dof x dof), its inverse D^-1, and the
  articulated **gain** Ga = P H D^-1. (The symbol Ga is deliberately distinct from the constraint
  Jacobian G of Sections 1/6/7 and the kinetic energy K of Sections 1/5/6 -- different objects.) The
  inertia the parent feels across the mobilizer is the articulated shift
  P^+ = P - Ga D Ga^T = P - P H D^-1 H^T P.

**The recursions -- each a single sweep, O(N).**

1. **Position -- outward (root -> leaves).** Compose the frame chain
   X_GB = X_GP X_PF X_FM(q) X_MB (Section 3.4); build H in Ground, the shift Phi, the body spatial
   inertia in Ground, and the per-atom stations r_a(Ground) = X_GB r_a.
2. **Velocity -- outward.** Shift the parent's spatial velocity and add the mobilizer term:

   ```
   V_GB[b] = (~Phi[b]) * V_GB[parent] + H[b] * u[b]
   ```

   The same sweep accumulates the two velocity-dependent bias terms the integrator needs. The
   gyroscopic force (Ground frame; I the unit inertia about Bo, c the COM offset, w = w_GB):

   ```
   b_gyro[b] = mass * ( w x (I w) ;  w x (w x c) )
   ```

   and the Coriolis bias acceleration, carried outward like the velocity (w_GP, v_GP the parent
   spatial velocity; w_PB_G, v_PB_G the cross-mobilizer velocity; w_FM the cross-mobilizer angular
   velocity in F; r_MB_F = R_FM * r_MB the Mo->Bo bond vector expressed in F):

   ```
   a_mob[b] = ( w_GP x w_PB_G ;
                w_GP x (v_GB - v_GP) + w_GP x v_PB_G + R_GF * ( w_FM x (w_FM x r_MB_F) ) )
   a_cor[b] = (~Phi[b]) * a_cor[parent] + a_mob[b]
   ```

   The centripetal term R_GF (w_FM x (w_FM x r_MB_F)) is nonzero for every body whose M frame is
   offset from Bo (every non-root Pin body, offset = bond length); dropping it silently pumps energy.
   It is the velocity-side companion of the X_MB pitfall of Section 3.4 -- r_MB = X_MB.p(), never
   X_BM.p().
3. **Articulated body inertia -- inward (leaves -> root).** Featherstone's backward pass:

   ```
   P[b]   = Mk_G[b] + sum_{children c} Phi[c] * P^+[c] * (~Phi[c])
   D[b]   = H[b]^T P[b] H[b],   Ga[b] = P[b] H[b] D[b]^-1,   P^+[b] = P[b] - Ga[b] D[b] Ga[b]^T
   ```

   A Weld (0 dof) contributes nothing through its (empty) mobility space, so P^+ = P.
4. **Forward dynamics (M^-1) -- inward then outward.** Given a generalized force f, the acceleration
   udot = M^-1 (f - f_bias) is produced by an inward residual pass (the per-body bias Z and the
   mobility residual eps = f - H^T Z) followed by an outward acceleration pass (A_GB and udot), using
   D^-1 and Ga. The bias force f_bias is assembled from b_gyro and a_cor of recursion 2 (the
   gyroscopic force plus the Coriolis acceleration mapped through the articulated inertia): with that
   f_bias this is the equation of motion solved each Verlet substep; with f_bias = 0 it is the pure
   operator a = M^-1 f used for the constraint solver and the kinetic energy.
5. **Force -- inward.** The transpose of forward kinematics: per-atom Cartesian forces are reduced to
   a generalized force by accumulating each body's spatial force and sweeping child -> parent through
   Phi^T, reading tau = H^T (transmitted spatial force). This is exactly the G^T assembly the
   constraint solver uses (Section 7.1).

**Mass-matrix operators -- all O(n), all matrix-free.** The same per-body factors {Phi, H, D, Ga}
yield, without ever forming M:

- **M^-1 f** -- forward dynamics (recursion 4 with f_bias = 0).
- **sqrt(M^-1) z** -- one outward sweep using the per-body D^-1/2 in place of D^-1. It maps a
  white-noise vector z to a generalized-speed vector with covariance M^-1, so the equipartition draw
  of Section 5.1 is realized as u = sqrt(R*T) * sqrt(M^-1) * z, giving u ~ N(0, R*T*M^-1) and hence
  p = M u ~ N(0, R*T*M) -- at O(N), with det M never appearing in the proposal.
- **ln det M = sum_b ln det(D_b)** -- the factorization makes the determinant the product of the
  per-body hinge determinants; this is the Fixman tree term of Section 6, read from the same
  D_b = H_b^T P_b H_b.
- **(1/2) u^T M u** -- the kinetic energy, from one application of M.

**The operator factorization.** In spatial-operator form the mass matrix factors through the same
local quantities,

```
M       = (I + H Phi Ga) D (I + H Phi Ga)^T
M^(1/2) = (I + H Phi Ga) D^(1/2)
v       = D^(1/2) (I + H Phi Ga)^T * theta_dot          (the whitening map)
```

with Ga the gain operator (the per-body Ga = P H D^-1 assembled across the tree) and v the whitened
velocity. Because each factor is local -- a body and its parent only -- M, M^-1, M^(1/2), and
ln det M all follow from one or two sweeps over the tree, which is the O(n) claim. The generalized
equipartition principle (Section 5.1) is then exact by construction: drawing the whitened v as unit
Gaussian noise and applying sqrt(M^-1) gives < u u^T > = R*T*M^-1, i.e. equipartition of kinetic
energy across the reduced coordinates despite the configuration-dependent metric.

### 3.6 Construction of the mobilizer frames from bond geometry

Section 3.4 gives the frame chain X_GB = X_GP X_PF X_FM(q) X_MB; this section specifies how its three
factors are actually built. X_FM(q) is the joint DOF (Section 3.4: Rot(q, Z) for a Pin, identity for
a Weld, Translate(q) for a Translation, quaternion+translation for a Free). The two **constant**
frames X_PF and X_BM are rebuilt once per configuration transfer (Section 9) from the actual Cartesian
geometry, in two stages.

**Stage 1 -- a frame per atom, from bond vectors.** Each atom a gets a Ground-frame pose frame[a]
and a cross-bond transform xpc[a], built from four atoms: the atom (self), its parent, its
grandparent, and a reference child (the four atoms that define the atom's internal coordinates).
With pSelf, pParent, pGparent, pRefChild their Cartesian positions:

```
R_G : X axis = unit(pSelf - pParent)                          // along the bond, parent -> self
      Y axis = component of unit(pGparent - pParent) perp to X   (Gram-Schmidt)
      Z axis = X x Y                                          // normal to the (self,parent,gparent) plane
G        = ( R_G , pParent )                                  // a frame anchored at the PARENT atom
theta    = dihedral(pGparent, pParent, pSelf, pRefChild)      // torsion about the parent->self bond
d        = |pSelf - pParent|                                  // bond length
B        = Rot(theta, X) * Translate(d along X) * Rot(180 deg, Y)
frame[a] = G * B * Rot(180 deg, X)
xpc[a]   = B
```

So the **bond vector is the local X axis**, the **grandparent fixes the perpendicular** (the
(self,parent,gparent) plane; its normal is Z), the **dihedral** rides on X, and the **bond length**
is the X-translation. The two 180-degree flips (about Y inside B, then about X) reproduce the
bond-center "about-face" convention of the original Molmodel build -- the parent's outboard
bond-center direction and the child's inboard one point in opposite senses along the bond -- so
frame[a] matches the original bond-center frame protocol to ~1e-14. (That equivalence is exactly what
the angle/chirality reconstruction cross-check verifies: rebuilding each atom's local bond-center
directions from bond angles + chirality and comparing against the position-derived frame localizes
any mismatch to the offending atom.) A **root** atom gets frame = (identity rotation, pSelf),
xpc = identity. An atom lacking a full (parent, grandparent) lineage falls back to X = bond
direction with xpc = Translate(d, X) * Rot(180 deg, Y) and no dihedral.

**Stage 2 -- mobilizer frames from the per-atom frames.** For body b with root (inboard) atom `root`,
and parent body p with root atom `proot`:

```
T_X_B = frame[root]                                          // the body frame B == its root atom's frame
if p is Ground:
    X_PF[b] = T_X_B                                          // F is the root-atom frame, in Ground
    X_BM[b] = identity                                       // M coincides with B
else:
    Proot_X_root = inverse(frame[proot]) * frame[root]       // pose of B's root atom in p's root-atom frame
    X_BM[b] = xpc[root] * X_to_Z                             // B -> M
    X_PF[b] = Proot_X_root * X_BM[b]                         // P -> F
where   X_to_Z = Rot(-90 deg, Y)                             // carries the bond (X) onto the mobilizer Z axis
```

`X_to_Z` is the only joint-convention-specific piece: a Pin rotates about the mobilizer **Z** axis,
but Stage 1 puts every bond on **X**, so the -90-degree rotation about Y carries X onto Z
(X -> +Z, Z -> -X). The M frame's Z axis therefore lies along the inboard bond, and X_FM(q) = Rot(q, Z)
turns the body about that bond. The reference dihedral is baked into xpc[root] (hence into X_BM and
X_PF), so at q = 0 the chain reproduces the **incoming** geometry exactly (Section 9), and the Pin
coordinate q is the torsion measured from that carried-over reference.

**The closure identity (the cheapest self-check).** With these definitions and X_FM(0) = identity,
the Section 3.4 composition must return the body frame:

```
X_GB(0) = X_GP * X_PF * X_FM(0) * X_MB
        = frame[proot] * (Proot_X_root * X_BM) * I * (X_BM)^-1
        = frame[proot] * Proot_X_root
        = frame[proot] * inverse(frame[proot]) * frame[root]
        = frame[root] = T_X_B.
```

Any change to the frame build must preserve this identity; comparing X_GB at q = 0 against frame[root]
for every body is the fastest regression test, and -- unlike an energy check -- it catches the
X_MB-vs-X_BM and bond->Z sign errors of Section 3.4 directly.

**Update frequency.** The three factors of the chain are recomputed on deliberately different
schedules, and that split is the point of the design:

| factor | recomputed | role |
|---|---|---|
| X_PF (P->F), X_BM (B->M, hence X_MB = ~X_BM) | once per Gibbs-block transfer (Section 9), plus once at world construction | the inboard bond's frozen geometry (length + angle + reference dihedral); held constant for the block's entire MD trajectory |
| X_FM(q) | every MD step | the live joint DOF; for a Pin, Rot(q, Z) with the current torsion q |
| X_GP (the parent body's X_GB) and X_GB itself | every MD step | the running composition; X_GP is constant only for a root body (parent = Ground) |

So X_PF and X_BM carry the per-block **frozen** geometry while X_FM, X_GP, and X_GB carry the
**per-step** motion. The same transfer that rebuilds X_PF and X_BM also resets q = 0, making the
incoming configuration the block's reference (Section 9); within a block only X_FM(q) and the
ancestor poses move, so X_GB(q) is the only thing the integrator recomputes each step.

**Relation to the Z-matrix (BAT tree).** The frame graph that Stage 1 consumes is filled from the
molecular bond graph with ring-closing bonds removed (parent = the bonded atom toward the root,
grandparent = its parent, refChild = the atom's first child). That acyclic tree is the same
BAT / Z-matrix internal-coordinate tree the topology builder produces as (z_i, z_j, z_k, z_l). The
correspondence is:

- tree topology -- each atom's parent (z_j) and grandparent (z_k) -- is identical to the Z-matrix;
- the bond length d = |self - parent| and the bond angle (the grandparent fixes the frame's
  perpendicular) are the Z-matrix bond/angle internal coordinates, and they are frozen into X_PF and
  X_BM (rebuilt each transfer from the actual Cartesian positions, Section 9);
- torsions are the live coordinates: a flexible dihedral about a tree bond becomes a Pin joint, so a
  Z-matrix torsion is q = X_FM(q). Its reference value is baked into X_BM via xpc[root] (which carries
  Rot(theta, X)), and X_FM(q) = Rot(q, Z) turns about that bond.

Two qualifications are load-bearing and must not be relaxed into an identity. First, the frame kernel
derives this tree and all geometry from the Cartesian positions; it shares the *tree* with the
separate Z-matrix builder, not the code, and it uses the atom's first child as the local dihedral
reference (pointing away from the root), whereas a textbook Z-matrix z_l references toward the root --
both describe the torsion about the same bond, but with a different fourth atom. Second, the frozen
bond lengths and angles are not hard-constrained to idealized values: they sit at whatever the
incoming geometry holds and are merely held fixed for the block (Gibbs conditioning, Section 9), never
a permanent constraint.

**Inboard and outboard atoms.** For body b with root atom root = bodyRootAtom[b] (the code's "inboard
(root) atom of this body"), the inboard bond is the tree edge (par[root], root):

- the **inboard** atom A = par[root] lives in the parent body; the F frame (X_PF) is fixed to the
  parent at this joint;
- the **outboard** atom B = root is the child body's root atom; the M frame (X_BM) and the body frame
  B are anchored here.

Which *physical* atom is inboard versus outboard is set per molecule by the Z-matrix root (a
mass-prioritized walk from a chosen root atom), not by any fixed convention such as N-before-CA: for a
given molecule it is determined, but it must be read off that molecule's tree rather than assumed.

---

## 4. The Gibbs scan

Blocks are executed in a chosen order -- a deterministic **systematic scan** or a randomized
**random scan**. The order is a free choice, not a constraint of the method: each block individually
leaves pi invariant (Sections 5, 13), and a composition of pi-invariant kernels leaves pi invariant
regardless of order. A systematic scan is generally **not reversible** (its time-reversal is the
reversed block order); a random scan that picks a block from a fixed distribution **is** reversible.
Either way pi is preserved. **Stationarity (pi*K = pi), not reversibility, is the requirement**
(Section 13).

Although torsional dynamics is the most efficient way to exploit the internal-coordinate
representation -- and the reason the representation exists -- it is **not** mandatory; a scan may mix
Cartesian, torsional, and other generalized-coordinate blocks freely.

- **Overlap.** The restricted spaces of different blocks may overlap. Redundant updates do not break
  stationarity, provided each block individually leaves pi invariant.
- **Ergodicity (precise statement).** Across one full round (all blocks), every mobile coordinate in
  the block set has been *offered* a move. Ergodicity is over the **sampled (restricted) subspace**
  and requires both that the blocks collectively cover that subspace and that each block achieves
  nonzero acceptance. Being offered a move is **not** the same as accepting one: a block stuck near
  zero acceptance contributes nothing even though aggregate statistics look healthy. Coordinates
  that are frozen in **every** block (e.g. bonds and angles in torsion-only operation) are not
  sampled; the sampled distribution is then the (Fixman-corrected) marginal over the union of the
  blocks' restricted spaces (Section 13.5).
- **Per-block diagnostics.** Acceptance rates are reported **per block** specifically to detect an
  under-sampled subspace, which would otherwise be masked by healthy global energy and acceptance
  statistics.

---

## 5. The HMC move inside a block

### 5.1 Momentum resampling and the internal-coordinate equipartition principle

At the start of each block the momenta are **fully** resampled (not partial):

```
p ~ N(0, R*T*M(q)),     equivalently    velocities u ~ N(0, R*T*M(q)^-1).
```

This is a Gibbs step on the momenta (always accepted). It realizes the **internal-coordinate
equipartition principle** of Jain, Park & Vaidehi (2012): in Cartesian coordinates every DOF carries
(1/2)*R*T of kinetic energy, but in generalized coordinates with a configuration-dependent metric
M(q) this does not hold per coordinate. The *generalized* equipartition theorem,
`< p_i (dH/dp_j) > = R*T delta_ij`, still holds, which forces the momentum covariance to be
`< p p^T > = R*T*M(q)`. Drawing p ~ N(0, R*T*M(q)) is exactly this distribution, so temperature is
assigned consistently despite the non-constant metric. The draw is matrix-free via the articulated
sqrt(M(q)^-1) operator (drawing velocities; momenta follow as p = M(q)*u), so det M(q) never appears
explicitly in the proposal; its configuration dependence is bookkept by the Fixman term in the
acceptance (Section 6).

In a **ring-closure block** the freshly drawn momenta are then RATTLE-projected onto the velocity
constraint surface G M^-1 p = 0 (Section 7) before integration begins, so the initial velocity is
consistent with the closed ring. That projection is exactly what injects the det(G M^-1 G^T)^(-1/2)
factor into the configurational marginal -- the factor the loop-closure term of F cancels (Section 6).
For acyclic molecules there is no projection and no such factor.

### 5.2 Propagation under the bare potential V

Dynamics are integrated with the **fixed-step Verlet integrator** (internal and Cartesian alike),
the velocity-Verlet-class method of Simbody's `VerletIntegrator`, using forces derived from the
**bare potential V only**. Concretely, "under V only" means: the integrator's forces are the
per-atom Cartesian forces -grad V from OpenMM, reduced to per-body spatial forces for the
articulated solver (Section 3); the gradients of the Fixman potential F and the Jacobian J (the
"torques") are **not** added to these forces. So V is what generates the trajectory; F and J enter
only the acceptance (Section 5.4). After each step, constraints (if any) are projected onto the
manifold (Section 7).

**The integrator is a second-order, semi-explicit predictor-corrector, not an explicit
kick-drift-kick leapfrog.** This distinction is load-bearing and is the source of the integrator's
actual conservation properties (Section 5.5): the position is advanced by an explicit second-order
Taylor step using the start-of-step acceleration, and the velocity is advanced by an **implicit
trapezoidal corrector** solved by functional iteration. The two coincide with textbook velocity
Verlet only for a separable, constant-mass Hamiltonian (the Cartesian case); for the
configuration-dependent metric M(q) and velocity-dependent (Coriolis/gyroscopic) forces of internal
coordinates they do not, which is exactly where the symplecticity caveat of Section 5.5 applies. The
integrator is run with a **fixed step** (Simbody: `setFixedStepSize`); in that mode Simbody's Verlet
is symplectic for separable systems and the step is **taken unconditionally** -- there is no
step-size adaptation, so a non-converged corrector does not shrink dt or reject the step (Section
5.5). It is the trajectory-level Metropolis test (Section 5.3), not per-step control, that supplies
correctness.

The per-step sequence (one fixed-step Verlet step of size dt; the SHAKE/RATTLE projections are no-ops
on acyclic molecules). qdot0 = N(q0) u0 for quaternion DOF and u0 otherwise; qdotdot0 is the
corresponding q-acceleration; udot = M^-1 (f - f_bias) is one O(n) forward-dynamics sweep
(Section 3.5, recursion 4), with f_bias carrying the gyroscopic and Coriolis terms:

```
1. position (Taylor, 2nd order)   q  <- q0 + dt*qdot0 + (dt^2/2)*qdotdot0   (all coords, incl. quaternion)
2. quaternion rotation (exp-map)   q  <- expmap(wHalf, dt) (x) q0 ; wHalf = u0 + (dt/2)*udot0   then q <- q/|q|
3. SHAKE                           project q so that sigma(q) = 0            (Section 7.1), refresh kinematics
4. velocity predictor (fwd Euler)  u  <- u0 + dt*udot0                      (forces from -grad V at q1)
5. velocity corrector (trapezoid)  repeat (<=10): u <- u0 + (dt/2)*(udot0 + udot1); recompute udot1
                                   until ||du||/||u|| <= tol (= min(1e-4, 0.1*accuracy))
6. RATTLE                          project u so that G M^-1 p = 0           (Section 7.1)
```

The corrector is **plain functional iteration**, which has a limited radius of convergence; if it
stops contracting it is abandoned and **the step is taken with the best u so far** (Simbody's
fixed-step behaviour -- it cannot shrink dt). A persistently non-converging corrector is the signal
that dt is too large for the block; the resulting large dH is rejected by the MH test, and only when
that test is disabled (e.g. an equilibration/AlwaysAccept phase) does the unconverged step manifest
as visible energy drift (Section 5.5). The orientation of a quaternion (Free/Ball) body is **not**
advanced by the linear Taylor formula of step 1; it is advanced by an **exact exponential-map
rotation** (step 2): the increment quaternion expmap(wHalf, dt) = (cos(theta), sin(theta) * wHalf/|wHalf|)
with theta = (1/2)|wHalf| dt and the midpoint angular velocity wHalf = u0 + (dt/2) udot0, applied as
a left Hamilton product onto q0, followed by renormalization that mops up only rounding. This is a
**deliberate divergence from Simbody's linear-Taylor-plus-renormalize quaternion update** (Section
3.4), and it is the propagator of record. Three properties make it the correct choice here: it keeps
|q| = 1 by construction (the renormalization is O(1e-16), not a real projection); it is reversible
and reduces to the linear update as dt -> 0, so it changes only the proposal, never the target
(Section 5.4); and -- decisively -- it advances the orientation **purely from the angular velocity**,
bypassing the quaternion second derivative qddot. The linear-Taylor update instead leans on qddot;
with the engine's reimplemented Free-joint quaternion kinematics that path **injects kinetic energy**
through the root body and, via the articulated recursion V_GB[b] = Phi V_GB[parent] + H u, through
the entire tree -- an observed one-sided KE pump of order 1e3 kJ/mol per trajectory even at dt = 1 fs.
The exp-map is robust to that latent inconsistency; the inconsistency itself (a not-yet-byte-faithful
N / Ndot / qddot for the Free joint relative to Simbody) is the open item flagged for golden-testing
in Section 3.4. The converged map is time-reversible and volume-preserving on the constraint manifold
(Section 5.5), which is what the MH test of Section 5.3 requires.

**Integration order and error estimate (the registered constants, and what they mean).** Simbody
registers this integrator as `AbstractIntegratorRep(handle, sys, 2, 3, "Verlet", true)`, i.e. with
**method min-order 2, method max-order 3, error-control enabled**, and inside the step it sets the
local error estimate's order to **errOrder = 3**. The class is documented as "a third order,
semi-explicit integrator" whose "velocities ... are only accurate to lower order." These three
numbers (2, 3, 3) are not interchangeable; precisely:

- *Definition.* A one-step method of **order p** has **global** error O(dt^p) over a fixed time
  interval and **local** (single-step) truncation error O(dt^(p+1)). The two differ by one power of
  dt because local errors accumulate over ~T/dt steps.
- **Method order = 2 (the governing number).** `getMethodMinOrder() = 2`. Verlet is a **second-order
  method**: trajectories and conserved quantities are accurate to **O(dt^2)** globally. This is *the*
  number that governs the science -- it is the h^2 in the energy bound (Section 5.5) and in the
  step-size dependence of every computed average (Section 5.6; Davidchack 2014). The practical rule
  follows directly: **halving dt quarters the discretization error and the energy fluctuation.** The
  velocities are specifically the second-order (lower-order) quantity -- Simbody's "only accurate to
  lower order" -- which is why velocity-dependent (Coriolis) forces degrade accuracy fastest in the
  internal-coordinate case (Section 5.5).
- **Max-order = 3 (the most accurate quantity, and the error-estimate order).** `getMethodMaxOrder()
  = 3` and `errOrder = 3`. The **position** Taylor step q1 = q0 + dt*qdot0 + (dt^2/2)*qdotdot0 has a
  **local** truncation error O(dt^3), and the integrator forms a corresponding **3rd-order local
  error estimate** (a per-step quantity, *not* the global trajectory order). "Third order" in the
  Simbody label refers to this local/error-estimate order, **not** to the global accuracy.
- **What errOrder = 3 is used for, and why it is inert here.** In variable-step mode the controller
  picks the next step from the error estimate as dt_new proportional to dt * (accuracy/err)^(1/errOrder)
  = (...)^(1/3). **Robosample runs the integrator in fixed-step mode (Section 5.5), so this estimate
  is computed but never acted upon** -- the step size cannot change. errOrder therefore governs
  nothing in the implemented sampler; it is documented here for fidelity to the reference and for the
  event a variable-step mode is ever enabled (which, per Section 5.5, would also forfeit
  symplecticity).

In one sentence for the reader who wants only the operative fact: **treat the integrator as second
order (errors and energy fluctuation scale as dt^2); the "third order" is a local/error-estimate
property that does not govern the sampled results and is unused in fixed-step operation.**

### 5.3 Acceptance

The proposal is accepted or rejected by an MH test on the full Hamiltonian function H (Section 6),
using the **exact** dH between the pre- and post-trajectory states:

```
accept with probability min(1, exp(-beta * dH)),    dH = H_new - H_old.
```

### 5.4 Guidance vs. acceptance Hamiltonian (Duane et al. 1987)

HMC (Duane, Kennedy, Pendleton & Roweth 1987) admits a **distinct guidance and acceptance
Hamiltonian**. The *guidance* Hamiltonian generates the trajectory (here: V's forces); the
*acceptance* Hamiltonian is the one whose Boltzmann distribution is actually sampled, used in the MH
test (here: the full H = V + K + F + J). Duane et al.'s key observation is that the trajectory
generator need **not** equal the target: as long as the proposal map is reversible and
volume-preserving and the MH test uses the *exact* acceptance Hamiltonian, the move samples the
acceptance distribution exactly, whatever the guidance was. Two consequences used here:

- **The Fixman torque is unnecessary for correctness.** Because F enters only the acceptance, the
  integrator need not compute grad F. (Robosample nonetheless implements the Fixman torque as a
  per-world option; including it changes only the proposal trajectory, never the target. In a
  pure-MD setting with no MH test, by contrast, the Fixman torque in the forces *is* required to
  reach the Boltzmann distribution -- hence its availability.)
- **An efficiency-oriented guidance potential is admissible.** A cheaper or smoother potential
  (e.g. an internal-coordinate force field; Section 12) may be used to *guide* proposals while the
  exact atomistic H is retained for acceptance, with no bias.

### 5.5 Properties of the fixed-step Verlet integrator

The validity of the HMC move rests on three properties of the (constrained) proposal map -- time
reversibility, volume preservation, and (for the energy bound) symplecticity. **These hold to
different degrees for the separable Cartesian case and the non-separable internal-coordinate case,
and the distinction is essential to read correctly.** What HMC strictly requires for *correctness*
is reversibility + volume preservation + exact dH (Section 5.4); symplecticity is what additionally
*bounds* the energy error and hence protects the acceptance rate.

- **Separable (Cartesian) blocks -- the strong case.** When M = M_3N is constant and the forces are
  velocity-independent, the implicit corrector converges in a single pass and the scheme reduces
  **exactly to Stoermer/velocity Verlet**, which is symplectic. Simbody's documentation states this
  for the fixed-step mode used here: with fixed time steps the integrator is symplectic and conserves
  energy extremely well. For such blocks the backward-error result applies: a symplectic integrator
  exactly conserves a nearby **shadow Hamiltonian** H~ = H + O(dt^2), so

  ```
  |H(t) - H(0)| <= C * dt^2
  ```

  is an O(dt^2) **bounded oscillation, not a secular drift**, over times exponentially long in 1/dt
  (smooth V, stable dt; near the stability limit or with stiff constraints the practical guarantee is
  the weaker O(dt^2) bounded error). The constant C grows with the stiffness of V.

- **Internal-coordinate blocks -- the qualified case.** Here the Hamiltonian is **non-separable**:
  K = (1/2) p^T M(q)^-1 p with the metric M depending on q, and the bias forces are
  velocity-dependent (Coriolis/gyroscopic). The corrector is then a genuine implicit solve and the
  method is the **implicit-trapezoid Verlet**, *not* textbook velocity Verlet. Its properties:
  - It is **time-reversible and volume-preserving when (and only when) the corrector converges** to
    the trapezoidal solution -- this is what HMC needs, and it holds at any dt small enough for the
    functional iteration to converge.
  - It is **not strictly symplectic** for the configuration-dependent metric, so the exact
    shadow-Hamiltonian guarantee above does **not** apply. Energy conservation is the weaker property
    that conservative non-symplectic, time-reversible methods enjoy: bounded for small dt, but a
    **secular energy drift is permitted** and grows with dt and with the magnitude of the
    velocity-dependent forces. (Simbody notes that with velocity-dependent forces the reported
    velocities are only accurate to lower order, reducing overall accuracy; Davidchack 2014 observes
    precisely this upward energy drift in rigid-body NVE integration once the step exceeds roughly
    60-70% of the stability threshold -- about 4.5 fs for rigid TIP4P water.) **A one-sided dH growth
    at large dt is therefore the expected behaviour of a correct implementation, not a defect**; it is
    the discretization error the MH test is designed to absorb.
  - **The MH test absorbs this error exactly in the sampled distribution** (Section 5.4): a drifting
    trajectory produces a large dH and is rejected, so the cost of a too-large dt is acceptance rate,
    never bias -- *provided the MH test is active*. When acceptance is disabled (an
    equilibration/AlwaysAccept phase), the bare integrator drift is what is observed, and it must not
    be given a physical interpretation.

- **Corrector convergence is a precondition, not an afterthought.** A non-converged corrector step is
  not the exact trapezoidal map and is therefore **not reversible**, which would invalidate the MH
  proposal ratio for that step. Simbody's fixed-step driver nonetheless takes the step (it cannot
  shrink dt); the implementation must therefore be run at a dt small enough that the corrector
  converges in practice, with the MH test rejecting the occasional bad trajectory. A persistently
  non-converging corrector signals that dt is too large for the block and should be reduced.

- **Manifold preservation (SHAKE/RATTLE).** When holonomic ring-closure constraints are active
  (Section 7), SHAKE projects positions and RATTLE projects velocities back onto the constraint
  manifold each step, and the constrained map remains reversible and measure-preserving **on that
  manifold** (and symplectic there in the separable case). A holonomic constraint sigma(q) = 0
  depends on positions q **only** (a single algebraic relation per constraint, not on velocities), so
  each independent constraint removes exactly one dimension from the configuration manifold
  (3N -> 3N - m for m constraints) and the conjugate momentum is restricted correspondingly. The
  projection costs O(m*n) per Newton iteration plus an O(m^3) solve (Section 7.1), small for the few
  ring closures of a macrocycle.

### 5.6 Time step is a property of the proposal, not the target (multiscale HMC)

A larger dt is a coarser Verlet map. Provided the map is reversible and volume-preserving and the
**exact dH** is used in the MH test, the move is exactly pi-stationary at any dt; only the energy
bound of Section 5.5 (hence the acceptance rate) changes. Assigning small dt to stiff blocks and
large dt to soft blocks therefore introduces **no bias** -- only a change in efficiency. This is the
sense in which Robosample is a multiscale HMC sampler.

The one precondition (Section 5.5) is that the implicit velocity corrector **converges** at the
chosen dt, since reversibility holds only for the converged trapezoidal map; an un-converged step is
not reversible and is the signal that dt is too large for that block. In practice this bounds the
usable dt from above for stiff blocks, but within that range the no-bias property is exact and
independent of dt.

### 5.7 Coordinate-space geometry, corrector convergence, and reversibility

This section makes precise three things the integrator's correctness rests on: why the propagator has
**two distinct paths** (one for flat coordinates, one for the orientation manifold), how a **too-large
dt is detected for free** from the corrector, and why the safe dt -- and hence reversibility -- is
**configuration dependent**, so that a single startup check cannot certify a whole run.

**Two propagation paths, set by the shape of each coordinate's space.** A mobilizer's generalized
coordinates live in one of two kinds of space, and the integrator advances each kind differently
because a straight-line (Taylor) step is only valid in a flat space.

- *Flat path (linear Taylor; qdot = u).* Pin torsions, Translation, and the Free joint's translation
  block. Their coordinates lie in a flat (zero-curvature) space: Translation and the Free-translation
  block in R^3; a Pin torsion on the circle S^1, which has **zero intrinsic curvature** (a circle is a
  line made periodic, locally indistinguishable from R) and is integrated as an unwrapped real. The
  velocity-to-coordinate map is the identity (qdot = u). A second-order Taylor step
  q1 = q0 + dt*qdot0 + (dt^2/2)*qddot0 lands on another valid point of the same flat space with local
  error O(dt^3); nothing to project. This path is genuine velocity Verlet.
- *Curved path (exponential map).* The orientation of a Free (or, when ported, Ball) body, stored as a
  unit quaternion living on **S^3, the unit 3-sphere in R^4** -- a compact manifold of **constant
  positive curvature** that **double-covers the rotation group SO(3)** (q and -q are the same physical
  rotation). The velocity-to-coordinate map qdot = N(q) w is configuration dependent and satisfies
  q . qdot = 0 (qdot is **tangent** to the sphere); the generalized speed w (the body angular velocity,
  3 numbers) lives in the tangent space / Lie algebra so(3), while the coordinate q (4 numbers) lives
  on the manifold. A straight-line step in the four ambient components leaves the sphere (|q| != 1) and
  must be repaired in one of exactly two ways:
  - *project back* (Simbody): take the flat Taylor step in R^4, then renormalize q <- q/|q|. The
    projection is not energy-neutral and depends on the quaternion second derivative qddot.
  - *move along the manifold* (the path used here): rotate q0 by the exact rotation
    expmap(wHalf, dt), the Lie-group exponential from so(3) to SO(3) lifted to the quaternion double
    cover, with the midpoint angular velocity wHalf = u0 + (dt/2) udot0 (giving a second-order,
    time-symmetric step). This never leaves S^3 -- the subsequent renormalization corrects only ~1e-16
    of rounding -- and it advances the orientation **purely from the angular velocity, bypassing
    qddot**.

  The exp-map is a **deliberate divergence from Simbody** (Section 3.4, Section 5.2). In a fully
  faithful engine the two repairs agree to O(dt^3) and either is fine; in this engine the project-back
  path leans on qddot, and the port's Free-joint quaternion kinematics (N, Ndot, qddot) are not yet
  byte-faithful to Simbody, so that path injects kinetic energy through the root body and -- via the
  articulated recursion V_GB[b] = Phi V_GB[parent] + H u -- through the entire tree (an observed
  one-sided KE pump of order 1e3 kJ/mol per trajectory even at dt = 1 fs). The exp-map is robust to
  that latent inconsistency because it never evaluates qddot. The inconsistency itself remains an open
  item to be closed by golden-testing the Free-joint kinematics against Simbody on a single free body
  (Section 3.4); until then the exp-map is the propagator of record. It changes only the proposal, not
  the target (Section 5.4).

  Note that what makes the orientation special is its **dimension and curvature**, not merely that it
  is a "rotation": the Pin torsion is also rotational but lives on the flat S^1, so it takes the flat
  path. Only the genuinely curved S^3 needs the manifold step.

**Corrector convergence as a free "dt-too-large" detector.** The velocity corrector solves the
implicit trapezoid u1 = u0 + (dt/2)(udot0 + udot(q1,u1)) by functional iteration. That iteration is a
contraction only when, roughly, (dt/2) * ||d udot / d u|| < 1, where d udot / d u collects the
velocity-dependent (Coriolis/gyroscopic) coupling and the inverse metric M(q)^-1. There is therefore a
**convergence radius in dt**: below it the iteration contracts in two or three sweeps; above it the
iterates stop contracting (the "change increased after iteration 1" branch) or fail to reach the
tolerance within the iteration cap. Detecting this costs nothing -- the integrator already computes the
relative change each sweep -- and it is a sharp, *local* indicator that dt exceeds what the current
configuration can integrate as a converged, reversible map.

The implementation uses this as a hard guard, deliberately **stricter than Simbody**. Two failure modes
are distinguished and handled differently:

- *Steric clash / runaway (non-finite force, or |u| amplifying past a cap).* A hard overlap gives
  OpenMM an infinite (r^-12) force; this is a **normal, transient** event in sampling. The step is
  **rejected as a move** (the proposal is discarded and the chain keeps its previous point), exactly as
  the Metropolis test would do with the resulting infinite dH. It does **not** abort the run.
- *Genuine non-convergence (finite, bounded iterates that will not reach the tolerance).* This is the
  pure "dt too large for this configuration" signal: the proposal map is no longer the converged,
  reversible trapezoidal map that Sections 5.5-5.6 require, so silently taking it would pump energy.
  The integrator **throws and terminates** with a diagnostic message reporting dt, the relative change
  reached, the tolerance, and the iteration count. At a correctly chosen dt this never fires (the
  corrector converges every step), so the throw is a tripwire for a misconfigured timestep, not a
  runtime branch the sampler relies on.

**Reversibility is configuration dependent; a startup check is a smoke test, not a guarantee.** It is
tempting to validate dt once, at the start, with a round-trip reversibility probe (integrate n steps
forward, flip the momenta, integrate n steps back, flip again, and measure how far the state returned:
~machine epsilon for a reversible map, O(trajectory size) when dt is too large; quaternion DOF compared
with the double-cover-aware distance min(|q-q0|, |q+q0|) because q and -q coincide on S^3). Such a probe
is implemented and is useful, but it certifies dt **only for the configuration it is run from**, for a
reason that must not be misread:

- It is **not** the curvature of S^3 that varies -- that curvature is *constant*, so the orientation
  manifold's intrinsic geometry is identical everywhere.
- What varies is the **mass metric M(q)** (the articulated inertia depends on every torsion) and the
  **local force stiffness d^2 V / d q^2** (a near-clash or compressed region is far stiffer than an open
  one). The largest reversible dt scales like 2 / omega_max with
  omega_max = sqrt( lambda_max( M(q)^-1 d^2V/dq^2 ) ), and **both factors are functions of q**.

Consequently a dt that is perfectly reversible in an open conformation can become non-reversible -- and
begin pumping energy -- when the chain visits a stiffer geometry. The startup reversibility probe answers
only "is dt in the right ballpark for the initial structure?"; it cannot certify the trajectory. The
genuine, ongoing guard is the **per-step corrector-convergence throw above**, which tests the *current*
configuration on *every* step. The reversibility probe is therefore best used as (i) a one-time startup
sanity check and (ii) an optional periodic probe from the current geometry, with the per-step throw as
the always-on safety net. Its cadence is exposed to the user (Python:
`world.set_reversibility_check(interval)`, backed by `SamplerConfig.reversibilityCheckInterval`) and
defaults to **off** (interval 0, zero overhead); a positive interval N runs the probe every N rounds,
starting at round 0 so the first invocation doubles as the startup check, integrating mdSteps
forward+back at the world's timestep and logging the relative round-trip residual. This is also why
the multiscale-dt freedom of Section 5.6 is bounded
per-configuration rather than by a single global number: the usable dt is a property of where in
configuration space the block currently is.

### 5.8 Nonequilibrium candidate moves for dense/caged systems (NCMC)

Sections 5.1-5.7 describe the single HMC move. This section adds a **distinct move
type** for the regime where that move fails not because of stiffness or the time
step, but because the proposed configuration is **caged** by other molecules.

**The problem (contact density, not solvent).** A torsional block mobilizes one
molecule's torsions while every other molecule is welded (Section 3.2). A single
isolated molecule therefore integrates on a smooth intramolecular surface and
tolerates a large dt (Section 5.6). The moment a *second* dense body is in contact
-- explicit solvent, a lipid bilayer/nanodisc, a binding partner, or another
subunit of a large assembly -- every torsional displacement drives the mobile
atoms into the **frozen** atoms of that body. The potential rises on an r^-12 wall:
the corrector convergence radius collapses (dt is forced down, Section 5.7) and any
step that does integrate produces a large dH and is rejected. The failure scales
with the number of contacting atoms and is identical for water, lipid, or a
neighbouring protein -- they are all dense frozen cages. **This is an efficiency
failure, not a correctness one** (Section 5.6): the MH test still targets the exact
canonical distribution; the proposals simply stop moving.

**The move (per-molecule alchemical decouple-move-recouple).** Following Nilmeier,
Crooks, Minh & Chodera (2011), the cage is removed *during the proposal* by a
nonequilibrium switch. The acceptance Hamiltonian's potential is made a function of
a coupling parameter lambda that scales the molecule's **intermolecular** nonbonded
with every other molecule (its intramolecular force field is untouched, so its own
torsional barriers and 1-4 contacts are preserved). One move is a protocol of
interleaved **perturbation** steps (change lambda at fixed configuration) and
**propagation** steps (one fixed-step internal-coordinate Verlet step at fixed
lambda, Section 5.2):

```
lambda: 1 -> 0  (decouple)   ... stride at lambda~0 (uncaged) ...   0 -> 1 (recouple)
```

Near lambda = 0 the molecule's torsions stride freely through the region occupied
by the frozen cage; as lambda returns to 1 the mobile DOF relax to a configuration
compatible with the environment. The frozen environment itself does not move within
this move -- it relaxes in the fully-flexible Cartesian world of the same Gibbs scan
(DOF coverage, Section 13.5, condition C4). Letting a *local environment shell*
co-move (a mixed Pin + Cartesian world) would let the cage yield inside the move
too; that is an efficiency refinement, not a correctness requirement.

**Acceptance (the load-bearing detail).** The propagation is the same fixed-step
Verlet of Section 5.2: deterministic, time-reversible, and volume-preserving on the
constraint manifold (Section 5.5). For such an integrator the conditional path
action vanishes (Nilmeier et al. 2011, Eq. 20; their bistable-dimer Eq. 28), so the
move is accepted on the **full Hamiltonian difference between the two lambda = 1
endpoints**, *not* on the protocol work:

```
accept with probability min(1, exp(-beta * (H_end - H_start))),
H = V + K + F + J   (Section 6),   both endpoints evaluated at lambda = 1.
```

Three consequences must be read precisely, because the natural-looking alternative
(adding a separate work term) is wrong:

- **Work is implicit, never added.** For a deterministic, reversible,
  volume-preserving propagator the protocol work w (sum over perturbation steps of
  the lambda-induced dV) equals H_end - H_start; the two are the *same* quantity, so
  using `H_end - H_start + w` double-counts. The NCMC benefit is already inside
  H_end - H_start: a slow protocol lets the mobile DOF relax as lambda returns, so
  the molecule lands in a low-energy, cage-compatible configuration and H_end is
  small; a fast protocol does not, and H_end is large. The accumulated w is retained
  only as a **diagnostic**, and its discrepancy from H_end - H_start is a free check
  on corrector convergence (it should be small).
- **F and J enter only at the endpoints.** The Fixman determinant F = (1/2) R T ln
  det M(q) and the external-rotation Jacobian J depend on geometry, not on lambda
  (lambda scales nonbonded V only; the articulated inertia M is unchanged). They are
  therefore absent from the protocol work and appear in H_start and H_end exactly as
  in the ordinary torsional move (Sections 6, 10). Likewise the alchemical
  (thermodynamic) perturbation changes no coordinate, so its phase-space Jacobian is
  unity (Nilmeier et al. 2011, "Thermodynamic Perturbation"), and the symmetric
  1 -> 0 -> 1 protocol is its own reverse, so the protocol-selection ratio is unity.
- **It reduces to the ordinary move.** With lambda held at 1 throughout, every
  perturbation dV = 0, the protocol is a plain Verlet trajectory, and the test
  collapses to the Section 5.3 acceptance min(1, exp(-beta dH)). This is the exact
  regression check for an implementation.

**Momentum handling.** Because momenta are fully resampled at the start of every
block (Section 5.1) -- itself a Gibbs move on the velocity marginal -- the NCMC
momentum-reversal that a momentum-persistent chain would require is unnecessary
here; Nilmeier et al. (2011) note exactly this option ("reinitialize velocities from
the Maxwell-Boltzmann distribution after each NCMC step"). The move is thus an HMC
move with a lambda-protocol trajectory, and inherits the Section 13.0
detailed-balance argument unchanged, with the proposal map being the protocol.

**Relation to the rest of the framework.** NCMC is the guidance/acceptance split of
Section 5.4 taken to its limit: the *proposal* is generated under a softened
(decoupled) potential, while the exact atomistic H is retained for acceptance, so
the sampled distribution is unbiased regardless of how aggressively the cage is
softened (Section 5.4). It composes with the mixed Gibbs scan exactly like any other
pi-invariant block (Section 13.4); a fully flexible Cartesian world must still appear
in the cycle to relax the environment and cover all DOF (Section 13.5). The current
implementation softens a molecule's intermolecular nonbonded for the non-periodic
(vacuum/implicit) case -- the large-assembly regime where caging appears even
without solvent; periodic/PME alchemy (reciprocal space couples all atoms) is a
separate extension. Chen & Roux (2014, 2015) give the closely related hybrid
nonequilibrium-MD/MC formulation and its symmetric-momentum-reversal acceptance.

---

## 6. The acceptance Hamiltonian function

```
H(q, p) = V(q) + (1/2) * p^T * M(q)^-1 * p + F(q) + J(q) + W_BAT(q)
```

The mass-metric M(q) enters H in **two** places: the **kinetic term** (1/2) p^T M^-1 p -- the
Riemannian cometric (inverse-metric) form of the kinetic energy on the restricted space -- and the
**Fixman potential** F, through ln det M (the log-volume of that metric). This is the only way the
configuration dependence of the metric reaches the acceptance. The last term W_BAT is the
internal-coordinate volume-element Jacobian; it is **identically zero (a constant that cancels) for
the implemented torsional and Cartesian worlds**, so for those worlds H = V + K + F + J as written
throughout this document. It is nonzero only for the planned Ball/Cyl worlds, which mobilize a hard
internal coordinate (see the W_BAT bullet below).

- **V(q)** -- potential energy (force field + chosen solvent model; Section 8), full atomistic model.
- **K = (1/2) p^T M(q)^-1 p** -- kinetic energy; M(q) is the reduced articulated mass-metric of the
  block (Section 3).
- **F(q) -- Fixman compensating potential.**

  ```
  F(q) = (1/2) * R * T * [ ln det( M(q) )  -  ln det( G(q) * M(q)^-1 * G(q)^T ) ]
  ```

  The first term cancels the det(M(q))^(1/2) factor that emerges when the Gaussian momenta are
  marginalized, so that the configurational marginal matches the target without the mass-metric
  artifact. The second term is present **only when holonomic ring-closure constraints are active**:
  the RATTLE momentum projection onto G*M^-1*p = 0 injects a det(G*M^-1*G^T)^(-1/2) factor into the
  marginal, which this term cancels. For acyclic molecules the second term vanishes, and F reduces
  to the tree determinant (ln det M = sum_b ln det(D_b)).

  *The reference, made exact.* F as written is a **difference of log-determinants relative to a
  reference**; the physically exact reference (Jain et al. 2013, eq. 11) is the **BAT** mass matrix
  M_B = J_BAT^T M_3N J_BAT, i.e. F_exact = (1/2)*R*T*[ ln det M - ln det M_B - ln det(G*M^-1*G^T) ].
  By Go-Scheraga (Jain et al. 2013, eq. 6-7), det M_B = \|J_BAT\|^2* det M_3N with
  \|J_BAT\| = (product_i r_i^2)(product_j sin theta_j)*const, so the reference splits into a truly
  constant Cartesian part ln det M_3N (cancels in dH unconditionally) and a configuration-dependent
  BAT-volume part -ln\|J_BAT\|^2. **Under torsional dynamics every r_i and theta_j is frozen, so
  -ln\|J_BAT\|^2 is also constant and drops** -- which is why F reduces to (1/2)*R*T*ln det M with
  the Cartesian matrix M_3N as the effective reference, and why ln det M_3N "cancels in dH." This
  shortcut is **valid only when no hard internal coordinate is mobile.** The moment a world mobilizes
  a bond or angle (the planned Ball/Cyl worlds), the corresponding factor of \|J_BAT\| is no longer
  constant and must be restored explicitly -- this is exactly the W_BAT term below.
- **W_BAT(q) -- internal BAT-volume Jacobian (Ball/Cyl worlds only).**

  ```
  W_BAT(q) = -R*T * [ sum_{mobile angles k} ln sin theta_k  +  sum_{mobile bonds i} ln r_i^2 ]
  ```

  This restores the configuration-dependent part of the BAT volume element \|J_BAT\| for any world
  that mobilizes a *hard* internal coordinate. The correct conditional of the flexible target over a
  world's mobile coordinates s is pi_flex(s | h*) proportional to exp(-beta*V) * \|J_BAT(s)\| ds
  (Section 13.3); after the Fixman term cancels det(M)^(1/2) a world targets exp(-beta*V) ds, so the
  residual factor \|J_BAT(s)\| must be supplied by exp(-beta*W_BAT) = \|J_BAT(s)\|, giving the form
  above. It is the **internal-coordinate analogue of J(q)** (which supplies the *external* rotational
  Haar factor): a Ball world mobilizing a bond angle theta_k carries -R*T*ln sin theta_k (equivalently
  -(1/2)*R*T*ln sin^2 theta_k, the same downward-divergent form as J, evaluated with the same eps
  floor near sin -> 0); a Cyl world mobilizing a bond r_i carries -R*T*ln r_i^2 = -2*R*T*ln r_i. For
  the **torsional** world all bonds and angles are frozen, so every term is constant and W_BAT
  contributes nothing to dH; for the **Cartesian** world there is no BAT chart at all (sampling is in
  the flat Cartesian measure dx directly, Section 13.1), so W_BAT = 0 there too. W_BAT is therefore
  identically inactive in the currently implemented sampler and becomes live only when Ball/Cyl are
  ported. Omitting it in a Ball/Cyl world makes that world target the wrong conditional and breaks the
  composition argument of Section 13.4 (condition C1).
- **J(q) -- external-rotation Jacobian.**

  ```
  J(q) = -(1/2) * R * T * sum_{free root bodies b} ln sin^2( gamma2_b )
  ```

  This is the volume-element Jacobian that restores the uniform (Haar) orientation measure for any
  robot whose root body is **free** (6 external DOF). gamma2_b is the polar ("pitch") angle of body
  b's orientation, read from its orientation quaternion (quaternion -> rotation -> pitch). To stay
  finite at the orientation singularity (pitch -> 0 or pi, where sin^2 -> 0), the implementation
  evaluates ln of max(sin^2(gamma2_b), eps) with a small floor eps (1e-12); this caps J near the
  pole without affecting the bulk measure. Translational external DOF contribute a unit Jacobian.
  J(q) is **identically zero when all roots are welded** (Section 3.3); it is nonzero only for
  free/floating bodies (e.g. a docked ligand), and is the validated form used by the program. The
  bond-length (r^2) and bond-angle (sin theta) factors of the BAT volume element are **constant**
  under torsional dynamics (those coordinates are frozen), so they cancel in dH and do not appear in
  H; if a world mobilizes them, they re-enter through W_BAT above.

**Cartesian special case.** In Cartesian coordinates M = M_3N is the constant diagonal atomic-mass
matrix, so F and J are configuration-independent constants (taken as zero) and the move reduces to
standard Cartesian HMC/MD.

The MH acceptance uses min(1, exp(-beta * dH)) with the **exact** dH between the pre- and
post-trajectory states, including F and J (and W_BAT in a Ball/Cyl world).

---

## 7. Constraints

Holonomic constraints are applied **only in internal coordinates** and **only on ring-closing
(cotree) bonds**, which cannot be represented as tree coordinates. A holonomic constraint depends on
positions q **only** -- sigma(q) = \|r_AB\|^2 - d0^2 = 0 is a single algebraic relation in q -- and
each independent one removes one dimension from the configuration manifold (Section 5.5). Each is
enforced at both levels of a RATTLE scheme:

- **Position (SHAKE):** after the Verlet position update, q is projected so that sigma(q) = 0.
- **Velocity (RATTLE):** after the velocity update, the momenta are projected so that
  G(q)*M(q)^-1*p = 0 (no relative velocity along the closed bond).

Both projections are required; either alone leaves a secular drift. RATTLE preserves the phase-space
measure **on the constraint manifold**, so the volume-preservation premise of the HMC move
(Section 5.5) holds with constraints active. The constraint contributes the det(G*M^-1*G^T) term to
F (Section 6).

The ring-closure distance d0 is part of the configuration carried between blocks: it tracks the
value left by the previous block (continuation), and is **not** reset to an idealized force-field
bond length.

### 7.1 The projection algorithm (exact specification)

Grounded in Simbody's loop-closure solver (RBDistanceConstraint + LengthConstraints, the
Newton-Raphson loop solver) and in SHAKE/RATTLE (Ryckaert et al. 1977; Andersen 1983). All three
users of the constraint Jacobian -- SHAKE, RATTLE, and the loop-closure Fixman determinant
(Section 6) -- assemble it identically, which is what guarantees they stay mutually consistent.

**Per-constraint Jacobian assembly (shared).** For each loop-closure bond c = (A, B) with current
r_AB = posG[A] - posG[B]:

1. Place the Cartesian force +r_AB on atom A and -r_AB on atom B. (This is the (1/2) dC/dr
   convention: the factor 2 in dC/dr = 2 r_AB is folded out, uniformly, everywhere.)
2. Reduce these atom forces to a generalized force by one inward rigid force-transmission sweep --
   the transpose of forward kinematics: accumulate each body's spatial force (linear = sum of atom
   forces; angular = sum of (r_atom - r_Bo) x f), sweep child -> parent shifting by the Phi^T
   operator, and read the per-DOF generalized force G_c^T = H_b . spatialForce_b. This is one row of
   G^T (an nu-vector).
3. Form M^-1 G_c^T with one O(n) articulated multiplyByMInv call.

Assemble the m x m coupling matrix W = G M^-1 G^T, with W[i][j] = G_i . (M^-1 G_j^T) and m the
number of loop closures (m is small; W is solved directly by symmetric elimination, which also
yields ln|det W|).

**RATTLE (velocity projection), after the velocity update.** Solve W * lambda = G u, where
(G u)_c = r_AB . (v_A - v_B) is the relative velocity along bond c (v_A, v_B read from the spatial
body velocities V_GB). Update

```
u  <-  u  -  sum_c  lambda_c * (M^-1 G_c^T).
```

After this, G u = 0: no relative velocity along any closed bond.

**SHAKE (position projection), Newton iteration after the position update.** Iterate (to a maximum
count) until max_c |C_c| < tol:

```
refresh forward kinematics (q -> X_GB -> posG)        (Section 3.4)
C_c   = |r_AB|^2 - d0_c^2                              (holonomic residual)
assemble G^T, M^-1 G^T, W                              (shared assembly, at current q)
solve  W * lambda = C
du    = - sum_c lambda_c * (M^-1 G_c^T)                (a velocity-space increment)
apply du to q:  dq = N(q) du for quaternion DOF (Section 3.4), dq = du otherwise; renormalize quats
```

This realizes the generalized-coordinate SHAKE step dq = -M^-1 G^T (G M^-1 G^T)^-1 C, iterated to
convergence. Both projections keep (q, p) on the constraint manifold -- the cotangent bundle of
sigma(q) = 0 -- and RATTLE preserves the measure on it (Section 5.5).

**Loop-closure Fixman term.** ln det(G M^-1 G^T) is read from the same W factorization and enters F
as -(1/2) RT ln det(G M^-1 G^T) (Section 6). Because every user assembles G in the identical
(1/2 dC/dq) convention, adopting the full dC/dq convention instead would scale W by 4^m and shift
the log-determinant by the configuration-independent constant m*ln 4, which cancels in dH; the
convention is therefore immaterial to sampling, but must be the **same** in all three users.

**Complexity.** Each multiplyByMInv and force-transmission sweep is O(n); there are m of them per
SHAKE iteration, plus the O(m^3) dense solve of W. For the small m of macrocyclic ring closures the
per-iteration cost is dominated by O(m*n).

---

## 8. Forces, energy, and solvent

- Energies and forces are evaluated on the **full atomistic model** (not coarse-grained) through
  **OpenMM**, using an **AMBER or CHARMM** force field. Per-atom Cartesian forces -grad V are summed
  into per-body spatial forces for the articulated-body solver (Section 3).
- The **solvent model is configurable and not mandated by the theory**: vacuum, implicit, or
  explicit. The current default implicit model is **GBSA-OBC (Onufriev-Bashford-Case, model II /
  OBC2)**, OpenMM's standard implementation. Choosing a different solvent (or an explicit box with
  its accompanying ensemble terms) changes V (and possibly the target weight, Section 2) but leaves
  the Gibbs/HMC construction untouched.
- Because V is evaluated on the Cartesian configuration, it is **independent of which generalized-
  coordinate chart a block uses**; the internal-coordinate machinery is a parameterization of the
  proposal, not of the energy (relevant to the inter-block handoff, Section 9).

---

## 9. Inter-block continuation

State is handed from one block to the next as the full configuration in **Cartesian coordinates**
(which encode every bond length, angle, and torsion). Cartesian coordinates are the shared
interchange representation precisely because the energy is chart-independent (Section 8).

On entering a block, that geometry becomes the block's reference. The block's generalized
coordinates Q are **set from the incoming Cartesian configuration** -- the operation analogous to
"set the default Q to fit a given configuration" (Simbody's `setQToFitTransform`/default-state idiom
in the earlier implementation) -- and the block's rigid bodies and mobilizer frames are rebuilt from
the **actual incoming Cartesian geometry**. Any coordinate the block freezes is thus held at its
**current carried-over value**, not an idealized one. The handed-over state is the previous block's
post-MH configuration -- the new point on acceptance, the retained point on rejection. Momenta are
not carried over; they are freshly resampled in each block (Section 5.1).

**Does this re-charting introduce bias?** No, provided the Cartesian -> Q assignment is exact for
the block's rigid-body definition, which it is here. Because the frames and rigid-body shapes are
rebuilt from the *actual* incoming coordinates (rather than reset to idealized bond lengths/angles),
the incoming configuration lies exactly on the new block's constraint manifold, so the round-trip
Cartesian -> Q -> Cartesian is the identity: no projection, no information loss, no change of
measure. The handoff is a **change of chart**, and since V is chart-independent and each block
applies the Fixman/J terms appropriate to its own chart, the target distribution is unchanged.

Bias *would* arise only if the handoff projected onto a **different** manifold -- e.g. if entering a
torsional block reset bonds/angles to idealized force-field values instead of reading the actual
incoming geometry. That would silently move the system onto an idealized slice and distort the
sampled distribution. Robosample avoids this by rebuilding from actual coordinates; this is the
implementation-level guarantee behind condition (C3) of Section 13.6.

---

## 10. Worked example: dH for one torsional block (acyclic, welded root)

This is the most common case -- the torsional special case of the general HMC move (Sections 5-6) --
with a welded root (so J = 0) and no ring closures (so the G*M^-1*G^T term of F is absent).

1. The block receives configuration q_old (Section 9) and evaluates V_old = V(q_old).
2. Momenta are resampled: p_old ~ N(0, R*T*M(q_old)). Then

   ```
   K_old = (1/2) * p_old^T * M(q_old)^-1 * p_old
   F_old = (1/2) * R * T * ln det M(q_old)
   H_old = V_old + K_old + F_old
   ```

3. The fixed-step Verlet integrator (Section 5.2) integrates the block's mobile torsions under V
   only (no Fixman torque), producing (q_new, p_new).
4. Evaluate

   ```
   V_new = V(q_new)
   K_new = (1/2) * p_new^T * M(q_new)^-1 * p_new
   F_new = (1/2) * R * T * ln det M(q_new)
   H_new = V_new + K_new + F_new
   ```

5. Accept with probability min(1, exp(-beta * dH)), where

   ```
   dH = (V_new - V_old) + (K_new - K_old) + (F_new - F_old)
   ```

Notes:

- det M(q) enters dH only through F, never through the momentum draw (Section 5.1).
- The bond-length, bond-angle, and (trivial) torsion Jacobian factors are constant and cancel, so
  they do not appear in dH.
- With a ring closure present, add the constraint contribution
  -(1/2)*R*T*[ ln det(G*M^-1*G^T)|_new - ln det(G*M^-1*G^T)|_old ] to dH.
- With a free root body, add J_new - J_old (Section 6) to dH.
- This example is torsional, so the BAT bond/angle factors are constant and W_BAT = 0. In a Ball/Cyl
  world (a hard internal coordinate mobile) add W_BAT_new - W_BAT_old (Section 6) to dH as well.

---

## 11. Mobilizer (joint) types and degree-of-freedom coverage

The world/block machinery is defined in Sections 3-4; this section enumerates the **joint
(mobilizer) types** a block may use, states which are wired into the current engine, and gives the
degree-of-freedom (DOF) coverage requirement. Each joint mobilizes a specific restricted space and
places its frames along bonds by default (placing them away from atom centers degrades transition
rates).

**Implemented in the current engine.** Four mobilizers are wired end to end -- the frame chain of
Section 3.4, the velocity hinge map H, the qdot map, and the model builder:

- **Weld** (0 DOF). Rigidly fixes the body to its parent; X_FM = identity. Used for rigid-body
  lumping (Section 3.2) and for a fixed root. Contributes det(D_b) = 1, i.e. nothing to the Fixman
  sum.
- **Pin / torsional** (1 DOF). Rotation about the mobilizer Z axis -- the bond, after the bond->Z
  alignment of Section 3.4 -- with bonds and angles welded. The torsional case; qdot = u. Large,
  low-frequency moves; M is the reduced articulated metric and Fixman is active.
- **Translation / Cartesian** (3 DOF). Three orthogonal slides; over the mobilized atoms M = M_3N
  (constant) so F = J = 0. The fully-flexible world is built from per-atom Translation joints;
  qdot = u.
- **Free** (6 DOF). A floating root: quaternion orientation plus translation. The orientation
  contributes J(q) (Section 6) and uses the parent-frame quaternion map of Section 3.4; qdot != u
  for the four quaternion components.

So the Cartesian (fully flexible) world is the Translation set and the torsional world is the Pin
set.

**Declared but not yet ported.** The joint enum also lists Slider, Cylinder, BendStretch, Ball
(spherical), SphericalCoords, and FreeLine, but the engine does not yet provide their X_FM / H_FM /
qdot kernels and does not build them, so they currently fall through to identity and must not be
selected. Two are the design target carried over from the 2020 Robosample paper and remain the
intended extension:

- **Spherical / "Ball"** (3 rotational DOF, quaternion). Would mobilize a torsion together with its
  adjacent **bond angle** -- **torsion-angle coupling**. Because it mobilizes a hard internal
  coordinate (the angle), its acceptance Hamiltonian **must** carry the internal BAT-volume term
  W_BAT = -R*T * sum_{mobile angles k} ln sin theta_k (Section 6); the Fixman term (1/2)*R*T*ln det M
  alone targets the wrong conditional for this world.
- **Cylindrical / "Cyl"** (2 DOF: one translation + one rotation about the translation axis). Would
  mobilize a torsion together with its adjacent **bond length** -- **torsion-bond coupling**.
  Likewise, mobilizing a hard bond coordinate requires the internal BAT-volume term
  W_BAT = -R*T * sum_{mobile bonds i} ln r_i^2 = -2*R*T* sum_i ln r_i (Section 6) in its acceptance
  Hamiltonian.

When ported, Ball (and the rotational part of any spherical joint) shares the Free joint's
quaternion kinematics, so it inherits the orientation-map convention of Section 3.4 and must be
golden-tested against Simbody. Conceptually the Ball and Cyl joints are Robosample's rigorous answer
to the bond/angle-torsion coupling that drives rigid-geometry distortion (Section 12): rather than
re-parameterizing the force field, they mobilize the coupled hard DOF in additional Gibbs worlds, so
the coupling is sampled exactly under the atomistic potential. The W_BAT terms above are part of
"sampled exactly": without them, mobilizing the hard coordinate would correct the piece-2 PMF
distortion (Section 12) only to introduce a fresh measure (piece-1-type) bias, and the world would no
longer leave pi_flex invariant (Section 13.4, condition C1).

**DOF coverage (ergodicity).** The chain targets the full flexible Boltzmann distribution iff the
**union** of mobile DOF across the worlds in a cycle covers every DOF (Section 13.5, condition C4).
No single constrained world covers everything -- a torsional world freezes bonds and angles, a Ball
world (when available) freezes bond lengths, a Cyl world freezes bond angles -- so the simplest
guarantee is one fully flexible Cartesian (Translation) world per cycle.

---

## 12. The rigid-geometry (metric) bias

Any method that samples in a reduced (e.g. torsional) space while holding bonds/angles fixed must
confront one fact: **the equilibrium distribution of the soft coordinates in a model with rigid
bonds/angles is not the same as their marginal in the fully flexible model** (Fixman 1974; van
Gunsteren & Berendsen 1977; Patriciu, Chirikjian & Pappu 2004; Echenique, Calvo & Alonso 2006). It
has two distinct pieces, frequently conflated:

1. **A measure (kinetic / metric) piece.** Marginalizing the Gaussian momenta of a constrained
   system leaves a configuration-dependent factor det(M(q))^(1/2). Uncorrected, the sampler targets
   exp(-beta*V)*det(M)^(1/2) instead of exp(-beta*V). This is what the **Fixman potential** removes
   (Section 6); it is analytic and force-field-agnostic. Patriciu et al. (2004) showed its
   importance grows with chain length; Echenique et al. (2006) found it matters beyond two residues.
   (Patriciu et al.'s *headline* recommendation -- ignore det(M) in the acceptance ratio -- applies
   to **purely positional** torsional Monte Carlo, where momenta are never sampled; they explicitly
   confine the determinant bias to torsional-space molecular dynamics in which the soft-mode momenta
   are computed explicitly and the frozen-coordinate momenta are set to zero. Robosample's HMC
   resamples momenta every block (Section 5.1), so it is squarely in that latter regime and **does**
   require F. The citation supports the chain-length scaling, not an argument against the correction.)

2. **A potential-of-mean-force (PMF) piece.** In a flexible molecule, bond angles *relax* in
   response to the torsional configuration (an angle widens to relieve a 1-4 clash at a given
   torsion). Freezing the angle removes that relaxation. This is an energetic effect -- the location
   of the bond/angle minimum moves with torsion -- and **no determinant correction captures it**. It
   is recovered only by letting the coupled hard coordinates move.

**How Robosample handles each piece.** Piece 1 is removed exactly by the Fixman potential in the
acceptance Hamiltonian. Piece 2 is recovered by **mobilizing the coupled hard DOF** -- through the
**fully flexible Cartesian world**, which alone covers every DOF and so suffices for correctness
(condition C4, Section 13.6), and, when ported, the Ball (angle) and Cyl (bond) worlds (Section 11)
as an efficiency refinement -- inside the Gibbs scan. This is the rigorous analogue of Kandel et
al.'s (2016) selective angle relaxation: where they relaxed angle constraints inside a single
trajectory (exact only if the relaxed DOF are independent of the rest), Robosample relaxes them in
**separate Gibbs worlds** corrected by an exact MH test, so no independence assumption is needed
(Section 13).

### 12.1 Why a torsional Gibbs block needs no internal-coordinate force-field correction

Chen, Im & Brooks (2005) showed that running torsion-angle dynamics with rigid covalent geometry
**grossly distorts** the torsional energy surface of Cartesian-parameterized force fields, and they
(like Katritch, Totrov & Abagyan 2003) repaired it by building a specialized **internal-coordinate
force field (ICFF)** -- modified torsion (CMAP-style) terms plus softened van der Waals and
electrostatics -- that approximates the source Cartesian field *without* a compensating potential.
Robosample needs no such correction, for two independent reasons:

- **The proposal is corrected by an exact acceptance test.** Chen et al.'s TAMD has no
  accept/reject step: its trajectory *is* the sample, so any distortion in the rigid-geometry
  forces enters the result, and the force field must be repaired to compensate. In Robosample the
  rigid-geometry dynamics is only a **proposal**; the MH test uses the full atomistic V (plus
  Fixman), so the sampled distribution is exact regardless of how distorted the proposal forces are
  (the distinct guidance/acceptance design, Section 5.4). At worst a poor proposal lowers
  acceptance; it cannot bias the result.
- **The frozen-geometry PMF distortion is integrated out.** The residual piece-2 distortion from
  holding bonds/angles fixed is removed not by patching V but by **mobilizing those coordinates** in
  the Cartesian/Ball/Cyl worlds of the scan (above). Pure TAMD cannot do this -- it never relaxes
  the frozen coordinates -- which is exactly why Chen et al. had to bake the relaxation into the
  force field.

Consequently Chen's ICFF is, for Robosample, **optional and orthogonal**: it could be used purely as
a *guidance* potential to ease 1-4-contact barriers at low temperature (with the atomistic V
retained for acceptance, Section 5.4) -- an anticipated efficiency extension, never a correctness
requirement. Note also that the native ff19SB / CHARMM backbone CMAP is simply part of V and is
applied unchanged in every world; it is a **different object** from the ICFF CMAP and does not
address rigidity.

---

## 13. Correctness of the mixed Gibbs scan (assumption-free)

This section proves that the scan leaves the fully flexible Boltzmann distribution invariant and
states every assumption. It is the rigorous form of the stationarity claim of Section 4 and the
marginal-matching claim of Section 2, and it defines intended behavior: an implementation that
violates a listed condition is outside the regime the theory guarantees.

### 13.0 Invariance, reversibility, detailed balance

Three properties of a Markov kernel K(x -> y), in increasing strength:

- **Invariance (stationarity).** pi is invariant for K if pi*K = pi, i.e. integral pi(x) K(x -> y)
  dx = pi(y). pi is a fixed point: once the chain is distributed as pi, it stays so. This is the
  **only property required** for the scan to sample pi.
- **Detailed balance.** K satisfies detailed balance with respect to pi if
  pi(x) K(x -> y) = pi(y) K(y -> x) for all x, y. This is a **sufficient but not necessary**
  condition for invariance (integrate both sides over x to recover pi*K = pi).
- **Reversibility.** A chain is reversible w.r.t. pi exactly when it satisfies detailed balance:
  forward and time-reversed transitions are equiprobable under pi. Hence **reversibility implies
  invariance, but invariance does not imply reversibility** -- a systematic-scan Gibbs sweep is
  invariant yet not reversible (its reverse is the reversed sweep).

**"Reversible by construction" (Section 4).** Each block is an HMC move built so that it satisfies
detailed balance w.r.t. pi: momentum resampling is an exact Gibbs draw on p (Section 5.1); the
proposal map is reversible and volume-preserving (Section 5.5); and the MH test uses the exact
acceptance Hamiltonian (Sections 5.3, 6). These three are precisely the ingredients of the standard
HMC detailed-balance proof, so the block leaves pi invariant by construction. The composite scan
inherits invariance (Section 13.4) but generally not reversibility.

### 13.1 Setup and notation

- Full configuration q as Cartesian coordinates x (for a welded-root single molecule, modulo overall
  translation/rotation; for a free root, including its 6 external DOF e).
- Target: pi_flex(x) proportional to exp( -beta * V(x) ) (Lebesgue measure dx).
- BAT split: hard coordinates h = (bonds b, angles theta) and soft coordinates tau (torsions), plus
  external e for a free root.
- BAT volume element factorizes and is **independent of the torsions**:

  ```
  dx = |J_BAT(b, theta)| * db * dtheta * dtau * de,
  |J_BAT| = ( product_i r_i^2 ) * ( product_j sin theta_j ) * const
  ```

### 13.2 What one constrained block targets (Fixman lemma)

A block B fixes its frozen coordinates at their incoming values, resamples p ~ N(0, R*T*M(.)),
integrates internal-coordinate fixed-step Verlet (Section 5.2) under V, and applies the MH test on
H_B = V + (1/2) p^T M^-1 p + F [+ J] [+ W_BAT], with F = (1/2)*R*T*ln det M. Internal-coordinate HMC with
H = V + (1/2) p^T M^-1 p leaves invariant the joint proportional to exp(-beta*H); marginalizing the
Gaussian p gives configurational density proportional to exp(-beta*V)*det(M)^(1/2). The Fixman term
contributes exp(-beta*F) = det(M)^(-1/2) (using beta*R*T = 1), which **cancels** the determinant.
Hence block B leaves invariant exp(-beta*V) in the flat measure of its own mobile coordinates (the
Fixman result; constrained marginal on the modified potential = unconstrained marginal). For a
torsional block at fixed h*, the target is exp(-beta*V(tau; h*)) dtau.

### 13.3 Identification with the conditional of pi_flex

Restricting pi_flex to { h = h* } and expressing it in tau:

```
pi_flex(tau | h*)  proportional to  exp( -beta * V(tau; h*) ) * |J_BAT(h*)|  dtau
                   proportional to  exp( -beta * V(tau; h*) )  dtau,
```

since |J_BAT(h*)| is constant in tau (Section 13.1). This is exactly what block B leaves invariant
(13.2), so **a torsional block is a valid Metropolis-within-Gibbs update of the torsions drawing
from the correct conditional pi_flex(tau | h)**.

For the Ball/Cyl worlds (when ported) the argument is **almost** the same but with one essential
difference: those worlds mobilize a hard coordinate (an angle or a bond), and |J_BAT| is **not**
constant in that coordinate. Restricting pi_flex to fixed bonds b* and expressing it in the Ball
world's mobile coordinates (theta, tau),

```
pi_flex(theta, tau | b*)  proportional to  exp(-beta*V) * |J_BAT(b*, theta)| dtheta dtau
                          proportional to  exp(-beta*V) * (product_k sin theta_k) dtheta dtau,
```

the sin theta_k factors over the mobilized angles **survive** (only the frozen r_i^2 factors are
constant and drop). After its Fixman term cancels det(M)^(1/2), a Ball world targets
exp(-beta*V) dtheta dtau, so it reproduces this conditional **only if** the residual factor
product_k sin theta_k is supplied -- which is precisely the role of W_BAT
(= -R*T sum_k ln sin theta_k; Section 6). The Cyl world is the bond-coordinate analogue, with W_BAT =
-2*R*T sum_i ln r_i restoring the product_i r_i^2 factor. With W_BAT included, each such world is a
valid Metropolis-within-Gibbs update of its (angle, torsion) or (bond, torsion) coordinates drawing
from the correct conditional; **without it, the world targets exp(-beta*V) in the flat internal
measure, which is the wrong conditional** (this is the bond-length-coordinate generalization Jain et
al. 2013 flag after their eq. 13, where the simple (1/2) ln det M form no longer suffices). For a
free root, J restores Haar measure on the orientation.

### 13.4 Composition

The Cartesian world C leaves pi_flex(x) invariant. A scan is an ordered composition of worlds, each
of which leaves pi_flex invariant (the constrained ones because leaving every conditional invariant
is equivalent to leaving the joint invariant while only the conditioned-on coordinates are held; C
because it leaves the full joint invariant). A composition of pi_flex-invariant kernels is
pi_flex-invariant. **Therefore the scan leaves pi_flex invariant**, in any (systematic or random)
order (Section 4) -- this is the Gibbs-sampling argument: constrained dynamics samples a conditional,
and alternating constraints samples the joint.

### 13.5 Why hard-DOF coverage is necessary, not optional

A constrained world never moves its frozen coordinates. If some hard coordinate is frozen in
**every** world of a cycle, the chain is confined to a slice of fixed value for it; its stationary
distribution is a conditional, **not** the marginal pi_flex(tau) = integral pi_flex(tau, h) dh. The
gap between them is precisely the rigid-geometry distortion of Chen et al. (2005) and Kandel et al.
(2016). Hence the **union of mobile DOF across the worlds in a cycle must cover every DOF**; the
simplest guarantee is one fully flexible Cartesian world per cycle. A single Ball or Cyl world is
not sufficient on its own (each freezes a different hard coordinate); coverage is a property of the
world **set**, not of any one world.

### 13.6 Conditions (the explicit assumptions)

The guarantee has three tiers, and the conditions map onto them. **C1-C3 give invariance:** each
world -- and hence the composite scan -- leaves pi_flex invariant. **C4 gives uniqueness:** it makes
pi_flex the *unique* stationary distribution; without it the chain is reducible (it still leaves
pi_flex invariant, but it never moves some coordinate, so it stays frozen on that coordinate's
initial value and converges to a conditional, not pi_flex). **C5 gives finite-time correctness:**
the covering world must be invoked often enough and mix well enough that the realized average is
unbiased in practice. All five:

- **(C1) Fixman present and correct in every constrained world.** F = (1/2)*R*T*ln det M, plus the
  -(1/2)*R*T*ln det(G*M^-1*G^T) term when ring closures are active, plus J for free roots, plus
  W_BAT for any world that mobilizes a hard internal coordinate (Ball: -R*T sum ln sin theta_k; Cyl:
  -2*R*T sum ln r_i; Section 6). For the implemented torsional and Cartesian worlds W_BAT is a
  constant and drops, but it is mandatory for Ball/Cyl. Omitting F makes a world target
  exp(-beta*V)*det(M)^(1/2) and biases the chain (a measure bias, distinct from the Chen PMF bias);
  omitting W_BAT in a Ball/Cyl world biases it the same way (a missing internal-volume factor). The
  Fixman *torque* is **not** required (Section 5.4).
- **(C2) One identical V across all worlds.** Same force field, atom set, parameters, and solvent
  model. Mixing e.g. ff19SB in one world and CHARMM in another defines no single pi. (An ICFF
  *guidance* potential is permitted, with the shared atomistic V used for acceptance -- Sections
  5.4, 12.1.)
- **(C3) Reversible, volume-preserving proposal with exact dH, and an exact Cartesian->Q handoff.**
  Each world's map is reversible and volume-preserving on its (constraint) manifold (Section 5.5),
  the MH test uses the exact dH including F and J, and the inter-block re-charting is exact (Section
  9). RATTLE supplies measure preservation on the ring-closure manifold (Section 7).
- **(C4) The union of mobile DOF across the worlds in a cycle covers all DOF** (Section 13.5).
- **(C5) Ergodicity / mixing (finite-time accuracy).** The DOF-covering world(s) must be invoked
  often enough, and mix the hard DOF well enough, that the realized chain explores the
  h-distribution; otherwise the finite-sample marginal is effectively conditioned on an
  under-sampled set of geometries -- a transient bias that mimics the Chen distortion even though it
  vanishes asymptotically. Per-block acceptance diagnostics (Section 4) are the intended detector.

A statement of the form "Robosample samples the fully flexible Boltzmann distribution" is true of
the **composite scan under (C1)-(C5)**. It is **false** of any single-constraint run (one world, or
a world set that does not cover all DOF): that samples a conditional. The Section 2 marginal-matching
claim should be read as a property of the composite scan.

---

## 14. Positioning among related torsional-sampling methods

Torsional sampling has a roughly thirty-year history; methods split mainly by how they treat the two
pieces of the rigid-geometry bias (Section 12). "Exact canonical (Boltzmann) target?" asks whether
the stationary distribution equals the canonical Boltzmann distribution of the *fully flexible*
atomistic model; "No (not intended)" marks tools built for refinement or optimization, where this is
not the design goal and is not a defect.

| Method / program (representative refs) | Coordinates / what is frozen | Primary goal | Treatment of rigid-geometry bias (pieces 1 + 2, Section 12) | Exact canonical target? | Relation to Robosample |
|---|---|---|---|---|---|
| **Fully flexible Cartesian MD / HMC** (reference) | All Cartesian DOF mobile; nothing frozen | Canonical sampling | No bias to correct; pays full stiff-mode cost | Yes (correct thermostat / MH test) | Robosample's Cartesian world is exactly this; one block in the scan |
| **TAMD for crystallographic refinement** -- X-PLOR / CNS (Rice & Brunger 1994; Stein, Rice & Brunger 1997) | Torsions mobile; bonds/angles rigid | Structure refinement vs. data | Ignored -- legitimately: target is a restraint pseudo-energy, not an ensemble | No (not intended) | Different problem (optimization); shares the recursive multibody machinery |
| **TAMD for NMR** -- DYANA / CYANA (Guntert et al. 1997); Xplor-NIH (Schwieters & Clore 2001) | Torsions mobile; bonds/angles rigid | Structure from restraints | Ignored (same rationale) | No (not intended) | Same O(n) equations of motion (Jain et al. 1993); not a Boltzmann sampler |
| **NEIMO / GNEIMO; angle-relaxed ICMD** (Jain, Vaidehi & Rodriguez 1993; Vaidehi, Jain & Goddard 1996; Kandel et al. 2016) | Torsions / relaxed angles mobile; rest rigid | Thermodynamic / kinetic torsional MD | Piece 1 via thermostat + (later) Fixman; piece 2 via selective angle relaxation (exact only if relaxed DOF independent) | Approaches it | Algorithmic ancestor; Robosample relaxes coupled DOF in separate Gibbs worlds (Ball/Cyl) with an exact MH test, dropping the independence assumption |
| **GNEIMO-Fixman** (Wagner et al. 2013; Jain et al. 2013) | Torsions mobile; bonds/angles rigid | Accurate torsional MD with correct Z | Piece 1 by exact Fixman for branched molecules (torque **in the forces**); piece 2 still rigid | Piece-1-exact; piece-2 approximate | Closest physics cousin. Robosample puts Fixman **only in MH acceptance** by default (Duane et al. 1987) and mobilizes coupled DOF for piece 2 |
| **ICFF / TAMD-CMAP** (Katritch, Totrov & Abagyan 2003; Chen, Im & Brooks 2005) | Torsions mobile; bonds/angles rigid | Efficient conformational sampling | Piece 2 **empirically**: modified torsion (CMAP) terms + softened nonbonded, replacing the compensating potential | Approximate; tied to one fixed constraint set | Robosample keeps unmodified V + Fixman (so it can switch constraints); an ICFF may serve as **guidance** with V retained for acceptance (Section 12.1) |
| **Biased-Probability MC / ICM** (Abagyan & Totrov 1994) | Torsions sampled; bonds/angles rigid | Global optimization / prediction | Bias deliberately introduced (phi-psi / rotamer zones); not the canonical ensemble | No (not intended) | MC not dynamics; biased proposals for search |
| **CBMC / concerted rotation / ICB** (Dodd et al. 1993; Ulmschneider & Jorgensen; Kramer 2004) | Torsions (segment regrowth) sampled | Equilibrium sampling of chains/macrocycles | Piece 1 by explicit Jacobian / detailed-balance bookkeeping per move | Yes, with correct Jacobian | Same measure problem solved move-by-move; Robosample folds it into one det(M) via HMC |
| **Torsional diffusion (ML)** (Jing, Corso et al. 2022); Boltzmann generators (Noe et al. 2019) | Torsions only; bonds/angles rigid | Fast conformer generation / amortized sampling | Pieces 1+2 absorbed implicitly by learning; not from a constrained Hamiltonian | Approximate, or per-molecule-trained | Data-driven/amortized vs. physics-based/exact; complementary |
| **Robosample / gMolmodel (this work)** (Spiridon & Minh 2017; Spiridon, Sulea, Minh & Petrescu 2020) | Mixed worlds: Cartesian (covering) + torsional (Ball + Cyl planned) | Canonical sampling, multiscale | Piece 1: exact Fixman in MH acceptance; piece 2: coupled DOF mobilized in the flexible Cartesian world (Ball/Cyl planned); ring closures via RATTLE | **Yes** (asymptotically; conditions in Section 13) | -- |

Notes:

- The refinement/NMR tools and ICM are what most people mean by "torsional dynamics." Robosample
  targets the canonical ensemble, which pre-empts apples-to-oranges comparisons.
- gMolmodel (Spiridon & Minh 2017) and Robosample (Spiridon et al. 2020) are the same method and
  lineage: gMolmodel first cast constrained-dynamics HMC (CDHMC) as a Gibbs move; Robosample is the
  refactored application adding the spherical/cylindrical joints (Section 11) and the optional
  mass-matrix-determinant gradient (Fixman torque). The current code base reimplements the multibody
  engine independently of Simbody/Molmodel (Section 3).

---

## 15. References

Grounded in the two Robosample papers (Spiridon & Minh 2017; Spiridon et al. 2020) and their cited
sources; landscape rows in Section 14 also draw on the related-methods literature.

- Fixman, M. Classical statistical mechanics of constraints. *Proc. Natl. Acad. Sci. USA* 1974, 71, 3050.
- Fixman, M. Simulation of polymer dynamics. I. General theory. *J. Chem. Phys.* 1978, 69, 1527.
- van Gunsteren, W. F.; Berendsen, H. J. C. Algorithms for macromolecular dynamics and constraint
  dynamics. *Mol. Phys.* 1977, 34, 1311.
- Go, N.; Scheraga, H. On the use of classical statistical mechanics in the treatment of polymer
  chain conformation. *Macromolecules* 1976, 9, 535.
- Jain, A.; Vaidehi, N.; Rodriguez, G. A fast recursive algorithm for molecular dynamics simulation.
  *J. Comput. Phys.* 1993, 106, 258.
- Featherstone, R. *Rigid Body Dynamics Algorithms*; Springer, 2008. (articulated body inertia;
  the O(N) articulated-body recursion of Section 3.5)
- Rodriguez, G.; Jain, A.; Kreutz-Delgado, K. A spatial operator algebra for manipulator modeling and
  control. *Int. J. Robot. Res.* 1991, 10, 371.
- Jain, A. Compensating mass matrix potential for constrained molecular dynamics. *J. Comput. Phys.*
  1997, 136, 289.
- Vaidehi, N.; Jain, A.; Goddard, W. A. Constant-temperature constrained MD: the Newton-Euler inverse
  mass operator method. *J. Phys. Chem.* 1996, 100, 10508.
- Jain, A.; Park, I.-H.; Vaidehi, N. Equipartition principle for internal coordinate molecular
  dynamics. *J. Chem. Theory Comput.* 2012, 8, 2581.
- Wagner, J. R.; Balaraman, G. S.; Niesen, M. J.; Larsen, A. B.; Jain, A.; Vaidehi, N. Advanced
  techniques for constrained internal coordinate MD. *J. Comput. Chem.* 2013, 34, 904.
- Jain, A.; Kandel, S.; Wagner, J.; Larsen, A.; Vaidehi, N. Fixman compensating potential for general
  branched molecules. *J. Chem. Phys.* 2013, 139, 244103.
- Vaidehi, N.; Jain, A. Internal coordinate molecular dynamics: a foundation for multiscale dynamics.
  *J. Phys. Chem. B* 2015, 119, 1233.
- Kandel, S.; Salomon-Ferrer, R.; Larsen, A. B.; Jain, A.; Vaidehi, N. Overcoming potential energy
  distortions in constrained internal coordinate MD simulations. *J. Chem. Phys.* 2016, 144, 044112.
- Patriciu, A.; Chirikjian, G. S.; Pappu, R. V. Analysis of the conformational dependence of
  mass-metric tensor determinants in serial polymers with constraints. *J. Chem. Phys.* 2004, 121, 12708.
- Echenique, P.; Calvo, I.; Alonso, J. L. Quantum mechanical calculation of the effects of stiff and
  rigid constraints in the conformational equilibrium of the alanine dipeptide. *J. Comput. Chem.*
  2006, 27, 1733.
- Katritch, V.; Totrov, M.; Abagyan, R. ICFF: a new method to incorporate implicit flexibility into
  an internal coordinate force field. *J. Comput. Chem.* 2003, 24, 254.
- Chen, J.; Im, W.; Brooks, C. L. III. Application of torsion angle molecular dynamics for efficient
  sampling of protein conformations. *J. Comput. Chem.* 2005, 26, 1565.
- Rice, L. M.; Brunger, A. T. Torsion angle dynamics: reduced variable conformational sampling
  enhances crystallographic structure refinement. *Proteins* 1994, 19, 277.
- Stein, E. G.; Rice, L. M.; Brunger, A. T. Torsion-angle molecular dynamics as a new efficient tool
  for NMR structure calculation. *J. Magn. Reson.* 1997, 124, 154.
- Abagyan, R.; Totrov, M. Biased probability Monte Carlo conformational searches and electrostatic
  calculations for peptides and proteins. *J. Mol. Biol.* 1994, 235, 983.
- Forrest, B. M.; Suter, U. W. Generalized coordinate hybrid Monte Carlo. *Mol. Phys.* 1994, 82, 393.
- Duane, S.; Kennedy, A. D.; Pendleton, B. J.; Roweth, D. Hybrid Monte Carlo. *Phys. Lett. B* 1987,
  195, 216.
- Nilmeier, J. P.; Crooks, G. E.; Minh, D. D. L.; Chodera, J. D. Nonequilibrium candidate Monte
  Carlo is an efficient tool for equilibrium simulation. *Proc. Natl. Acad. Sci. USA* 2011, 108 (45),
  E1009. (NCMC; the decouple-move-recouple proposal and the deterministic-propagation acceptance of
  Section 5.8. Note D. D. L. Minh is shared with the Robosample lineage.)
- Chen, Y.; Roux, B. Efficient hybrid non-equilibrium molecular dynamics -- Monte Carlo simulations
  with symmetric momentum reversal. *J. Chem. Phys.* 2014, 141, 114107; and Generalized Metropolis
  acceptance criterion for hybrid non-equilibrium molecular dynamics -- Monte Carlo simulations.
  *J. Chem. Phys.* 2015, 142, 024101. (Hybrid neMD/MC; momentum-reversal acceptance referenced in
  Section 5.8.)
- Geman, S.; Geman, D. Stochastic relaxation, Gibbs distributions, and the Bayesian restoration of
  images. *IEEE Trans. Pattern Anal. Mach. Intell.* 1984, PAMI-6, 721.
- Chodera, J. D.; Shirts, M. R. Replica exchange and expanded ensemble simulations as Gibbs sampling.
  *J. Chem. Phys.* 2011, 135, 194110.
- Leimkuhler, B.; Reich, S. *Simulating Hamiltonian Dynamics*; Cambridge University Press, 2004.
  (symplectic integrators, backward error analysis, shadow Hamiltonian)
- Davidchack, R. L. Discretization errors in molecular dynamics simulations with deterministic and
  stochastic thermostats. arXiv:1412.7067, 2014. (h^2 discretization error, shadow Hamiltonian, and
  the onset of upward energy drift in rigid-body NVE beyond ~60-70% of the stability threshold)
- Ryckaert, J.-P.; Ciccotti, G.; Berendsen, H. J. C. Numerical integration of the Cartesian
  equations of motion of a system with constraints: molecular dynamics of n-alkanes (SHAKE).
  *J. Comput. Phys.* 1977, 23, 327.
- Andersen, H. C. Rattle: a "velocity" version of the SHAKE algorithm for molecular dynamics
  calculations. *J. Comput. Phys.* 1983, 52, 24.
- Sherman, M. A.; Seth, A.; Delp, S. L. Simbody: multibody dynamics for biomedical research.
  *Procedia IUTAM* 2011, 2, 241. (algorithmic lineage; not a dependency of the current code)
- Flores, S. C.; Sherman, M. A.; Bruns, C. M.; Eastman, P.; Altman, R. B. Fast flexible modeling of
  RNA structure using internal coordinates (Molmodel). *IEEE/ACM Trans. Comput. Biol. Bioinf.* 2011,
  8, 1247. (algorithmic lineage; not a dependency of the current code)
- Eastman, P.; et al. OpenMM 7: rapid development of high-performance algorithms for molecular
  dynamics. *PLoS Comput. Biol.* 2017, 13, e1005659.
- Spiridon, L.; Minh, D. D. L. Hamiltonian Monte Carlo with constrained molecular dynamics as Gibbs
  sampling. *J. Chem. Theory Comput.* 2017, 13 (10), 4649.
- Spiridon, L.; Sulea, T. A.; Minh, D. D. L.; Petrescu, A.-J. Robosample: a rigid-body molecular
  simulation program based on robot mechanics. *Biochim. Biophys. Acta Gen. Subj.* 2020, 1864 (8), 129616.
- Dodd, L. R.; Boone, T. D.; Theodorou, D. N. A concerted rotation algorithm for atomistic Monte
  Carlo sampling of polymer melts. *Mol. Phys.* 1993, 78, 961.
- Kramer, A. Thermodynamic sampling of molecular conformations. arXiv:physics/0401036, 2004.
- Jing, B.; Corso, G.; Chang, J.; Barzilay, R.; Jaakkola, T. Torsional diffusion for molecular
  conformer generation. *NeurIPS* 2022, arXiv:2206.01729.
- Noe, F.; Olsson, S.; Kohler, J.; Wu, H. Boltzmann generators. *Science* 2019, 365, eaaw1147.
