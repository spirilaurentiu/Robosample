# Robosample

## Coding rules

This is scientific software that requires special reasoning and coding constraints.

### Rule 0 - Adhere to theoretical background

This document explains scientific definitions, motivation, reasoning and default behavior.
A change may deliberately depart from it, but only if the departure is argued and its impact documented.

### Rule 1 - Think Before Coding

State **ALL** assumptions explicitly. Ask rather than guess.
Break the implementation into logical, bite-sized tasks or milestones.
Push back when a simpler approach exists. Stop when confused.

### Rule 2 - Simplicity First

Minimum code that solves the problem. Nothing speculative.
No abstractions for single-use code.

### Rule 3 - Surgical Changes

Touch only what you must. Don't improve adjacent code.
Match existing style. Don't refactor what isn't broken.

### Rule 4 - Goal-Driven Execution

Define success criteria. Loop until verified.
If in doubt, iteratively ask questions until requirements and edge cases are covered.
Strong success criteria let Claude loop independently.

### Rule 5 - Use the model only for judgment calls

Use for: classification, drafting, summarization, extraction.
Do NOT use for: routing, retries, deterministic transforms.
If code can answer, code answers.

### Rule 6 - Surface conflicts, don't average them

If two patterns contradict, pick one (more recent / more tested).
Explain why. Flag the other for cleanup.

### Rule 7 - Read before you write

Before adding code, read exports, immediate callers, shared utilities.
If unsure why existing code is structured a certain way, ask.

### Rule 8 - Tests verify intent, not just behavior

Tests must encode **WHY** behavior matters, not just WHAT it does.
A test that can't fail when business logic changes is wrong.
Output is valid only if you can name the theory your change **preserves** or **disrupts**.

### Rule 9 - Checkpoint after every significant step

Summarize what was done, what's verified, what's left.
Don't continue from a state you can't describe back.
Determinism is a feature, not a side effect.

### Rule 10 - Conventions in this codebase are load-bearing

Conformance > taste inside the codebase.
If you think a convention is harmful, surface it. Don't fork silently.
For instance, Frame conventions (`F` vs `M`, parent-frame vs body-frame angular velocity), index orderings (angular over linear in spatial vectors), sign conventions (`Phi` vs `~Phi` for forces vs velocities) are not stylistic.

### Rule 11 - Fail loud

"Completed" is wrong if anything was skipped silently.
"Tests pass" is wrong if any were skipped.
Always surface uncertainty.
Defaults that can be overridden when justified.

### Rule 12 - Self reflect

Review your changes before finishing the answer.
Act as your own hostile reviewer.

### Rule 13 - Programming languages, libraries and hardware

- C++17 on GCC 12.0 - 15.0 compiled with `x86_64-v3` (`SSE4.2`, `POPCNT`, `AVX`, `AVX2`, `BMI1` / `BMI2`, `FMA`, `LZCNT`, `MOVBE`) and LTO.
- Python 3.12
- CUDA 12.0 - 13.2
- OpenMP
- OpenBLAS and LAPACK available. CUDA counterparts available on CUDA capable devices (virtually all targeted machines)
- OpenMM 8.5 compiled directly into our source tree.

Code will run on consumer grade hardware (gamer PCs).

## Scientific background

### What problem are we solving?

Suppose a molecular system has coordinates `x`. The equilibrium probability of observing a configuration is `pi(x) proportional to e^-(beta * U(x))` where `U` is the potential energy and `beta = 1/(R*T)` is the Boltzmann constant. The entire purpose of the software is to generate a sequence of molecular configurations whose long-term distribution is `pi(x)`. Concretely, we are trying to maximize **basin hopping** and **intra-basin exploration**.

### Molecular system size

Up to 1M atoms clustered in up to 100k rigid bodies. Examples:

- Alanine dipeptide in vacuum
- Deca alanine in explicit solvent
- FFAR1 (GPCR) in implicit solvent, but with explicit membrane nanodisc
- Spliceosome in explicit solvent

## Definitions for default behavior

### **HAMILTONIAN GEOMETRY**

The space of all possible states `(p,q)` (position and momentum) is called the phase space. The Hamiltonian `Hamiltonian(q,p)` is the total energy written as a function of position and momentum and it dictates the complete motion through 2 equations:

- `dq/dt =  d(H)/dp`: Position is steered by how energy varies with momentum.
- `dp/dt = -d(H)/dq`: Momentum by how energy varies with position, but with a flipped sign.

The system does not move towards lower energy, but stays along iso-Hamiltonian levels. At every phase point, the symplectic form maps the energy gradient to a flow that conserves `H`.

### **MARKOV CHAIN**

The sampler produces `s0 -> s1 -> s2 -> s3 -> s...`. The next state depends on the current one: `prob(s(n+1)|s(n),s(n-1),...) = prob(s(n+1)|s(n))`. This is the Markov property. The Markov chain is the central object of the entire theory.

### **STATIONARY DISTRIBUTION**

A distribution `pi` is stationary if one step leaves it unchanged: `pi*TransitionKernel = pi`. If states are already distributed according to `pi`, another sampling step still gives `pi`.

### **DETAILED BALANCE**

Let `a` and `b` be two states and `TransitionKernel` the transition kernel between them. Then, `pi(a) * TransitionKernel(a->b) = pi(b) * TransitionKernel(b->a)`. Interpretation: The probability current flowing from `a` to `b` equals the current flowing backward. Detailed balance implies stationarity. Stationarity does not require detailed balance.
Example: let there be two conformations trans and gauche. At equilibrium:

- 100 transitions/day trans -> gauche
- 100 transitions/day gauche -> trans

The populations remain constant.

### **ERGODICITY**

The chain is ergodic if, starting from any `s0`, trajectory averages converge to expectations under `pi`: `(1/N) * sum f(s_n) -> integral f(x) * pi(x) * dx` where `pi` is defined over all degrees of freedom (not a subspace). Stationarity says the kernel preserves `pi`. Ergodicity says the chain actually reaches `pi`. Both are needed: preservation without reachability is useless. One-line counterexample: the identity kernel `TransitionKernel(s,s)=1` preserves every distribution but never moves. Stationary, not ergodic. In practice, the failure is reducibility: barriers the proposal cannot cross in the run budget. The chain samples one basin correctly and ignores the others. Example: deca-alanine has an alpha-helical basin and an extended beta/PPII basin separated by several kcal/mol. With short-trajectory Cartesian HMC at 300 K:

- 1000 intra-basin transitions/day in alpha
- 1000 intra-basin transitions/day in beta
- 0 inter-basin transitions/day
Detailed balance holds locally. Starting in alpha, the beta basin contributes zero weight instead of its true Boltzmann weight. The chain is non-ergodic with respect to the full pi.

### **HAMILTONIAN**

Let the Hamiltonian energy function be `Hamiltonian(q,p) = U(q) + KE(q,p) + U_Fixman(q) + U_Jacobian(q)`.

For generalized coordinates:

- `U(q)` is the potential energy computed on full-atom model.
- `KE(q,p) = 1/2* p^T *M(q)^-1* p` is the kinetic energy.
- `U_Fixman(q)` is a correction term (see below) and is null in Cartesian space.
- `U_Jacobian(q)` is `ln( sin^2(gamma2) )` summed over all `Free`-rooted molecules. Note that `gamma2` is the pitch of each root's absolute orientation. This term is null in Cartesian space.

For cartesian coordinates, `KE = 1/2 * p^T * M^-1 * p` where `p^T` is the transpose momenta vector.
The Hamiltonian determines the equilibrium distribution: `pi(q,p) proportional to e^(-beta*H(q,p))`

### **ACCEPTANCE**

A proposal is a candidate move: `s -> s'`. Example: draw momenta, integrate dynamics for 100 steps. The proposal is not yet accepted.

The Metropolis rule: `A=min(1,e^(-beta*deltaH))`. If `deltaH == 0`, then accept with probability `1`. If energy increased, accept probabilistically. Acceptance restores detailed balance.

Simple HMC example: draw momenta, integrate Hamilton's equations, compute `deltaH`, accept or reject. This produces a Markov kernel `T(s->s')`

### **REVERSIBILITY**

A proposal is reversible if the reverse path exists:

- Forward: `(q, p)` -> `(q', p')`.
- Reverse: `(q',-p)` -> `(q, -p)`.

Without reversibility, the backward transition probability cannot be computed.

### **VOLUME PRESERVATION**

The proposal must preserve phase-space volume `(dq,dp)`. No compression. No expansion. Otherwise the proposal density changes. The Jacobian would appear in the acceptance probability.

### **SYMPLECTICITY**

Symplectic integrators preserve Hamiltonian geometry. Consequences: phase-space volume preservation, near-energy conservation, shadow Hamiltonian. Energy errors remain bounded: `H(t) - H(0) = O(dt^2)`. This gives high acceptance rates. Symplecticity improves efficiency. It is not strictly required for correctness.

### **GIBBS SAMPLING**

A Gibbs step updates only part of the variables. Examples:

- `q=(q1, q2, q3)`. Update only `q2`. The other coordinates remain fixed.
- torsional dynamics: only torsions move; bond lengths and bond angles are fixed. This constitutes a torsional Gibbs block.

A Gibbs block defines: mobile coordinates, frozen coordinates, mass matrix, Hamiltonian, timestep. Each block is itself an HMC sampler.

### **COMPOSITION**

If `TransitionKernel1` and `TransitionKernel2` are transition kernels which preserve `pi`, then `TransitionKernel2 * TransitionKernel1` also preserves `pi`. This is why Cartesian blocks and torsional blocks can be mixed without introducing biases.

### **GENERALIZED COORDINATES**

Instead of Cartesian coordinates `(x,y,z)`, we use `q=(bond_length, bond_angle, torsion/dihedral)`. The mass matrix then becomes `M(q)`. Based on the accessible degrees of freedom, the space can be:

- Flat (linear Taylor):
  - `Pin` torsions live on the circle `S^1`, which has **zero intrinsic curvature** (a circle is a line made periodic, locally indistinguishable from `R`) and is integrated as an unwrapped real.
  - `Translation` and `Free`  coordinates lie in the flat, zero-curvature space `R^3`.
- Curved path (exponential map): The **orientation** of a `Free` or `Ball` body, stored as a unit quaternion living on **S^3, the unit 3-sphere in R^4** which is a compact manifold of **constant positive curvature** that **double-covers the rotation group SO(3)** (`q` and `-q` are the same physical rotation).

From a physics point of view, bonds and angles vary rapidly and contribute little to the overall RMSD. On the contrary, torsions are softer modes and heavily influence conformational transitions, but these transitions are dependent on all degrees of freedom. For example, torsional dynamics cannot sample alone cis to trans isomerization of alanine dipeptide which is gated by a 1-4 clash that requires angle relaxation. The transition occurs only when this torsions are sampled together with bond lengths.

### **MASS METRIC**

The generalized mass matrix acts as a metric `ds^2 = dq^T * M(q) * dq`. Nearby coordinate changes can have different physical meanings depending on configuration. This is the origin of: Fixman correction to correct for marginalizing over momenta, generalized equipartition and velocity-dependent forces.

### **EQUIPARTITION**

The correct momentum distribution is `p ~ N(0, R*T*M(q))`, not `p ~ N(0, I)`. Otherwise the kinetic temperature is wrong.

### **FIXMAN POTENTIAL**

The coordinate transformation introduces a metric bias. The correction is `U_Fixman(q) = 1/2 * R * T * ln(det(M(q)))`. It removes the sampling distortion caused by generalized coordinates.

Generalized-coordinate (e.g. torsional) sampling is designed so that the marginal distribution of configurations matches the marginal of the **fully-flexible Cartesian** Boltzmann distribution restricted to the sampled subspace; the Fixman term removes the mass-metric artifact that would otherwise distort this marginal, and the mixed Gibbs scan supplies the relaxation of the frozen coordinates.

### **CONSTRAINTS**

Constraints define a manifold `sigma(q)=0`. Examples: ring closure. SHAKE projects positions. RATTLE projects velocities.

## Internal-coordinate representation (robot model)

Each molecule is represented as a **robot** (kinematic tree of rigid bodies connected by joints, or *mobilizers*). For multiple molecules, this becomes a forest rooted at a shared single ground frame.

An **O(n) articulated-body recursion** -- Featherstone's articulated-body algorithm, in spatial-operator form -- propagates positions, velocities, accelerations, and forces up and down the tree without forming or inverting the mass matrix M. **Per-atom** Cartesian forces (`-grad(U)`) computed by OpenMM are reduced to per-body spatial forces (net force + torque about the body origin) for the articulated-body solver.

The recursion is realized in a custom, data-oriented (structure-of-arrays) multibody engine; the global `M`, `M^-1`, `sqrt(M)`, and `det M` are never assembled, only applied as `O(n)` operators.

## Z-matrix

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

## The mobilizer frame chain (Simbody convention)

Each mobilized body `B` is positioned by a chain of frames between ground (`G`) and `B`. Following the Simbody convention, the frames are:

- **Frame_G** -- Ground (inertial) frame.
- **Frame_P** -- parent body frame (the inboard body's frame).
- **Frame_F** -- the inboard *fixed* mobilizer frame, rigidly attached to the parent `P`.
- **Frame_M** -- the outboard *moving* mobilizer frame, rigidly attached to the child body `B`.
- **Frame_B** -- the child body frame, anchored at the body's root (inboard) atom; every atom of the body is stored as a constant *station* `r_a` in `B`.

Two transforms are constant (built once at setup from the Z-matrix geometry) and one carries the joint motion:

- `X_PF` : `Frame_P` -> `Frame_F`, the inboard frame's placement on the parent (constant)
- `X_BM` : `Frame_B` -> `Frame_M`, the outboard frame's placement on the child body (constant)
- `X_FM(q)` : `Frame_F` -> `Frame_M`, the mobilizer transform, the ONLY configuration-dependent factor (joint DOF)

The body's pose in Ground is the ordered composition `X_GB = X_GP * X_PF * X_FM(q) * X_MB` with `X_MB = (X_BM)^-1`.

Atom positions follow as `r_a(Ground) = X_GB . r_a`. Only `X_FM(q)` changes as `q` moves; `X_PF` and `X_BM` are fixed at construction. For a `Pin` (torsion) joint the rotation is about the mobilizer Z axis, so the setup bakes a fixed bond->Z alignment (a -90 deg rotation about Y in this engine) into `X_PF`/`X_BM`, placing Z along the rotated bond.

For joints with a quaternion-parameterized rotation (`Free` and `Ball`), the four orientation coordinates are the quaternion `q = (qw, qx, qy, qz)` representing `R_FM`, and the three rotational generalized speeds are the angular velocity `w_FM = (wx, wy, wz)` of `Frame_M` relative to `Frame_F` **expressed in the parent F frame** (Simbody convention). The rotation map, the quaternion derivative, its inverse, and the second derivative are exactly Simbody's (Rotation.h: calcUnnormalizedNForQuaternion / NInv):

```text
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

Component form of the derivative -- the form the integrator and the SHAKE `q`-update must use:

```text
qdot_w = (1/2) * ( -qx*wx - qy*wy - qz*wz )
qdot_x = (1/2) * (  qw*wx + qz*wy - qy*wz )
qdot_y = (1/2) * ( -qz*wx + qw*wy + qx*wz )
qdot_z = (1/2) * (  qy*wx - qx*wy + qw*wz )
```

This is the parent-frame map. The body-frame map `qdot = (1/2) q (x) (0, w)` negates the off-diagonal (cross-coupling) terms of `N` -- giving `qdot_x = (1/2)(qw*wx - qz*wy + qy*wz)`, and likewise for `qdot_y`, `qdot_z` -- and expects `w` expressed in `Frame_M`. It **must not** be used here, because `w_FM` is expressed in `F` and `R_FM(q)` above is the standard map: pairing the standard `R_FM` with the body-frame `N` integrates the orientation with the wrong-handed angular velocity.

## Construction of the mobilizer frames from bond geometry

`X_FM(q)` is the joint DOF (`Rot(q, Z)` for a `Pin`, identity for a `Weld`, `Translate(q)` for a Translation, quaternion+translation for a `Free`). The two **constant** frames `X_PF` and `X_BM` are rebuilt once per configuration transfer from the actual Cartesian geometry, in two stages.

### Stage 1 -- a frame per atom, from bond vectors

Each atom a gets a Ground-frame pose `frame[a]` and a cross-bond transform `xpc[a]`, built from four atoms: the atom (self), its parent, its
grandparent, and a reference child (the four atoms that define the atom's internal coordinates). With `pSelf`, `pParent`, `pGparent`, `pRefChild` their Cartesian positions:

```text
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

So the **bond vector is the local X axis**, the **grandparent fixes the perpendicular** (the (self,parent,gparent) plane; its normal is Z), the **dihedral** rides on X, and the **bond length** is the X-translation. The two 180-degree flips (about Y inside B, then about X) reproduce the bond-center "about-face" convention of the original Molmodel build -- the parent's outboard bond-center direction and the child's inboard one point in opposite senses along the bond -- so `frame[a]` matches the original bond-center frame protocol to ~1e-14. That equivalence is exactly what the angle/chirality reconstruction cross-check verifies: rebuilding each atom's local bond-center directions from bond angles + chirality and comparing against the position-derived frame localizes
any mismatch to the offending atom. A **root** atom gets frame = (identity rotation, pSelf), xpc = identity. An atom lacking a full (parent, grandparent) lineage falls back to X = bond direction with xpc = Translate(d, X) * Rot(180 deg, Y) and no dihedral.

### Stage 2 -- mobilizer frames from the per-atom frame

For body b with root (inboard) atom `root` and parent body p with root atom `proot`:

```text
T_X_B = frame[root]                                          // the body frame B == its root atom's frame
if p is Ground:
    X_PF[b] = T_X_B                                          // Frame_F is the root-atom frame, in Ground
    X_BM[b] = identity                                       // Frame_M coincides with B
else:
    Proot_X_root = inverse(frame[proot]) * frame[root]       // pose of B's root atom in p's root-atom frame
    X_BM[b] = xpc[root] * X_to_Z                             // Frame_B -> Frame_M
    X_PF[b] = Proot_X_root * X_BM[b]                         // Frame_P -> Frame_F
where   X_to_Z = Rot(-90 deg, Y)                             // carries the bond (X) onto the mobilizer Z axis
```

`X_to_Z` is the only joint-convention-specific piece: a `Pin` rotates about the mobilizer **Z** axis, but Stage 1 puts every bond on **X**, so the -90-degree rotation about Y carries X onto Z (X -> +Z, Z -> -X). The `Frame_M` frame's Z axis therefore lies along the inboard bond, and `X_FM(q) = Rot(q, Z)` turns the body about that bond. The reference dihedral is baked into `xpc[root]` (hence into `X_BM` and `X_PF`), so at `q=0` the chain reproduces the **incoming** geometry exactly, and the `Pin` coordinate `q` is the torsion measured from that carried-over reference.

With these definitions and `X_FM(0) = Identity`, the composition must return the body frame:

```text
X_GB(0) = X_GP * X_PF * X_FM(0) * X_MB
        = frame[proot] * (Proot_X_root * X_BM) * I * (X_BM)^-1
        = frame[proot] * Proot_X_root
        = frame[proot] * inverse(frame[proot]) * frame[root]
        = frame[root] = T_X_B.
```

Any change to the frame build must preserve this identity; comparing X_GB at q = 0 against frame[root] for every body is the fastest regression test, and -- unlike an energy check -- it catches the X_MB-vs-X_BM and bond->Z sign errors.

The three factors of the chain are recomputed on deliberately different schedules, and that split is the point of the design:

- `X_PF (P->F), X_BM (B->M, hence X_MB = ~X_BM)`: once per Gibbs-block transfer, plus once at world construction. Role: the inboard bond's frozen geometry (length + angle + reference dihedral); held constant for the block's entire MD trajectory.
- `X_FM(q)`: every MD step. Role: the live joint DOF; for a `Pin`, `Rot(q, Z)` with the current torsion `q`.
- `X_GP` (the parent body's X_GB) and `X_GB` itself: every MD step. Role: the running composition; `X_GP` is constant only for a root body (parent = Ground)

So `X_PF` and `X_BM` carry the per-block **frozen** geometry while `X_FM`, `X_GP`, and `X_GB` carry the **per-step** motion. **The same transfer that rebuilds `X_PF` and X_BM `also` resets `q=0`**, making the incoming configuration the block's reference. Within a block only `X_FM(q)` and the ancestor poses move, so `X_GB(q)` is the only thing the integrator recomputes each step.

The frame graph that Stage 1 consumes is filled from the molecular bond graph with ring-closing bonds removed (parent = the bonded atom toward the root,
grandparent = its parent, refChild = the atom's first child). That acyclic tree is the same BAT / Z-matrix internal-coordinate tree the topology builder produces as `(z_i, z_j, z_k, z_l)`. The correspondence is:

- Tree topology: each atom's parent `z_j` and grandparent `z_k` is identical to the Z-matrix.
- The bond length `d = |self - parent|` and the bond angle (the grandparent fixes the frame's perpendicular) are the Z-matrix bond/angle internal coordinates, and they are frozen into `X_PF` and `X_BM`.
- Torsions are the live coordinates: a flexible dihedral about a tree bond becomes a `Pin` joint, so a Z-matrix torsion is `q = X_FM(q)`. Its reference value is baked into `X_BM` via `xpc[root]` (which carries `Rot(theta, X)`), and `X_FM(q) = Rot(q, Z)` turns about that bond.

Two qualifications are load-bearing and must not be relaxed into an identity. First, the frame kernel derives this tree and all geometry from the Cartesian positions; it shares the *tree* with the separate Z-matrix builder, not the code, and it uses the atom's first child as the local dihedral reference (pointing away from the root), whereas a textbook Z-matrix `z_l` references toward the root -- both describe the torsion about the same bond, but with a different fourth atom. Second, the frozen bond lengths and angles are not hard-constrained to idealized values: they sit at whatever the incoming geometry holds and are merely held fixed for the block (Gibbs conditioning), never a permanent constraint.

For body `b` with root atom `root = bodyRootAtom[b]` (the code's "inboard (root) atom of this body"), the inboard bond is the tree edge `(par[root], root)`:

- The **inboard** atom `A = par[root]` lives in the parent body. The `F` frame `X_PF` is fixed to the parent at this joint.
- The **outboard** atom `B = root` is the child body's root atom. The `M` frame `X_BM` and the body frame `B` are anchored here.

Which physical atom is inboard versus outboard is set per molecule by the Z-matrix root, not by any fixed convention such as N-before-CA: for a given molecule it is determined, but it must be read off that molecule's tree rather than assumed.

## The multibody dynamics engine (articulated-body recursions)

Every kinematic and mass-matrix quantity is computed by `O(n)` recursions over the n bodies of the tree via **Featherstone's articulated-body algorithm**, in spatial-operator form. The global mass matrix `M`, its inverse `M^-1`, its square root, and its determinant are **never assembled as dense matrices**. Each is applied as an operator by one or two sweeps. This linear cost is the reason the internal-coordinate representation exists, and it is what makes per-block sampling scale to large systems.

Each body `b` carries, in the Ground frame, a spatial velocity `V_GB = (w ; v)` (angular **over** linear), a spatial force `F = (torque ; force)` (same ordering), and a spatial inertia `Mk_G[b]` (the body's rigid-body spatial inertia about its own origin `Bo`, re-expressed in Ground; mass, first moment `mass*c`, and the inertia tensor). The spatial inner product is the straight component pairing `< V, F > = w * torque + v * force`, so a generalized force is `tau = H^T F` under that pairing. Four operators relate a body to its parent p:

- **Phi (rigid shift / transmission).** With `r=p_PB` the vector from the parent origin `Po` to the body origin `Bo` (in Ground), `Phi` is the 6x6 built from 3x3 blocks `[[1, r_x],[0, 1]]`, `r_x` the cross-product matrix of `r`. `~Phi` (its transpose) shifts velocities and accelerations **outward** (parent -> child); `Phi` shifts forces and inertias **inward** (child -> parent). It carries no mass - pure geometry.
- **H (hinge / joint map).** Maps a body's dof generalized speeds to its cross-mobilizer spatial velocity, `V_PB_G = H*u`. `H` is built in the `F` frame and rotated to Ground.
- **P (articulated body inertia).** The apparent spatial inertia a body exhibits as the free base of the subtree hanging below it. Unlike a rigid-body inertia, it is a full symmetric 6x6 (the apparent mass depends on the direction of push, and there is no single center of mass).
- **D, Ga.** The mobility-space (hinge) inertia `D = H^T * P * H` (dof x dof), its inverse `D^-1`, and the articulated **gain** `Ga = P * H * D^-1` assembled down the tree. The symbol `Ga` is deliberately distinct from the constraint Jacobian `G` and the kinetic energy `KE`. The inertia the parent feels across the mobilizer is the articulated shift `P^+ = P - Ga * D * Ga^T = P - P * H * D^-1 * H^T * P`.

**The recursions -- each a single sweep, O(N):**

1. **Position -- outward (root -> leaves).** Compose the frame chain `X_GB = X_GP * X_PF * X_FM(q) * X_MB`. Build `H` in Ground, the shift `Phi`, the body spatial inertia in Ground, and the per-atom stations `r_a(Ground) = X_GB * r_a`.
2. **Velocity -- outward.** Shift the parent's spatial velocity and add the mobilizer term `V_GB[b] = (~Phi[b]) * V_GB[parent] + H[b] * u[b]`.
   The same sweep accumulates the two velocity-dependent bias terms the integrator needs. The gyroscopic force (Ground frame; `I` the unit inertia about `Bo`, `c` the COM offset, `w=w_GB`): `b_gyro[b] = mass * ( w x (I w) ;  w x (w x c) )`.

   and the Coriolis bias acceleration, carried outward like the velocity (`w_GP`, `v_GP` the parent spatial velocity; `w_PB_G`, `v_PB_G` the cross-mobilizer velocity; `w_FM` the cross-mobilizer angular velocity in `F`; `r_MB_F = R_FM * r_MB` the `Mo` -> `Bo` bond vector expressed in `F`):

   ```text
   a_mob[b] = ( w_GP x w_PB_G ;
                w_GP x (v_GB - v_GP) + w_GP x v_PB_G + R_GF * ( w_FM x (w_FM x r_MB_F) ) )
   a_cor[b] = (~Phi[b]) * a_cor[parent] + a_mob[b]
   ```

   The centripetal term `R_GF (w_FM x (w_FM x r_MB_F))` is nonzero for every body whose `M` frame is offset from `Bo` (every non-root `Pin` body, offset = bond length); dropping it silently pumps energy. It is the velocity-side companion of the `X_MB` pitfall -- `r_MB = X_MB.p()`, never `X_BM.p()`.
3. **Articulated body inertia -- inward (leaves -> root).** Featherstone's backward pass:

   ```text
   P[b]   = Mk_G[b] + sum_{children c} Phi[c] * P^+[c] * (~Phi[c])
   D[b]   = H[b]^T P[b] H[b],   Ga[b] = P[b] H[b] D[b]^-1,   P^+[b] = P[b] - Ga[b] D[b] Ga[b]^T
   ```

   A `Weld` (0 dof) contributes nothing through its (empty) mobility space, so `P^+ = P`.
4. **Forward dynamics (M^-1) -- inward then outward.** Given a generalized force `f`, the acceleration `udot = M^-1 (f - f_bias)` is produced by an inward residual pass (the per-body bias `Z` and the mobility residual `eps = f - H^T * Z`) followed by an outward acceleration pass (`A_GB` and `udot`), using `D^-1` and `Ga`. The bias force `f_bias` is assembled from `b_gyro` and `a_cor` of recursion 2 (the gyroscopic force plus the Coriolis acceleration mapped through the articulated inertia): with that `f_bias` this is the equation of motion solved each Verlet substep; with `f_bias = 0` it is the pure operator `a = M^-1 * f` used for the constraint solver and the kinetic energy.
5. **Force -- inward.** The transpose of forward kinematics: per-atom Cartesian forces are reduced to a generalized force by accumulating each body's spatial force and sweeping child -> parent through `Phi^T`, reading `tau = H^T` (transmitted spatial force). This is exactly the `G^T` assembly the constraint solver uses.

**Mass-matrix operators -- all O(n), all matrix-free.** The same per-body factors `{Phi, H, D, Ga}`
yield, without ever forming M:

- **M^-1 * f:** Forward dynamics (recursion 4 with `f_bias = 0`).
- **sqrt(M^-1) z:** Outward sweep using the per-body `D^-1/2` in place of `D^-1`. It maps a white-noise vector `z` to a generalized-speed vector with covariance `M^-1`, so the equipartition draw is realized as `u = sqrt(R*T) * sqrt(M^-1) * z`, giving `u ~ N(0, R*T*M^-1)` and hence `p = M u ~ N(0, R*T*M)` at O(N), with det `M` never appearing in the proposal.
- **ln(det(M)) = sum_b ln(det(D_b)):** The factorization makes the determinant the product of the per-body hinge determinants - the Fixman tree term, read from the same `D_b = H_b^T * P_b * H_b`.
- **(1/2) u^T M u:** Kinetic energy, from one application of `M`.

**The operator factorization.** In spatial-operator form the mass matrix factors through the same local quantities"

```text
M       = (I + H * Phi * Ga) * D * (I + H * Phi * Ga)^T
M^(1/2) = (I + H * Phi * Ga) D^(1/2)
v       = D^(1/2) (I + H * Phi * Ga)^T * theta_dot       (the whitening map)
```

with `Ga` the gain operator and `v` the whitened velocity. Because each factor is local (a body and its parent only) `M`, `M^-1`, `M^(1/2)`, and
`ln det(M)` all follow from one or two sweeps over the tree, which is the O(n) claim. The generalized equipartition principle is then exact by construction: drawing the whitened `v` as unit Gaussian noise and applying `sqrt(M^-1)` gives `< u u^T > = R * T * M^-1`, i.e. equipartition of kinetic energy across the reduced coordinates despite the configuration-dependent metric.

## The Gibbs scan

Blocks are executed in a chosen order (a deterministic systematic scan) or randomized (random scan). The order is a free choice, not a constraint of the method: each block individually leaves `pi` invariant, and a composition of `pi`-invariant kernels leaves `pi` invariant regardless of order. A systematic scan is generally not reversible (its time-reversal is the reversed block order). A random scan that picks a block from a fixed distribution is reversible.  Either way `pi` is preserved. **Stationarity (`pi*TransitionKernel = pi`), not reversibility, is the requirement**.

Although torsional dynamics is the most efficient way to exploit the internal-coordinate representation -- and the reason the representation exists -- it is **not** mandatory; a scan may mix Cartesian, torsional, and other generalized-coordinate blocks freely.

- **Overlap.** The restricted spaces of different blocks may overlap. Redundant updates do not break stationarity, provided each block individually leaves pi invariant.
- **Ergodicity.** Across one full round (all blocks), every mobile coordinate in the block set has been **offered** a move. Ergodicity is over the **sampled (restricted) subspace** and requires both that the blocks collectively cover that subspace and that each block achieves nonzero acceptance. Being offered a move is **not** the same as accepting one: a block stuck near zero acceptance contributes nothing even though aggregate statistics look healthy. Coordinates that are frozen in **every** block (e.g. bonds and angles in torsion-only operation) are not sampled; the sampled distribution is then the (Fixman-corrected) marginal over the union of the blocks' restricted spaces.
- **Per-block diagnostics.** Acceptance rates are reported **per block** specifically to detect an under-sampled subspace, which would otherwise be masked by healthy global energy and acceptance statistics.

## The HMC move inside a block

### Momentum resampling and the internal-coordinate equipartition principle

At the start of each block the momenta are **fully** resampled (not partial): `p ~ N(0, R*T*M(q))`. Equivalently, `u ~ N(0, R*T*M(q)^-1)`.

This is a Gibbs step on the momenta (always accepted). It realizes the **internal-coordinate equipartition principle** : in Cartesian coordinates every DOF carries `(1/2) * R * T` of kinetic energy, but in generalized coordinates with a configuration-dependent metric `M(q)` this does not hold per coordinate. The *generalized* equipartition theorem, `< p_i (dH/dp_j) > = R * T * delta_ij`, still holds, which forces the momentum covariance to be`< p p^T > = R * T * M(q)`. Drawing `p ~ N(0, R*T*M(q))` is exactly this distribution, so temperature is assigned consistently despite the non-constant metric. The draw is matrix-free via the articulated `sqrt(M(q)^-1)` operator (drawing velocities; momenta follow as `p = M(q)*u`), so `det M(q)` never appears explicitly in the proposal; its configuration dependence is bookkept by the Fixman term in the acceptance.

In a ring-closure block the freshly drawn momenta are then RATTLE-projected onto the velocity constraint surface `G * M^-1 * p = 0` before integration begins, so the initial velocity is consistent with the closed ring. That projection is exactly what injects the `det(G M^-1 G^T)^(-1/2)` factor into the configurational marginal, the factor the loop-closure term of `F` cancels. For acyclic molecules there is no projection and no such factor.

### Propagation under the bare potential `U`

Dynamics are integrated with the **fixed-step Verlet integrator** (internal and Cartesian alike), using forces derived from the bare potential `U` only. Concretely, the integrator's forces are the per-atom Cartesian forces `-grad(U)` from OpenMM, reduced to per-body spatial forces for the articulated solver; the gradients of the Fixman potential `U_Fixman` and the Jacobian `U_Jacobian` (the "torques") are not added to these forces. So `U` is what generates the trajectory; `U_Fixman` and `U_Jacobian` enter only the acceptance. After each step, constraints (if any) are projected onto the manifold.

The integrator is a second-order, semi-explicit predictor-corrector, not an explicit kick-drift-kick leapfrog. This distinction is load-bearing and is the source of the integrator's actual conservation properties: the position is advanced by an explicit second-order Taylor step using the start-of-step acceleration, and the velocity is advanced by an **implicit trapezoidal corrector** solved by functional iteration. The two coincide with textbook velocity Verlet only for a separable, constant-mass Hamiltonian (the Cartesian case); for the configuration-dependent metric `M(q)` and velocity-dependent (Coriolis/gyroscopic) forces of internal coordinates they do not, which is exactly where the symplecticity caveat applies.

The integrator is run with a **fixed step**, making it symplectic and reversible for separable systems and the step is **taken unconditionally** since there is no step-size adaptation. Thus, a non-converged corrector does not shrink `dt` or reject the step. It is the trajectory-level Metropolis test, not per-step control, that supplies correctness.

The per-step sequence (one fixed-step Verlet step of size `dt`; the SHAKE/RATTLE projections are no-ops on acyclic molecules). `qdot0 = N(q0) * u0` for quaternion DOF and `u0` otherwise; `qDotDot(q)` is the corresponding q-acceleration; `udot = M^-1 * (f - f_bias)` is one O(n) forward-dynamics sweep, with `f_bias` carrying the gyroscopic and Coriolis terms:

#### 1. Position (2nd-order Taylor expansion)

- `q <- q0 + dt*qdot0 + (dt^2/2)*qdotdot0` - not applied to quaternions.

---

#### 2. Quaternion Rotation (exponential map update)

- `q <- advanceQuatExp(wHalf, dt) (x) q0`
- where `wHalf = u0 + (dt/2)*udot0`
- then normalize:
  - `q <- q / |q|`

---

#### 3. SHAKE (position constraint projection)

- Project positions to satisfy constraints `sigma(q) = 0`
- Followed by kinematic refresh after projection

---

#### 4. Velocity Predictor (forward Euler step)

- `u <- u0 + dt*udot0` where forces are evaluated from `-grad(U) at q1`

---

#### 5. Velocity Corrector (trapezoidal iteration)

- Iterate up to 10 times:

  - `u <- u0 + (dt/2)*(udot0 + udot1)`
  - recompute `udot1`

- Convergence criterion:
  - `||du|| / ||u|| <= tol`

- where:
  - `tol = min(1e-4, 0.1*accuracy)`
  - where `accuracy` is a user-defined parameter.

---

#### 6. RATTLE (velocity constraint projection)

The RATTLE velocity correction uses fixed-point iteration to enforce `G * M^-1 p = 0`. If the iteration fails to converge within the allowed iterations, the best available velocity estimate is accepted. Free-body quaternions are not advanced by the linear Taylor formula of step 1; they are advanced by an **exact exponential-map rotation** (step 2): the increment quaternion `advanceQuatExp(wHalf, dt) = (cos(theta), sin(theta) * wHalf/|wHalf|)` with `theta = (1/2)|wHalf| dt` and the midpoint angular velocity `wHalf = u0 + (dt/2) udot0`, applied as a left Hamilton product onto q0, followed by renormalization that mops up only rounding. This preserves unit norm and avoid dependence on `qddot`. The resulting constrained integrator remains second-order accurate globally, with third-order local error estimates that are computed but unused in the fixed-step sampler.

## Acceptance

The proposal is accepted or rejected by an MH test on the full Hamiltonian function H, using the **exact** dH between the pre- and post-trajectory states: `probability min(1, exp(-beta * dH))` where `dH = H_new - H_old`.

HMC (Duane, Kennedy, Pendleton & Roweth 1987) admits a **distinct guidance and acceptance Hamiltonian**. The *guidance* Hamiltonian generates the trajectory (here: `U` forces); the *acceptance* Hamiltonian is the one whose Boltzmann distribution is actually sampled, used in the MH test. Duane et al.'s key observation is that the trajectory generator need **not** equal the target: as long as the proposal map is reversible and volume-preserving and the MH test uses the *exact* acceptance Hamiltonian, the move samples the acceptance distribution exactly, whatever the guidance was. Two consequences used here:

- The Fixman torque is unnecessary for correctness because `U_Fixman` enters only the acceptance and the integrator need not compute `grad(U_Fixman)`
- An efficiency-oriented guidance potential is admissible: a cheaper or smoother potential (e.g. an internal-coordinate force field) may be used to *guide* proposals while the exact atomistic `H` is retained for acceptance, with no bias.

The holonomic constrained proposal map must be time-reversible and volume-preserving with respect to the constrained measure (or carry its Jacobian in the ratio); with exact `dH` of the acceptance Hamiltonian, the move is then exactly `pi`-stationary. Symplecticity is a strengthening of volume preservation that additionally yields a conserved shadow Hamiltonian and hence bounded energy error (an efficiency property, not a correctness one). For the separable Cartesian case, exact-arithmetic Verlet is reversible, volume-preserving, and symplectic. For the non-separable internal-coordinate case the truncated corrector and non-converged RATTLE acceptance can break reversibility and volume preservation themselves (the correctness properties) unless the truncation is shown to be a time-symmetric, volume-preserving involution. Until then, π-stationarity of the internal-coordinate move is an assumption, not a theorem.

A larger `dt` is a coarser Verlet map. Provided the map is reversible and volume-preserving and the **exact dH** is used in the MH test, the move is exactly `pi`-stationary at any `dt`. Assigning small `dt` to stiff blocks and large `dt` to soft blocks therefore introduces **no bias**, only a change in efficiency. This is the sense in which Robosample is a multiscale HMC sampler.

Reversibility is configuration dependent. It is **not** the curvature of S^3 that varies since that is **constant**, so the orientation manifold's intrinsic geometry is identical everywhere. What varies is the **mass metric M(q)** and the **local force stiffness d^2(U) / d(q)^2** (a near-clash or compressed region is far stiffer than an open one). A `dt` that is perfectly symplectic in an open conformation can become non-symplectic and begin pumping energy when the chain visits a stiffer geometry.

Any method that samples in a reduced (e.g. torsional) space while holding bonds/angles fixed must confront one fact: **the equilibrium distribution of the soft coordinates in a model with rigid bonds/angles is not the same as their marginal in the fully flexible model**. It has two distinct pieces, frequently conflated:

1. **A measure (kinetic / metric) piece.** Marginalizing the Gaussian momenta of a constrained system leaves a configuration-dependent factor `det(M(q))^(1/2)`. Uncorrected, the sampler targets `exp(-beta*U)*det(M)^(1/2)` instead of `exp(-beta*U)`. This is what the **Fixman potential** removes.
2. **A potential-of-mean-force (PMF) piece.** In a flexible molecule, bond angles *relax* in response to the torsional configuration (an angle widens to relieve a 1-4 clash at a given torsion). Freezing the angle removes that relaxation. This is an energetic effect: the location of the bond/angle minimum moves with torsion -- and **no determinant correction captures it**. It is recovered by **mobilizing the coupled hard DOF** through the
**fully flexible Cartesian world**, which alone covers every DOF and so suffices for correctness. Robosample relaxes them in **separate Gibbs worlds** corrected by an exact MH test, so no independence assumption is needed.

Running torsion-angle dynamics with rigid covalent geometry **grossly distorts** the torsional energy surface of Cartesian-parameterized force fields, but that can be repaired by building a specialized **internal-coordinate force field (ICFF)** (modified torsion (CMAP-style) terms plus softened van der Waals and
electrostatics) that approximates the source Cartesian field *without* a compensating potential. Robosample needs no such correction, for two independent reasons:

- **The proposal is corrected by an exact acceptance test.** Moreover, in Robosample, the rigid-geometry dynamics is only a **proposal**; the MH test uses the full atomistic `U` (plus Fixman), so the sampled distribution is exact regardless of how distorted the proposal forces are (the distinct guidance/acceptance design). At worst a poor proposal lowers acceptance; it cannot bias the result.
- **The frozen-geometry PMF distortion is integrated out.** The residual piece-2 distortion from holding bonds/angles fixed is removed not by patching `U` but by **mobilizing those coordinates** in the `Cartesian` worlds of the scan.
