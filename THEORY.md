# Robosample -- Theory of Operation

Robosample performs **Gibbs sampling coupled with Hamiltonian Monte Carlo (HMC)**. A
configuration is updated by a fixed-order sweep of *Gibbs blocks*; each block proposes a new
configuration by integrating constrained molecular dynamics and accepts or rejects it with a
Metropolis-Hastings (MH) test on a Hamiltonian. Blocks run either in **internal coordinates**
(bond-angle-torsion; currently only torsions are mobile) or in **Cartesian coordinates**
(ordinary MD). The target is the canonical (constant-temperature) Boltzmann distribution;
energies are evaluated on the full atomistic model with the OpenMM force field and GBSA-OBC2
implicit solvent.

Notation convention: "~" means "distributed as"; "proportional to" is written out; superscripts
use "^" (e.g. M^-1, p^T, sin^2), and "1/2" denotes one half.

---

## 1. Notation

| Symbol | Meaning |
|---|---|
| q | generalized coordinates **mobile in the current block** (a subset of torsions, plus the 6 external DOF of any free root body) |
| p | generalized momenta conjugate to the mobile q |
| M(q) | configuration-dependent generalized mass-metric (articulated mass matrix) of the block's reduced multibody tree |
| D_b | per-body articulated factor (the dof x dof block produced by the articulated-body algorithm at body b); ln det M(q) = sum_b ln det(D_b), computed in O(n) without forming M |
| M_3N | constant 3N x 3N Cartesian mass matrix (diagonal atomic masses); the reference in the Fixman term and the metric used in Cartesian blocks |
| V(q) | potential energy (OpenMM force field + GBSA-OBC2), evaluated on the full atomistic model |
| K(p,q) | kinetic energy, (1/2) p^T M(q)^-1 p |
| F(q) | Fixman compensating potential (Section 6) |
| J(q) | external-rotation Jacobian for free root bodies (Section 6) |
| T, R | temperature; molar gas constant R = 8.3144626e-3 kJ/(mol*K). Energies are per mole (kJ/mol); beta = 1/(R*T) |
| sigma(q), d0 | ring-closure constraint sigma(q) = |r_AB|^2 - d0^2 = 0, with target distance d0 (the carried-over closure distance) |
| G(q) | constraint Jacobian d(sigma)/dq of the holonomic ring-closure constraints sigma(q) = 0 |
| gamma2_b | polar ("pitch") Euler angle of free root body b, extracted from its orientation quaternion |
| pi | target distribution; Z is the canonical partition function |

---

## 2. Target distribution

The target is the canonical Boltzmann distribution

```
pi(q) = Z^-1 * exp( -beta * V(q) )
```

(constant N, T; with implicit solvent there is no simulation box, so there is no volume or PV
term -- "NVT" here means constant-temperature canonical sampling). Internal-coordinate sampling
is designed so that the marginal distribution of configurations matches the marginal of the
**fully-flexible Cartesian** Boltzmann distribution restricted to the sampled (torsional)
subspace; the Fixman term (Section 6) removes the mass-metric artifact that would otherwise
distort this marginal.

---

## 3. Internal-coordinate representation (robot model)

Each molecule is represented as a robot (kinematic tree), and the Featherstone articulated-body
algorithm propagates positions, velocities, forces, and accelerations up and down the tree in
O(n). Per-atom Cartesian forces (-grad V) computed by OpenMM are reduced to per-body spatial
forces (net force + torque about the body origin) for the articulated-body solver.

**Freezing by rigid-body lumping (no zero-velocity holds, no Schur complement).** A block that
freezes a coordinate does **not** hold that DOF at zero velocity. Instead, rigidly connected
atoms are lumped into composite rigid bodies (welds), so the frozen coordinates simply do not
exist in that block's reduced tree. Consequently each block has a **single reduced articulated
mass matrix M(q)**, and the same M(q) drives both momentum sampling and the Fixman determinant.
Because (i) momenta are fully resampled at the start of every block -- nothing is carried across
blocks to condition on -- and (ii) each block's M(q) is the genuine articulated mass of its own
reduced tree rather than a sub-block of a larger matrix, **no Schur complement arises**.

---

## 4. The Gibbs scan

Blocks execute in a fixed (deterministic) order; an MH test is applied after each block. Each
block's HMC move is reversible with respect to pi by construction and therefore leaves pi
invariant. The composite scan is generally **not** reversible -- its time-reversal is the
reversed block order -- but it still leaves pi invariant, because composing kernels that each
leave pi invariant preserves pi. **Stationarity (pi*K = pi), not reversibility, is the
requirement.**

- **Overlap.** Torsion subsets may overlap between blocks. Redundant updates do not break
  stationarity, provided each block individually leaves pi invariant.
- **Ergodicity (precise statement).** Across one full round (all blocks), every mobile coordinate
  in the block set has been *offered* a move. Ergodicity is over the **sampled (torsional)
  subspace** and requires both that the blocks collectively cover that subspace and that each
  block achieves nonzero acceptance. Being offered a move is **not** the same as accepting one: a
  block stuck near zero acceptance contributes nothing even though aggregate statistics look
  healthy. Bond lengths and angles are not mobile in torsion-only operation, so they are not
  sampled; the sampled distribution is the (Fixman-corrected) torsional marginal.
- **Per-block diagnostics.** Acceptance rates are reported **per block** specifically to detect an
  under-sampled subspace, which would otherwise be masked by healthy global energy and acceptance
  statistics.

---

## 5. The HMC move inside a block

1. **Momentum resampling.** At the start of each block the momenta are **fully** resampled (not
   partial): p ~ N(0, M(q)), equivalently velocities u ~ N(0, R*T*M(q)^-1). This is a Gibbs step
   on the momenta (always accepted) and follows the internal-coordinate equipartition principle.
   The draw is matrix-free via the articulated sqrt(M(q)^-1) operator, so det M(q) never appears
   explicitly in the proposal; its configuration dependence is bookkept by the Fixman term in the
   acceptance ratio (Section 6).
2. **Propagation.** Dynamics are integrated with **velocity Verlet** (internal and Cartesian
   alike) under the **bare potential V only**. Constraints are then projected onto the manifold
   (Section 7).
3. **Acceptance.** The proposal is accepted or rejected by MH on the full Hamiltonian H
   (Section 6).

**Guidance vs. acceptance Hamiltonian.** The Fixman potential and the external-rotation Jacobian
(their gradients / "torques") do **not** enter the forces -- the leapfrog is guided by V alone.
They enter **only** the MH acceptance, through H. This is the Duane-et-al. distinct
guidance/acceptance design; the Fixman torque is not required for correctness.

**Time step is a property of the proposal, not the target (multiscale HMC).** A larger dt is a
coarser leapfrog map. Provided the map is reversible and volume-preserving and the **exact dH** is
used in the MH test, the move is exactly pi-stationary at any dt. Assigning small dt to stiff
blocks and large dt to soft blocks therefore introduces **no bias** -- it changes only efficiency.
This is the sense in which Robosample is a multiscale HMC sampler.

---

## 6. The acceptance Hamiltonian

```
H(q, p) = V(q) + (1/2) * p^T * M(q)^-1 * p + F(q) + J(q)
```

- **V(q)** -- potential energy (force field + GBSA-OBC2), full atomistic model.
- **K = (1/2) p^T M(q)^-1 p** -- kinetic energy; M(q) is the reduced articulated mass-metric of
  the block (Section 3).
- **F(q) -- Fixman compensating potential.**

  ```
  F(q) = (1/2) * R * T * [ ln det( M(q) )  -  ln det( G(q) * M(q)^-1 * G(q)^T ) ]
  ```

  The first term cancels the det(M(q))^(1/2) factor that emerges when the Gaussian momenta are
  marginalized, so that the configurational marginal matches the target Boltzmann distribution
  without the mass-metric artifact. The second term is present **only when holonomic ring-closure
  constraints are active**: the RATTLE momentum projection onto G*M^-1*p = 0 injects a
  det(G*M^-1*G^T)^(-1/2) factor into the marginal, which this term cancels. For acyclic molecules
  the second term vanishes (no constraints), and F reduces to the standard tree determinant
  (ln det M = sum_b ln det(D_b)). Additive constants (e.g. the constant Cartesian reference
  ln det M_3N) cancel in dH.
- **J(q) -- external-rotation Jacobian.**

  ```
  J(q) = -(1/2) * R * T * sum_{free root bodies b} ln sin^2( gamma2_b )
  ```

  This is the volume-element Jacobian of the orientation parameterization of any molecule whose
  root body is **free** (6 external DOF: 3 translation + 3 rotation), where gamma2_b is the polar
  ("pitch") Euler angle of body b taken from its orientation quaternion. It restores the uniform
  (Haar) measure on the body's orientation. The translational external DOF contribute a unit
  Jacobian. J(q) is **identically zero when all roots are welded** (e.g. single-molecule torsional
  sampling in a fixed frame); it is nonzero only for free/floating bodies (e.g. a docked ligand).
  The bond-length (r^2) and bond-angle (sin theta) factors of the full bond-angle-torsion volume
  element are **constant** under torsional dynamics -- those coordinates are constrained -- so they
  cancel in dH and do not appear in H.

**Cartesian special case.** In Cartesian coordinates M = M_3N is the constant diagonal
atomic-mass matrix, so F and J are configuration-independent constants (taken as zero) and the
move reduces to standard Cartesian HMC/MD.

The MH acceptance uses min(1, exp(-beta * dH)) with the **exact** dH between the pre- and
post-trajectory states, including F and J.

---

## 7. Constraints

Holonomic constraints are applied **only in internal coordinates** and **only on ring-closing
(cotree) bonds**, which cannot be represented as tree coordinates. Each is enforced at both levels
of a RATTLE scheme:

- **Position (SHAKE):** after the Verlet position update, q is projected so that
  sigma(q) = |r_AB|^2 - d0^2 = 0.
- **Velocity (RATTLE):** after the velocity update, the momenta are projected so that
  G(q)*M(q)^-1*p = 0 (no relative velocity along the closed bond).

Both projections are required; either alone leaves a secular drift. RATTLE preserves the
phase-space measure **on the constraint manifold**, so the volume-preservation premise of the HMC
move (Section 5) holds with constraints active. The constraint contributes the det(G*M^-1*G^T)
term to F (Section 6).

The ring-closure distance d0 is part of the configuration carried between blocks: it tracks the
value left by the previous block (continuation), and is **not** reset to an idealized force-field
bond length.

---

## 8. Forces, energy, and solvent

- In internal coordinates, per-atom Cartesian forces -grad V are computed by OpenMM and summed into
  per-body spatial forces for the articulated-body solver.
- Energies are evaluated on the **full atomistic model** (not coarse-grained).
- Implicit solvent is **GBSA-OBC (Onufriev-Bashford-Case, model II / OBC2)**, OpenMM's standard
  implementation.

---

## 9. Inter-block continuation

State is handed from one block to the next as the full configuration (Cartesian positions, which
encode every bond length, angle, and torsion). On entering a block, that geometry becomes the
block's reference: rigid-body shapes are rebuilt from the **actual** incoming coordinates, so any
coordinate the block freezes is held at its **current carried-over value**, not an idealized one.
The handed-over state is the previous block's post-MH configuration -- the new point on
acceptance, or the retained point on rejection. Momenta are not carried over; they are freshly
resampled in each block (Section 5).

---

## 10. Worked example: dH for one torsional block (acyclic, welded root)

For a single internal-coordinate block with a welded root (so J = 0) and no ring closures (so the
G*M^-1*G^T term of F is absent), one move proceeds as follows.

1. The block receives configuration q_old (Section 9) and evaluates V_old = V(q_old).
2. Momenta are resampled: p_old ~ N(0, M(q_old)). Then

   ```
   K_old = (1/2) * p_old^T * M(q_old)^-1 * p_old
   F_old = (1/2) * R * T * ln det M(q_old)
   H_old = V_old + K_old + F_old
   ```

3. Velocity Verlet integrates the block's mobile torsions under V only (no Fixman torque),
   producing (q_new, p_new).
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

- det M(q) enters dH only through F, never through the momentum draw (Section 5).
- The bond-length, bond-angle, and (trivial) torsion Jacobian factors are constant and cancel,
  so they do not appear in dH.
- With a ring closure present, add the constraint contribution
  -(1/2)*R*T*[ ln det(G*M^-1*G^T)|_new - ln det(G*M^-1*G^T)|_old ] to dH.
- With a free root body, add J_new - J_old (Section 6) to dH.
