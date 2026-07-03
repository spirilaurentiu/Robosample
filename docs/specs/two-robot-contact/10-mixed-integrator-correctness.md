# Two-robot contact world — Mixed-integrator correctness (CORE)

Terms from `00-...` §2. This file is the normative core: what must hold for a map
that advances R in generalized coordinates and E in Cartesian, sharing one
timestep and one force evaluation, to sample the correct joint Boltzmann density.

## 1. The target density

The mixed move SHALL sample the joint flexible-model Boltzmann density over
`(φ, x_E)`:

```
π(φ, x_E) ∝ exp(−β U(φ, x_E))                                            (T)
```

i.e. the target robot's marginal is the *flexible* (Gō–Scheraga eq 1) density with
NO residual mass-metric weight, matching the Cartesian environment's native
flexible density. `U` is the single shared OpenMM potential over all atoms (R's
atoms are functions of `φ`; E's atoms are `x_E`). This is the mixed-manifold
generalization of Spiridon & Minh 2017 eq:2→eq:3: HMC in `φ` naturally produces
the *rigid*-model weight `|M_φ|^{1/2} e^{−βU}` (eq:2), and `U_F` reweights it to
the flexible `e^{−βU}`.

## 2. Block-diagonality of the mixed mass metric (the load-bearing lemma)

**CLAIM M1.** With R and E atom-disjoint, the joint mass metric in `(φ, x_E)` is
block-diagonal:

```
M_mixed = Jᵀ M J = [ M_φ(φ)   0  ]      M_φ = J_Rᵀ M_R J_R,   M_E = const.
                   [   0     M_E ]
```

*Derivation.* Let `Q = (x_R, x_E)` be the full Cartesian coordinates, `M =
diag(M_R, M_E)` the diagonal atomic mass matrix (Spiridon eq:mass-metric-tensor,
`M` position-independent). The generalized coordinates are `(φ, x_E)`. The
Jacobian `J = ∂Q/∂(φ, x_E)` is block structured because R's atoms depend only on
`φ` and E's atoms are themselves the coordinates:
`∂x_R/∂φ = J_R`, `∂x_R/∂x_E = 0`, `∂x_E/∂φ = 0`, `∂x_E/∂x_E = I`. Hence
`M_mixed = Jᵀ M J` has zero off-diagonal blocks and `M_E` is the (constant)
diagonal Cartesian mass matrix of E's atoms. □

**PRECONDITION P1 (construction rule + ENFORCED guard).** R's atoms and E's atoms
are disjoint and no generalized coordinate of R moves any E atom. E's atoms are
flagged `cartSolventAtoms` AND — because the engine assigns *every* atom to some
articulated body (`buildModel`, `src/World.cpp:391-504`), so "E has no body" is NOT
realizable — E SHALL be constructed as one or more **0-DOF `Weld`-rooted, all-`Rigid`
bodies** whose atom membership is exactly the `cartSolventAtoms` mask. A 0-DOF body
contributes 0 to `calcLogDetM` (skipped by the `dof==0` guard) and 0 to the
articulated `ke` (`V_GB = 0`), while its atoms move per-atom in flat Cartesian via
`verletStep` (positions owned by `posG`, skipped by `fillAtomPositionsFromBodies`) —
so E is fully flexible in Cartesian yet invisible to the articulated metric, which is
what makes M1/M2/M3 hold. **If any flagged atom instead sits in a body with DOF > 0**
(e.g. a `Free`-rooted second robot — nothing in the engine forbids it), its KE is
**double-counted** (articulated `ke` + `keSolvent`), its momentum **double-drawn**
(`multiplyBySqrtMInv` + `drawSolventVelocities`), and `ΔU_F` picks up a spurious
E-dependent term: a silent, non-crashing Boltzmann corruption. This SHALL be enforced,
not assumed — `setCartesianSolvent` (`src/World.cpp:1559`) SHALL `throw` (fail-loud) if
any flagged atom's body has `bodyNU > 0`, or if a flagged body is only partially covered
by the mask. The zero cross-Jacobians in M1 depend on this. Bonded R–E contact (a shared
covalent link) violates P1 and is out of scope.

## 3. Fixman on the mixed manifold (question 3, decisive sub-answer)

**CLAIM M2 — the Cartesian subset contributes NO configuration-dependent Fixman
term; the existing robot-only Fixman is the complete and correct correction.**

*Derivation.* The mixed-manifold Fixman is Spiridon eq:3 with `M_mixed` in place
of the constrained metric and the full Cartesian `M_{3N} = diag(M_R^{cart}, M_E)`
as reference:

```
U_F^mixed = ½RT ln( |M_mixed| / |M_{3N}| )
          = ½RT ln( |M_φ|·|M_E| / (|M_R^{cart}|·|M_E|) )      (M1: det block-diag)
          = ½RT ln( |M_φ| / |M_R^{cart}| ).                    (|M_E| cancels)
```

The environment block `|M_E|` appears identically in numerator and denominator and
cancels **exactly**, because E is flat Cartesian (its constrained and reference
metrics are the same constant `M_E`). What remains is precisely the single-robot
Fixman `½RT ln(|M_φ|/|M_R^{cart}|)`. There is NO cross term. □

**M2a — reconciliation with the code as written.** `calcFixman` (`World.cpp:880`)
returns `½RT·(lnDetM − lnDetZ − lnDetMCartesian_)` where `lnDetM = calcLogDetM`
sums `ln det(D_b)` over the articulated bodies but **E's welded 0-DOF bodies are
skipped by the `dof==0` guard, contributing 0** (per P1), so `lnDetM = ln|M_φ(R)|`,
while
`lnDetMCartesian_` (`World.cpp:680-685`) sums `3·ln m_a` over **all atoms**,
i.e. `ln|M_R^{cart}| + ln|M_E|`. Therefore the code computes

```
U_F^code = ½RT( ln|M_φ| − ln|M_E| − ln|M_R^{cart}| ) = U_F^mixed − ½RT ln|M_E|.
```

The discrepancy is exactly `−½RT ln|M_E|`, a **constant** (E's atom set and masses
are fixed across a trajectory). It cancels in every acceptance ratio
(`ΔU_F = U_F(φ_end) − U_F(φ_start)`), so it does not bias sampling. The code is
therefore correct for a *static* environment atom set with no change required.

**INVARIANT INV-FIX (must-check).** For any fixed `E`, `U_F` computed with vs
without E's atoms in the Cartesian reference differ by a constant independent of
`(φ, x_E)`; consequently `ΔU_F` over any trajectory is bitwise identical whether or
not E is in `lnDetMCartesian_`. A test that perturbs `(φ, x_E)` and checks
`ΔU_F` invariance under adding/removing E from the reference SHALL pass. NOTE: if
a future design lets E's membership change *within* a run, the constant changes and
this invariant guards against silently importing that drift.

## 4. The momentum draw (question 3)

**CLAIM M3.** Drawing the two momentum blocks INDEPENDENTLY samples the correct
joint kinetic Gaussian `N(0, RT·M_mixed⁻¹)`.

*Derivation.* By M1 `M_mixed⁻¹ = diag(M_φ⁻¹, M_E⁻¹)` — block-diagonal, zero
cross-covariance. Hence `u ~ N(0, RT·M_φ⁻¹)` (Spiridon eq:momenta-draw, realized
by `multiplyBySqrtMInv`, `World.cpp:1686`) and `v_E ~ N(0, RT·M_E⁻¹)` drawn
per-component as `σ_a = √(RT/m_a)` (`drawSolventVelocities`, `World.cpp:1580`) are
independent draws of the correct joint law. The kinetic energy in `H`,
`ke + keSolvent = ½uᵀM_φu + ½Σ_E m|v|² = ½pᵀM_mixed⁻¹p`, matches. □

The same block-diagonality (M1) that kills the Fixman cross-term (M2) makes the
independent draw exact (M3) — one fact, two consequences.

## 5. Detailed balance / measure of the composed map (question 3)

**CLAIM M4.** The mixed move samples (T) exactly provided:
1. **(momentum refresh)** `(u, v_E)` are redrawn from `N(0, RT·M_mixed⁻¹)` per M3
   before each trajectory (full-refresh policy) — a Gibbs update of the momentum
   marginal that leaves `π(φ,x_E,u,v_E) ∝ e^{−βH'}` invariant;
2. **(proposal map)** the `L`-step composed integrator `T_L` on
   `(φ, u, x_E, v_E)` is deterministic, volume-preserving, and F-reversible
   (an involution up to momentum sign);
3. **(acceptance)** accept on the modified (acceptance) Hamiltonian
   `H' = H + U_F` — i.e. `min(1, e^{−β(H'_end − H'_start)})`.

*Why acceptance-only Fixman is exact here.* This is the CDHMC guidance/acceptance
split (Spiridon eq:modified-hamiltonian; Duane et al. 1987
`references/papers/duane_1987_hmc` guidance-vs-acceptance separation; Brubaker
et al. 2012 `references/papers/brubaker_2012_chmc` dual Hamiltonians). The
proposal is generated by the *unmodified* dynamics (no Fixman torque — R's ABA
dynamics already carry the position-dependent `M_φ` and the velocity-dependent
Coriolis/gyroscopic `f_p`, Spiridon eq:vv-momentum). A deterministic,
volume-preserving, F-reversible proposal accepted on `H'` leaves `e^{−βH'}`
invariant regardless of which Hamiltonian *generated* it. The `(φ,x_E)`-marginal
of `e^{−βH'}` is `|M_φ|^{1/2}|M_E|^{1/2} e^{−β(U+U_F)} ∝ e^{−βU}` (M2), which is
(T). □

**M4-risk (shared force symmetry).** Volume-preservation and F-reversibility of
`T_L` require ONE shared force evaluation per substep at the *joint* configuration
(`bridge_.evaluate` on the combined state — already the case, `verletStep`
`evalDerivs`), and the standard symmetric Verlet structure on the union
(position drift → shared force → trapezoid velocity for both blocks). The existing
`verletStep` satisfies this: E's position is frozen through the robot corrector, so
its trapezoid velocity `v1 = v0 + (h/2)(a0+a1)` (`RobotIntegrator.hpp:423-429`) uses
the post-corrector `frcG` — the *same* joint-configuration `a1` the robot corrector
converged to.

**INVARIANT INV-REV (must-check).** `checkReversibility` (`RobotIntegrator.hpp`,
called at `World.cpp:1443`) SHALL be evaluated on the JOINT (R+E) map at the
two-robot working `dt`, not on R alone. Round-trip residual `≤ 1e-6` (relative).
A residual that grows only when E is present localizes a coupling asymmetry in the
shared-force handoff. This is the runtime certificate for M4 condition 2.

**PRECONDITION P2 (runtime guard).** Under full refresh, no momentum flip is
needed (the redraw erases direction). Under partial refreshment (`20-...` §4) the
map SHALL negate BOTH `u` and `v_E` on rejection; a flip of one block only breaks
M4 condition 2 for the mixed metric. The corrector-nonconvergence path
(`RobotIntegrator.hpp:400-407`, "take the step anyway") is F-reversible-uncertain;
under full-refresh HMC its heat is Metropolis-rejected (self-correcting), but under
partial refreshment it SHALL be a reject/reduce-`dt` condition (inherits F3 from
`docs/specs/ncmc-explicit-solvent/20-inner-integrator.md`).

## 6. Touch list
- `World::setCartesianSolvent` (`:1559`): ADD the P1 construction guard — `throw`
  (fail-loud) if any flagged atom's body has `bodyNU > 0`, or if a flagged body is
  only partially covered by the mask. This is the enforcement P1 promises; today
  there is NONE, so a `Free`-rooted E corrupts sampling silently. This is the one
  required *new* code path in the correctness core.
- `World::calcFixman` (`:851`) / `lnDetMCartesian_` (`:680`): no change to the
  Fixman *arithmetic* (M2a — correct for a welded-E construction); add INV-FIX as a
  test. Convention at risk: `calcLogDetM` (skips E's 0-DOF bodies via `dof==0`) vs
  `lnDetMCartesian_` (all atoms) — document, do not "fix" by narrowing the reference
  (removes only a canceling constant and risks desync with `keSolvent`'s atom set),
  and do NOT "fix" `calcLogDetM` to explicitly exclude E — the `dof==0` skip already
  yields `ln|M_φ(R)|`.
- `World::reinitialize` (`:1608`) / `drawSolventVelocities` (`:1580`): the
  independent block draw is already correct (M3); no change.
- `RobotEngine::verletStep` (`:125`): shared-force symmetry already satisfied
  (M4-risk); the quaternion exp-map root advance (`:224-258`) SHALL stay
  reversible — do NOT "correct" it toward linear-Taylor (that reintroduces the KE
  pump, `30-...`).
- Conventions at risk: Frame F/M for the root joint (exp-map `w_FM` in F,
  `:241-252`); angular-over-linear ordering in the 6-DOF Free block; `M_φ` vs
  `M_R^{cart}` scope in Fixman.
