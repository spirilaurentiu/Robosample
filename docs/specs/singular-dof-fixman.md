# Spec: Singular / phantom generalized-coordinate DOFs and the Fixman potential

Status: IMPLEMENTED (Steps 2-3) + REVIEWED (merge). Step 4 (weld) DEFERRED. See Implementation outcome.
Owner: coder

## Implementation outcome (2026-07)

Implemented against the Review outcome (below), LEMMA-first, in `src/RobotEngine.cpp`:
- **§8 LEMMA ran first** — the per-phantom `ln D_b` swing across two configs measured **0 / 1.42e-14 / 0** (machine-eps), empirically REFUTING C4: it is a run-constant, no acceptance bias. Kept as `Cyclic1APQPhantomLogDetIsRunConstant`.
- **Step 2 (primary):** `logDetSymPD`→`pseudoLogDet` sharing one relative null-lock (`nullLockTol`/`eigDecompAndTol`) with `invertDense` (CC1/CC4/N1); a null direction contributes 0. `symSqrtInv` null direction → 0 (S3). `solveSmallSpd`/`symSqrt` floors left alone (S2).
- **Step 3 (fail-loud gate):** `realizeArticulatedBodyInertias` throws if a locked `D_b` is on a non-leaf-Torsion body. A locked *leaf Torsion* is a proven structural phantom (`D_b = n̂ᵀ I_M n̂`, a run-constant) and is allowed — this is the provably-safe substitute for the deferred build-time detector (avoids the B2 hazard).
- **Step 4 (build-time weld): DEFERRED** — not required for correctness (dynamics already lock the phantom unbiased; `pseudoLogDet` fixes `ln|M|`; the gate catches every non-leaf-Torsion lock). Only a wasted (deterministic, unbiased) RNG draw and `nq/nu` tidiness are lost.
- Reviewed: merge-ready on the science; full statistical gate green (25/25 slow tests), no tolerance loosened.

DEFERRED / follow-up (documented, not urgent):
- The build-time weld (Step 4).
- The gate currently HARD-THROWS (vs the §6 SHOULD-Metropolis-reject) on a *transient* config-dependent singularity of a REAL non-Torsion angle joint (`BendStretch`/`SphericalCoords`/…) mid-trajectory — unreachable today (all flexibilities are Torsion), but must become a reject-and-log once non-Torsion flexibilities ship.
- The config-dependent angle-joint singularity (KE↔Fixman divergence cancellation) — a separate spec (§6 NOTE).

---

Status (pre-implementation): REVIEWED — NOT ready for implementation; re-scope required (see Review outcome)
Owner (pre-impl): coder (after re-scope + the §8 LEMMA is run)

## Review outcome (hostile pass)

The theory core (C1, C3) and the fix *direction* (C6: weld structural phantoms / make `calcLogDetM` share the `invertDense` lock) are CONFIRMED sound. But the spec's headline justification is over-claimed and two `SHALL`s are unsafe as written:

- **C4 is REFUTED for the 1APQ case (was the central motivation).** For a leaf single-atom on-axis Torsion, `P_b` is exactly config-independent and `q_phantom` is frozen for the whole run (`u_phantom=0` always ⇒ `qdot=0`), so `ln D_phantom` is a **run-constant (~1e-16)** that **cancels exactly in `ΔU_F`** — there is NO O(RT) Metropolis-acceptance bias. The finding is therefore **demoted**: it is (a) a differential-oracle inconsistency (`calcLogDetM` returns a nonsense finite `≈−76` that can't be compared cross-engine) and (b) a *latent* robustness fragility (the residue only becomes config-dependent if a phantom's `q` can move — a null eigendirection of a multi-dof joint, or a config-dependent `P_b`). Re-ground the rationale on correctness-of-oracle + robustness, NOT acceptance bias. **The §8 "spurious Fixman term" LEMMA MUST be run first** to confirm the swing is ~machine-eps (predicted) before any code.
- **B2 — the §6 detector `SHALL` can weld a REAL config-dependent DOF.** Bullets (ii) collinear-subtree and the reference-geometry "within tol" test would weld an angle-flexible joint (`BendStretch`/`SphericalCoords`/…) that is merely collinear at the reference (θ→0/π pole) — the exact structural-vs-config-dependent conflation §1 forbids. Restrict the weld to (i) leaf single-atom on-axis **Torsion** and (iii) exactly-zero-mass subtree; gate on joint type; require a structural certificate that `D_b=0` for ALL `q`, never `D_b(reference)<lock` alone.
- **S1** CC4 is internally inconsistent (the nm "within tol" is a second threshold vs the amu·nm² inertia lock). **S2** do NOT force the hinge lock onto `Constraints::solveSmallSpd` (dimensionally distinct: `G M⁻¹ Gᵀ` is length²/mass). **S3** the `symSqrtInv` null convention is underspecified — it cannot be both pseudo-inverse (0) on a null dir AND preserve the exact-inverse identity `multiplyBySqrtM` (NMA route B) relies on; state which wins. **N1** the lock is RELATIVE (`max(1e-12, scale·1e-12)`), not absolute `1e-12`; the gate/`calcLogDetM` must reuse the relative logic.
- Confirmed safe: the fail-loud gate does NOT false-positive on light-but-real joints (H-torsion `D_b~1e-2 ≫ 1e-12`).

Net: the underlying inconsistency is real and worth fixing on oracle+robustness grounds, but it is **not an active sampling-correctness bug** for the tested systems. Re-scope C4, fix B2/S1/S2/S3/N1, and run the LEMMA before implementation.

---

Status (original): ready for review
Owner: coder (after review)
Scope: how the flexible internal-coordinate (Featherstone ABA) engine SHALL handle a generalized DOF whose hinge inertia `D_b = ~H_b P_b H_b` is null, and the consequences for `calcLogDetM` / `calcFixman`, the `sqrt(M^-1)` momentum draw, and Metropolis correctness. Read-only research artifact.

Triggering finding: the `1APQ` cyclic flexible model has a singular tree mass matrix (`minEigD ≈ 1.019e-33`, `tests/fixtures/robotics_oracle_molecules/1APQ_cyclic.moldyn.npz`) because it contains leaf single-atom `Torsion` bodies whose one atom lies ON the hinge axis, so `D_b = ~H_b P_b H_b = 0`.

NOTE (the reframing): the actual defect is NOT that the `1e-300` clamp fires (it never engages at `1e-33`). It is that two independent paths treat the singular direction inconsistently — `invertDense`'s `1e-12` null-space lock REMOVES it from the dynamics (unbiased), while the un-locked `logDetSymPD` KEEPS it in `ln det M` (adds a config-dependent residue). That inconsistency injects non-cancelling floating-point noise into the Fixman acceptance term.

## 1. Problem restatement

A *phantom DOF* is a generalized speed `u_b` (a column `H_b` of the hinge map) whose hinge inertia `D_b = ~H_b P_b H_b` (RobotEngine.cpp:761-768, 1085-1092) is null: the mobilizer motion moves no mass in its outboard subtree. In Fixman/mass-metric terms the internal-coordinate mass matrix `M_phi = J^T M J` is rank-deficient, `det M_phi = 0`, because the coordinate is a **gauge fiber**: the map `phi -> Cartesian` is constant along `u_b`, its Jacobian column is zero, and `phi -> q` is non-injective. Equivalently, one factor of the O(n) determinant identity `det M = prod_b det D_b` (Jain et al. 2013 eq:16) is zero.

The engine reaches this quantity by two INDEPENDENT and currently INCONSISTENT paths:

- **Dynamics path** — `invertDense` (RobotEngine.cpp:317-353) applies a **1e-12 null-space lock**: for `D_b ≈ 1.019e-33 < 1e-12` it returns `DI = 0` (a pseudo-inverse, not `1/D_b`). Downstream: `G_b = P_b H_b · DI = 0`, drawn `u_b = 0` in `multiplyBySqrtMInv` (`symSqrt(DI=0)=0`), `a_b = 0`, and `P^+_parent = P - G ~PH` unchanged. The phantom is *frozen and decoupled* — de-facto welded at runtime.
- **Determinant path** — `calcLogDetM` -> `logDetSymPD` (RobotEngine.cpp:73-105, 1071-1095) uses a **1e-300 Cholesky clamp** that does NOT engage at `1e-33`, so it ADDS `ln(1.019e-33) ≈ -75.9` to `ln|M_tree|`, treating the direction as present with a tiny eigenvalue.

Same singular direction: *removed* by the dynamics, *kept* by the determinant. The sampling-visible symptom is spurious non-cancelling floating-point noise in `Delta U_F` (§3, Claim C4).

**Reading chosen.** The reported case is a **structural** phantom: `D_b` is null for ALL `q` (single terminal atom on the rotation axis, fixed bond geometry). This is materially different from a **configuration-dependent** near-singularity (a real DOF whose `D_b(q) -> 0` only at isolated `q`, e.g. an angle-flexible joint at a collinear pole where `det J_B ∝ prod_i sin theta_i` vanishes; Jain 2013 eq:6). This spec handles the structural case; the config-dependent case is a separate failure mode (§6/§8, out of scope, flagged). Conflating them is the canonical correctness bug here and is explicitly avoided.

## 2. Binding definitions (theory + codebase)

- **Target marginal** (Spiridon & Minh 2017 eq:2): `rho(phi) ∝ |M_phi|^{1/2} e^{-beta U}`.
- **Fixman compensating potential** (Spiridon eq:3): `U_F(phi_f) = (1/2) beta^{-1} ln(|M_{N_f}| / |M_{3N}|)`, making the sampled marginal `∝ |M_{3N}|^{1/2} e^{-beta U}`; since `|M_{3N}|` is torsion-independent (Go-Scheraga; Jain eq:7), the flexible-torsion marginal becomes `∝ e^{-beta U}`.
- **Modified Hamiltonian** (Spiridon): `H' = H + U_F`, used ONLY in accept/reject; guidance dynamics use unmodified `H`.
- **Momentum draw** (Spiridon): `u ~ N(0, kT·M^{-1})`, realized `u = M^{-1/2} g` via `multiplyBySqrtMInv`.
- **O(n) determinant identity** (Jain eq:16-17): `det M = prod_b det D_b`, so `ln|M_tree| = sum_b ln det D_b` — exactly `calcLogDetM`.
- **Codebase Fixman assembly** (World.cpp:842-864): `U_F = (1/2) RT (ln|M_tree| - ln det(G M^-1 G^T) - ln|M_{3N}|)`, with `ln|M_{3N}| = 3 sum_{a: m_a>0} ln m_a` (World.cpp:664-667) — already EXCLUDES massless atoms.
- **Hinge inertia** (binding): `D_b = ~H_b P_b H_b`, `P_b` the articulated-body inertia. "Null" == `invertDense` treats a direction as null: `|eig(D_b)| <= max(1e-12, scale·1e-12)`.

## 3. Claims

- **C1 (geometry).** `D_b = 0` (structural phantom) iff the mobilizer motion `H_b` moves no outboard mass, for all `q`. Triggers: (i) a leaf single-atom `Torsion` body whose atom lies on the rotation axis (the `1APQ` case); (ii) a leaf/subtree all of whose atoms are collinear with the hinge axis; (iii) a massless/virtual site carrying a flexible DOF (`P_b` rank-deficient from zero mass). A phantom DOF moves no atom, changes no atom position, hence no `U`; it is a redundant gauge coordinate carrying zero physical information.
- **C2 (density).** With a structural phantom, `|M_{N_f}| ≡ 0` for all `phi`, so `rho(phi) ∝ |M_{N_f}|^{1/2} e^{-beta U} ≡ 0` is not normalizable and `ln|M_{N_f}| = -inf`. The theory-correct object is the pushforward to Cartesian, which omits the gauge coordinate. Correct behavior: **remove** the coordinate (weld to 0-dof); then `|M_{N_f}| > 0` and the marginal is well-defined. Removing it changes no observable (C1). Pruning is bias-free.
- **C3 (dynamics already unbiased).** The `1e-12` lock makes the dynamics use `M^+`: `DI_phantom = 0` => `u_phantom = 0`, `a_phantom = 0`, `G_phantom = 0` (no ancestor leakage), and KE via `V_GB` is unaffected (an on-axis point mass contributes `1/2 m|v_atom|^2 = 0`). The draw/accelerations/forces/KE are NOT biased by the reported phantom.
- **C4 (the determinant path DOES bias acceptance).** `calcLogDetM` does not share the lock: it adds `ln(D_phantom)` with `D_phantom ≈ 1.019e-33` a floating-point RESIDUE whose value depends on the current frame `R_GB` and thus DIFFERS between trajectory endpoints. Hence `Delta U_F` acquires `(1/2) RT [ln D_phantom(new) - ln D_phantom(old)]`, an O(RT) term that does NOT cancel, entering `min(1, exp(-beta Delta H'))` every step and distorting acceptance of the PHYSICAL DOFs. Not a clean rejection, not a controlled approximation.
- **C5 (the clamp comment is wrong on three counts).** (i) The `1e-300` clamp does not engage at `1e-33`, so it protects nothing here; (ii) a hugely negative `ln|M|` drives `U_F` hugely NEGATIVE, which FAVORS the singular direction, not rejects it; (iii) the operative mechanisms are the `invertDense` lock (freezes dynamics, good) and the un-locked log-det (injects C4 noise, bad), which disagree. "Reject cleanly" is not the behavior.
- **C6 (fix = consistency + removal).** Make both paths treat the gauge direction identically and match `|M_{3N}|`'s existing massless-atom exclusion: remove the phantom from BOTH the state and the determinant (weld at build), or equivalently compute the pseudo-determinant `prod_{non-null} det D_b` using the SAME null-space lock as `invertDense`. Build-time weld preferred (structural fact, no runtime threshold coupling, no frozen coordinate, no wasted RNG draw); the pseudo-determinant is the safety net.

## 4. Derivation sketch (codebase notation)

Leaf single-atom body, atom at body origin `Bo` on torsion axis `n̂` through `Bo`.
- Point-mass inertia (`c=0`): `P_b = [[0,0],[0,mI_3]]` in `[angular;linear]` blocks.
- Torsion hinge in Ground: `H_b = [omega; v]`, `omega = R n̂`, linear part `v = ` velocity of `Bo` under unit rotation about `n̂` through `Bo` `= 0`. So `H_b = [R n̂; 0]`.
- `P_b H_b = [0; 0] = 0` => `D_b = ~H_b (P_b H_b) = 0` for all `q`. Structural. (Off-axis atom => `v != 0`, `D_b > 0`.)

Consequences: `det M_tree = prod_b det D_b` has a zero factor => `ln|M_tree| = -inf` analytically, numerically `logDetSymPD` returns `≈ -75.9`. `invertDense(D_phantom)`: `DI = 0` => `G_b = 0`, `P^+_parent` unchanged. Draw: `u_phantom = 0`. KE: contributes 0. Acceptance: `Delta U_F` inherits the non-cancelling O(RT) residue difference (C4).

Config-dependent contrast (out of scope, do NOT conflate): angle-flexible `det J_B = sin theta_ex · d_2^2 prod_{i>=3} d_i^2 sin theta_i` (Jain eq:6) makes `det M_B ∝ prod sin^2 theta_i` vanish smoothly as `theta_i -> 0/pi`. There `rho ∝ |M_{3N}|^{1/2} e^{-beta U}` stays FINITE (Fixman cancels the suppression), so the treatment is the KE<->Fixman divergence cancellation, NOT welding.

## 5. Correctness conditions

- **CC1.** `calcLogDetM` and `invertDense` MUST agree on which directions of every `D_b` are null. A direction removed by the dynamics pseudo-inverse MUST contribute factor 1 (`0` to `ln det`), never `ln(residue)`.
- **CC2.** After the build pass, no flexible body is a structural phantom: `min-eig(D_b) > lock threshold` for every body with `bodyNU > 0` at the reference geometry.
- **CC3.** Welding a structural phantom changes no atom position, no `U`, no physical KE at any `q` (bias-free).
- **CC4.** The null-direction convention (the `invertDense` lock threshold `max(1e-12, scale·1e-12)`, MD units amu·nm^2) is the SINGLE source of truth shared by `invertDense`, `calcLogDetM`, and the runtime gate. No second independent threshold (the `1e-300` clamp, the `symSqrt`/`symSqrtInv` floors, `solveSmallSpd`) may define "singular".
- **CC5.** `|M_{3N}|` already excludes massless atoms; the tree determinant MUST be consistent (a massless-site DOF MUST NOT contribute a null/`-inf`/residue factor).

## 6. Recommended engine behavior (SHALL / SHOULD / NOTE)

- **SHALL — build-time structural-phantom detection and weld (primary fix).** During `World::build`/`buildModel`, every flexible body whose hinge inertia is null for all `q` SHALL be pruned to 0-dof (`JointType::Rigid`, `bodyNU=0`). Detection SHALL be structural: classify as phantom iff the outboard subtree carries no mass the hinge moves — (i) a leaf single-atom body with its atom on the hinge axis within `tol`, (ii) a leaf/subtree with all outboard atoms on the hinge line within `tol`, or (iii) an outboard subtree of total mass 0. WHY: removes the gauge coordinate from BOTH paths consistently (CC1, C6); config-independent; matches the massless-atom exclusion in `|M_{3N}|` (CC5); bias-free (CC3).
- **SHALL — make `calcLogDetM` consistent with the dynamics inverse.** `logDetSymPD`/`calcLogDetM` SHALL apply the SAME null-space lock as `invertDense` (CC1, CC4): a null direction contributes 0 to `ln det`, giving `sum_{non-null} ln det D_b`. Removes the spurious `ln(residue)` (C4) for anything that escapes build-time welding, and makes the differential oracle `logDetM` reference-comparable.
- **SHALL — replace the silent clamps/floors with the shared lock.** The `1e-300` clamp in `logDetSymPD` (:90-92), the floors in `symSqrt`/`symSqrtInv` (:426,436,449,459), and in `Constraints::solveSmallSpd` (Constraints.cpp:216-218) SHALL NOT define "singular" independently — they never engage at the physical scale (`1e-33`) and only convert underflow into a meaningless finite value that hides the error (C5). Subordinate them to the shared lock (CC4).
- **SHALL — runtime fail-loud gate.** At `realizeArticulatedBodyInertias`, if `invertDense` null-locks any direction of a `D_b` on a body NOT marked phantom-welded, fail loud (throw with body id, `JointType`, atom indices, `min-eig(D_b)`), not silently lock-and-continue. Reuses the `invertDense` threshold, so it fires exactly when a direction is locked (CC4) and does NOT false-positive on a light-but-real torsion (`D_b ~ O(1e-2) >> 1e-12`).
- **SHOULD — distinguish transient from structural.** A transient config-dependent near-singularity mid-trajectory SHOULD force a Metropolis REJECTION + log (unbiased: measure-zero in the target); a build-time structural phantom SHOULD be welded (never reached at runtime); a structural phantom surviving build detection SHOULD hard-throw.
- **NOTE.** Guidance dynamics use unmodified `H`; no Fixman torque path to correct — only the acceptance `ln|M|` and the draw/KE. Draw/KE already unbiased (C3); only `ln|M|` needs the fix.
- **NOTE (out of scope).** Config-dependent singularities of angle-flexible joints (`BendStretch`, `SphericalCoords`, `Cartesian`, `FreeLine`), where `det J_B ∝ prod_i sin theta_i` vanishes at a collinear pole (Jain eq:6/7), are NOT structural phantoms and MUST NOT be welded; they need the KE<->Fixman divergence-cancellation analysis (separate spec, §8).

## 7. Touch list

- `src/RobotEngine.cpp`: `logDetSymPD` (:73-105) + `calcLogDetM` (:1071-1095) remove `1e-300` clamp, apply the `invertDense` lock => pseudo-determinant; `invertDense` (:317-353) expose per-body null-count so `calcLogDetM`+gate share it (CC4); `symSqrt`/`symSqrtInv` (:424-464) subordinate floors to the lock; `realizeArticulatedBodyInertias` D/DI (:761-768) add the fail-loud gate; `calcKineticEnergy` (:1101-1110) no change (document `V_GB` path robustness, C3).
- `src/World.cpp`: `buildModel` (near :655-670) add the build-time phantom detector + weld; reuse the massless-atom exclusion at :664-667 as the consistency target (CC5); `calcFixman` (:842-864) no formula change (correct once `ln|M_tree|` is a pseudo-determinant / phantom welded).
- `src/Constraints.cpp`: `calcConstraintLogDet`/`solveSmallSpd` (:139-171, :216-218) align floor with the shared lock (singular loop-closure `G M^-1 G^T` is an analogous separate case).
- `include/RobotModel.hpp`: weld sets `bodyJoint=Rigid`, `bodyNU=0`; detector consumes `bodyRootAtom`, `bodyAtoms*`, `atomMass`, hinge geometry.

Conventions at risk: F/M vs Ground for the hinge (on-axis condition = Ground-frame linear part of `H_b` at `Bo` vanishes, RobotEngine.cpp:516-523); the `1e-12` lock as single source of truth across `invertDense`/`calcLogDetM`/gate; `det M = prod_b det D_b` becoming a pseudo-determinant; the massless-atom exclusion already in `|M_{3N}|`.

## 8. Verification plan

- **PRECONDITION (runtime guard).** At `realizeArticulatedBodyInertias`, for every body `bodyNU>0` not phantom-welded: `invertDense` locks NO direction of `D_b`; else throw with `(body, JointType, atoms, min-eig)`. (The §6 fail-loud gate.)
- **INVARIANT — build removes structural phantoms (CC2).** Build `1APQ_cyclic`; assert every body the oracle flags at `minEigD ≈ 1e-33` has `bodyNU==0` in the built `RobotModel`, and `calcLogDetM` is finite with all per-body `D_b` above the lock. Assert on the BUILT model, not the input.
- **INVARIANT — pruning is bias-free (CC3).** For each detected phantom, set its would-be coordinate to distinct values; assert the atom Cartesian position, `U`, and per-body physical KE are identical to machine precision across them (that invariance IS the definition). If any value moves the atom / changes `U`, it was NOT a phantom and the test must fail.
- **INVARIANT — path consistency (CC1, C6).** On an UNPRUNED model with a phantom (detection disabled), assert `calcLogDetM` (locked) equals `sum_{non-null} ln det D_b` and excludes `ln(residue)`; and that drawn `u`, `A_GB`, KE, ancestor `P` are bitwise unchanged whether the phantom is welded or lock-pruned (they must agree — C3).
- **LEMMA — the spurious Fixman term (C4).** On the UNPRUNED, UN-locked-`logDetM` path, evaluate `ln|M_tree|` at two configs of `1APQ_cyclic`; the phantom contributes `ln D_phantom(config)` (config-dependent FP residue) => a non-cancelling `O(2-5)` spurious swing not attributable to any physical DOF. After the fix: swing 0. NOTE: the magnitude of this swing is asserted from reasoning, not yet measured — this LEMMA measures it (confirm/bound rather than take on faith).
- **LEMMA — log-det gap on the reference.** `ln|M_tree|_fixed - ln|M_tree|_buggy ≈ -ln(D_phantom) ≈ +75.9` per phantom body.
- **INVARIANT — gate fires correctly, no false-positive.** Throws on an unpruned phantom; does NOT fire on the pruned `1APQ_cyclic` nor on acyclic references whose smallest physical `D_b` is `O(1e-2)` (ala-dipeptide, butane). Tolerance relative to body inertia scale, never absolute `1e-300`.
- **INVARIANT — Boltzmann/detailed-balance on the pruned model.** The existing `FixmanBoltzmann`/Shirts gate MUST pass, exercising the draw (`multiplyBySqrtMInv`), KE (`calcKineticEnergy`), AND Fixman `ln|M|` (`calcLogDetM`) TOGETHER (the defect is a per-term inconsistency).

## 9. Open questions

- **Q1.** Are massless/virtual sites ever assigned a flexible joint in the target force fields? If so, the mass-based detector branch (C1.iii) is load-bearing (CC5).
- **Q2.** Drop the pruned coordinate from `q`/`u` (smaller state) or keep and hard-lock? Preferred: drop (SHALL, §6); if IO/analysis assumes fixed `nq`/`nu`, a lock-in-place variant may be needed (does not change correctness).

## References

- Spiridon & Minh 2017 (`references/papers/spiridon_2017_cdhmc_gibbs/equations.md`): eq:2, eq:3, eq:momenta-draw, eq:modified-hamiltonian, eq:6.
- Jain et al. 2013 (`references/papers/jain_2013_fixman_branched/equations.md`): eq:6, eq:7, eq:16, eq:17.
- Fixman 1974 (`references/index.yaml` key `fixman_1974`).
- Codebase: `src/RobotEngine.cpp` (logDetSymPD :73-105, invertDense :317-353, symSqrt/symSqrtInv :424-464, realizeABI D/DI :761-768, multiplyBySqrtMInv :992-1018, calcLogDetM :1071-1095, calcKineticEnergy :1101-1110); `src/World.cpp` (lnDetMCartesian_ :664-667, calcFixman :842-864); `src/Constraints.cpp` (:139-171, :216-218); `include/RobotModel.hpp` (:27-177); `docs/specs/robotics-oracle-differential.md` (§2b null-space lock, row `logDetM`).
