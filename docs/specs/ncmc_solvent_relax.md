# SPEC: Solvent-relaxing NCMC — co-integrated mixed internal/Cartesian HMC

## Claim

Add a sampling mode in which, during an NCMC-style solute move, the **contact
environment (explicit solvent, or a co-mobilized shell of a contacting
molecule) relaxes inside the proposal** rather than being welded to an
infinite-mass wall. The mode samples the exact joint canonical distribution

```
pi(q, p_q, x_s, p_s) ∝ exp(-beta · H_joint)
H_joint = U(q, x_s) + ½ p_qᵀ M(q)⁻¹ p_q + ½ Σ_s |p_s|²/m_s
          + U_Fixman(q) − ½ R T · logSineSqr(q)
```

where `q` = solute internal (torsional/Free-root) coordinates with mass-metric
`M(q)`, `p_q` its generalized momenta, and `(x_s, p_s)` = **Cartesian**
positions/momenta of solvent (or shell) atoms that have LEFT the robot tree and
are integrated by OpenMM. `U(q,x_s)` is the single full-system OpenMM potential
(at λ=1). The new term versus today is the solvent Cartesian kinetic energy
`½ Σ_s |p_s|²/m_s`; `U_Fixman` and `logSineSqr` are **unchanged** (they depend on
`q` only).

This is the joint mixed-coordinate HMC the user described ("torsional dynamics
for solute + Cartesian OpenMM integration for solvent, in one Gibbs block").
Alchemical softening becomes an **optional** layer on top, not a prerequisite —
which matters because for GPCR-in-membrane / receptor–ligand contact, alchemy
does little and the real fix is letting the partner move inside the move.

## Architectural facts confirmed by reading the code

- `verletStep` (`RobotIntegrator.hpp`) reaches OpenMM only through
  `bridge.evaluate(s)` = `setAtomPositionsInGround(ALL atoms)` +
  `getForcesFromOpenMM`. Every OpenMM force eval **overwrites all OpenMM
  positions from the robot `atomPosG` array** — solvent positions are slaved to
  the robot tree.
- `calcLogDetM` (`RobotEngine.cpp:1066-1071`) **skips zero-DOF (Weld) bodies**
  (`det=1 → log 0`), so welded solvent contributes nothing to `ln|M_tree|`.
- `lnDetMCartesian_` (`World.cpp:664-668`) = `3·Σ_a ln m_a` over **all atoms incl.
  solvent** — a config-independent constant.
- `calcKineticEnergy` sums `½ uᵀMu` over robot bodies only; welded solvent
  contributes 0 KE.
- The exact-`dH` acceptance argument (`World.cpp:1934-1951`) rests entirely on the
  protocol being a deterministic, volume-preserving, F-reversible map `T` with
  `F·T·F = T⁻¹`, closed by the λ-palindrome pinned at λ=1.

## The two failure modes, kept separate (Rule 3)

1. **Measure/metric piece** (`det M(q)^{1/2}`, removed by Fixman). Solvent is flat
   Cartesian: constant diagonal metric ⇒ **zero** Fixman contribution and
   **zero** `logSineSqr`. Removing solvent from the tree does not perturb the
   solute `det M(q)` (proof §4). Already correct; stays correct.
2. **PMF piece** (frozen hard-DOF relaxation). Exactly what is broken today: the
   welded solvent cage never relaxes, so the λ=1 endpoint pays the full unrelaxed
   insertion energy and self-rejects. The fix **mobilizes** those frozen Cartesian
   DOF inside the move so the environment accommodates the moved solute. This is
   the entire point; do not conflate with the Fixman piece.

---

## Path comparison

### Path 1 — deterministic joint NVE Verlet, keep the exact-`dH` acceptance (RECOMMENDED)

Co-integrate solute internal coordinates and solvent Cartesian coordinates as ONE
deterministic, volume-preserving, momentum-flip-reversible map over the joint
phase space, sharing one OpenMM force evaluation per step. Solvent moves under
**plain NVE Verlet (no thermostat)** with SETTLE/rigid constraints. Accept on the
joint `dH` at λ=1 — the existing machinery, extended only by the solvent KE term.

**Why this is the correct minimal increment:** it reuses the proven exact-`dH`
argument verbatim. NVE Verlet conserves the joint `H_joint` to O(dt²); the
solute's insertion energy flows into solvent KE instead of spiking `U`, so `dH`
stays small and acceptance recovers. With λ≡1 (no alchemy) it is a clean joint
HMC that solves the docking/membrane case directly.

**Correctness conditions (Path 1):**

- **Per-substep F-reversibility.** Each joint substep operator `V_s` must satisfy
  `F·V_s·F = V_s⁻¹`. Build each as a **symmetric (Strang) palindrome of
  F-reversible factors**:
  ```
  V_s = S½ ∘ Q ∘ S½
  ```
  where `Q` = the existing solute internal Verlet step (`RobotEngine::stepTo`,
  already F-reversible — certified by `checkReversibility`), and `S½` = a half
  solvent Verlet op (half-kick on `p_s` + half position drift), each F-reversible.
  A palindromic composition of F-reversible factors is F-reversible:
  `F(S½ Q S½)F = S½⁻¹ Q⁻¹ S½⁻¹ = (S½ Q S½)⁻¹`. **Do NOT use a Lie–Trotter split**
  (`A_solvent_full ∘ B_solute_full`): `F(AB)F = A⁻¹B⁻¹ ≠ (AB)⁻¹`, silently breaks
  the exact-`dH` test. Hazard #1.
- **Outer λ-palindrome unchanged.** With `V_s` F-reversible and `λ_s = λ_{N-1-s}`
  pinned at 1 (current `protocolLambda`), the whole protocol `T = V_{N-1}…V_0`
  satisfies `F·T·F = T⁻¹`. Proof closes identically to today's.
- **Volume preservation.** `Q` is Verlet (vol-preserving); solvent kicks/drifts are
  unit-Jacobian shears. Joint map is volume-preserving.
- **No thermostat inside the move.** A stochastic thermostat on solvent destroys
  determinism and the deterministic-map `dH` test; thermostatting ⇒ Path 2.

### Path 2 — work-based NCMC (Nilmeier/Crooks/Minh & Chodera 2011)

At each fixed λ, run a `π_λ`-preserving propagation that includes native OpenMM MD
on solvent (SETTLE, full steps, optionally thermostatted), interleaved with solute
internal propagation; accept on accumulated generalized work.

**Correctness conditions (Path 2):**

- Each fixed-λ propagation kernel must **exactly preserve `π_λ`** (joint canonical
  measure at that λ). Plain Verlet does NOT (shadow drift); you must either
  (a) Metropolize/GHMC the propagation, or (b) account the **shadow work** (the
  `H_joint` change during deterministic propagation segments) in the accumulator.
  The current code's `work` (`World.cpp:1856-1868`) accumulates ONLY the
  fixed-config perturbation terms and is explicitly diagnostic-only — **not** a
  valid Path-2 accumulator as written.
- **Acceptance** (generalized work incl. kinetic terms, since solvent velocities
  are now state):
  ```
  P_accept = min(1, exp(-beta · w_total))
  w_total = Σ_perturb  [ U_{λ_{t+1}}(q_t, x_{s,t}) − U_{λ_t}(q_t, x_{s,t}) ]
          + Σ_propagate[ shadow-H drift if deterministic Verlet is used ]
  ```
  with `U_Fixman(q)` and `logSineSqr(q)` folded into `U_λ` consistently at every
  λ-perturbation evaluation (state functions of `q`; their drift across the
  protocol is part of `w_total`).
- Momentum handling per Nilmeier 2011 velocity-randomization: resample
  `(p_q, p_s)` from Maxwell–Boltzmann at the start of each move; within the
  protocol momenta are deterministic state and must be flipped on reverse if not
  resampled.

**Cost:** higher-variance acceptance; must implement correct shadow-work or
Metropolized propagation — substantial new machinery, departs from the exact-`dH`
philosophy. **Benefit:** native OpenMM MTS/Langevin/SETTLE robustness for large
boxes where NVE shadow drift over a long protocol erodes acceptance; the only
viable route for **dual-thermostat (Drude)** systems.

---

## 2. Exact acceptance expression (chosen path = Path 1)

At the λ=1 endpoints:

```
Hstart = pe0 + ke_q0 + ke_s0 + fixman0 − ½·RT·logSineSqr0
Hend   = pe1 + ke_q1 + ke_s1 + fixman1 − ½·RT·logSineSqr1
dH     = Hend − Hstart
accept iff  u ~ U(0,1) < exp(−beta·dH)    (min(1,·) form)
```

Relative to today's `currentTotalEnergy()`:
- `pe` = full OpenMM PE at λ=1 over **all** atoms (solute Cartesian image of `q` +
  live solvent `x_s`). Same call (`bridge_.calcPotentialEnergy()`), but `x_s` now
  reflects the relaxed solvent.
- `ke_q` = `RobotEngine::calcKineticEnergy` (solute generalized `½uᵀMu`). Unchanged.
- **`ke_s` = `½ Σ_s |p_s|²/m_s` (NEW)** = solvent Cartesian KE, from OpenMM
  `State::Velocities` of the solvent subset.
- `fixman`, `logSineSqr` — **unchanged**, function of `q` only.

**What cancels:** the constant solvent-mass term in `lnDetMCartesian_` (and any
constant frame offset in `logSineSqr`) is identical at start/end ⇒ cancels in
`dH`. **What does NOT cancel and IS the productive signal:** `ke_q + ke_s + pe`
— solvent KE absorbs the insertion energy the welded wall used to dump into `pe`.
**Gate:** with solvent frozen (no `p_s` drawn, `ke_s≡0`) and λ≡1, this reduces
EXACTLY to the current torsional-HMC Metropolis test (`World.cpp:1971`).

## 3. KE / momentum bookkeeping

- **Solvent velocity draw:** at move start (alongside `reinitialize`'s solute `u`
  seeding, `World.cpp:1543-1626`), draw `v_s ~ N(0, RT/m_s)` via
  `OpenMMContext::setVelocitiesToTemperature(T, seed)` **restricted to the solvent
  subset** (solute frozen in the solvent context, §5). Compute `ke_s0`, add to
  `Hstart`.
- **Where it enters:** add `state_.energy.ke_solvent` (NEW field) to both `Hstart`
  and `Hend`. Keep it separate from `state_.energy.ke` (solute generalized KE) so
  the existing NMA `ke_mix` correction (`World.cpp:1662-1664`) and Fixman
  bookkeeping are untouched.
- **Momentum-flip reversal:** the F-map reversal must flip BOTH `u` (solute) AND
  `p_s` (solvent). For per-move acceptance, the existing "resample-each-block ⇒ no
  explicit flip needed" argument (`World.cpp:1847-1853`) extends to `p_s` provided
  `p_s` is redrawn every move and never read across moves. For the
  `checkReversibility` probe on the JOINT map (§7), an explicit `p_s → −p_s` flip
  is required.
- **Interaction with `sqrt(M⁻¹)` seeding:** none — solute seeding
  (`multiplyBySqrtMInv`) is over `q` DOF only; solvent has a trivial diagonal
  metric ⇒ independent Gaussian draw.

## 4. Fixman & Jacobian implications (confirmation + the one flag)

- **Solute `det M(q)` unchanged by removing solvent from the tree.** `calcLogDetM`
  (`RobotEngine.cpp:1066-1071`) skips zero-DOF bodies; welded solvent already
  contributes 0 to `ln|M_tree|`. Deleting those bodies changes the sum by exactly 0.
- **No Fixman, no `logSineSqr` from solvent.** Flat Cartesian metric ⇒ constant
  determinant; no Free-root quaternion ⇒ no pitch `γ2` term.
- **The one det-M flag:** `lnDetMCartesian_` (`World.cpp:664-668`) sums `3·ln m_a`
  over ALL atoms incl. solvent — a config-independent constant ⇒ cancels in `dH`
  regardless. **Decision:** KEEP solvent atoms in `model_.atomMass`/`atomPosG`
  (needed for the OpenMM force/position contract and a constant Cartesian
  reference); remove them only from the **body/tree** topology. If an implementer
  drops solvent from `model_.atomMass` instead, the reference constant changes but
  still cancels — acceptable, but document it. Either way, assert `fixman` is
  **byte-identical** before/after on a solute-only frozen-solvent run (§7).

## 5. Minimal change surface

- **Insertion point: `World.cpp:1912` (step ii) is the ONLY one.** Replace the
  single `RobotEngine::stepTo(...)` with the symmetric joint substep
  `S½ ∘ Q ∘ S½`:
  1. `S½`: half-kick solvent velocities + half-drift solvent positions in OpenMM,
     using forces at the current config.
  2. `Q`: `RobotEngine::stepTo(...)` — solute internal step (unchanged). Its
     `bridge.evaluate` force calls must see the **current solvent positions**;
     holds automatically because solvent positions live in `atomPosG` and are
     written back after each `S½`.
  3. `S½`: trailing half-kick using forces at the post-drift config; read solvent
     velocities/positions back into the Cartesian buffers.
  After the trailing `S½`, write updated solvent `x_s` into `state_.atomPosG()`
  solvent slots so the next substep's `setAtomPositionsInGround` pushes the
  relaxed solvent (not a stale frozen copy).
- **OpenMM-side mechanism to advance ONLY solvent — use a dedicated frozen-solute
  context, NOT mass-toggling the live context.** Build once a second OpenMM
  `System`/`Context` whose **solute particle masses are 0** (OpenMM freezes
  zero-mass particles, `OpenMMContext.cpp:141`), solvent masses physical, SETTLE
  intact, a `VerletIntegrator`. Each substep: `setPositions(all)` (solute image
  from robot, solvent from buffer) → integrate solvent → read back solvent
  positions+velocities. Rationale: live-context mass changes force a reinitialize
  per substep (prohibitive); position restraints on solute are not cleanly
  reversible and bias the result. Frozen-mass Verlet on the moving subset is
  deterministic, volume-preserving, F-reversible (flip solvent velocities).
  - **Half-kick exposure (hazard #1):** OpenMM's monolithic `VerletIntegrator::
    step()` does full kick-drift-kick internally and does **not** expose
    half-kicks, so naive `step(1)` between solute steps yields a Lie–Trotter
    (non-symmetric, non-reversible) split. Implement `S½` via a `CustomIntegrator`
    with separately addressable half-kick / drift `addComputePerDof` ops (the
    codebase already builds a `CustomIntegrator` for MTS, `OpenMMContext.cpp:30`).
- **Solute → OpenMM sync between substeps:** `RobotEngine::
  fillAtomPositionsFromBodies` already writes solute Cartesian positions to
  `atomPosG`; push those to the solvent context each substep via `setPositions`.
- **Force reduction must skip solvent:** with solvent removed from the tree,
  `getForcesFromOpenMM` (`ForceBridge.hpp:79-98`) must NOT scatter solvent atom
  forces onto any robot body. Map solvent atoms to a sentinel body or skip by an
  `isSolventAtom` mask.

## 6. Reversibility hazard list (each with its test)

1. **Asymmetric (Lie–Trotter) splitting** → `F·T·F ≠ T⁻¹`. Test: §7 joint
   `checkReversibility` residual must stay ~machine-eps; Lie–Trotter gives
   O(dt·force) residual.
2. **OpenMM not exposing half-kicks** → de-facto asymmetric split. Test: same as
   #1; unit-test the `CustomIntegrator` `S½` against an analytic harmonic
   oscillator for self-reverse.
3. **Thermostat / velocity reseeding mid-protocol** (Path 2 hazard leaking into
   Path 1) → breaks determinism. Test: assert no `setVelocitiesToTemperature`
   between move start and acceptance.
4. **SETTLE/constraint reversibility:** SETTLE is analytic+reversible; a one-sided
   SHAKE tolerance is not. Test: joint reversibility probe with constrained water.
5. **`enforcePeriodicBox` wrapping** of coords fed back to the robot engine
   (`OpenMMContext` getState uses `enforcePeriodicBox`, e.g. `:344`). If solvent
   positions are box-wrapped on one leg but not the reverse, the map is not
   reversible. Test: reversibility probe across a periodic boundary. **Mitigation:**
   feed UNWRAPPED solvent coords to the robot-side force eval; wrap only for output.
6. **Solvent momentum not flipped** on reverse → probe fails. Test: §7 probe flips
   `p_s`.
7. **Solvent positions overwritten by `setAtomPositionsInGround`** — if the solvent
   buffer is not written back into `atomPosG` after `S½`, the solute force eval
   re-freezes solvent. Test: assert solvent `atomPosG` actually changes across a
   substep at λ=1.

## 7. Test plan (encodes intent)

- **T0 — exact reduction (must hold first).** λ≡1, solvent frozen (`ke_s≡0`, no
  solvent draw): assert `ncmcMove()` `dH`, accept decision, and `fixman` are
  **bit-identical** to the current welded path. Proves no regression in the
  measure piece.
- **T1 — joint reversibility probe.** Extend `RobotEngine::checkReversibility` (or
  add a joint variant) to flip `(u, p_s)` and integrate the JOINT `S½ Q S½` map
  forward/back. Residual ~machine-eps at safe dt, O(1) at too-large dt —
  distinguishes the correct symmetric split from Lie–Trotter (hazard #1/#2). Wire
  to `set_reversibility_check`.
- **T2 — single explicit water, analytic.** One TIP3P/SPC water + a rigid/one-
  torsion solute, small box. (a) λ≡1: verify `⟨ke_s⟩ → (3/2)·N_solvent_dof·RT`
  (equipartition) and that the solute marginal matches a long reference Cartesian
  HMC. (b) Two-particle harmonic surrogate with known joint Gaussian — assert
  sampled covariance matches analytic.
- **T3 — Fixman/Boltzmann invariance.** Reuse `TestFixmanBoltzmann.cpp` /
  `TestMassScaleInvariance.cpp`: solute configurational marginal unchanged by the
  presence of mobilized solvent.
- **T4 — acceptance trends (the actual goal).** On 2ala in a small explicit box:
  acceptance-vs-dt and acceptance-vs-`ncmc_steps` must now **improve** (curve
  shifts up, dt ceiling rises) vs the welded baseline. Falsifiable success
  criterion; a biased impl would show high acceptance but fail T0/T1/T5.
- **T5 — ensemble validation (bias gate).** Log-probability-ratio slope test
  (`references/index.yaml`, ensemble-validation ~240-251): two runs at slightly
  different β give log-ratio of energy histograms linear with slope `−(β2−β1)`.
  Catches a biased-but-accepting sampler that T4 alone would not.
- **T6 — physical barrier needing solvent relaxation.** A contact-gated move
  (solute reorientation against a packed solvent shell) rejected under welding must
  now be accepted, AND the resulting populations pass T5.

## 8. Recommendation

**Adopt Path 1 (deterministic joint NVE Verlet, exact-`dH`) as the primary:**
(i) reuses the proven exact-`dH` acceptance, adds only one energy term + a
symmetric splitting — smallest correct increment; (ii) at λ≡1 it is a general
joint mixed-coordinate HMC that fixes docking/membrane contact WITHOUT alchemy
(alchemical softening layers on top as optional efficiency); (iii) NVE
co-integration keeps joint `H` near-conserved, which is exactly why acceptance
recovers — low, deterministic variance.

**Keep Path 2 (work-based) as the documented fallback** for (a) large solvent
boxes where NVE shadow drift over a long protocol erodes acceptance and native
OpenMM Langevin/MTS+SETTLE is more robust, and (b) **Drude / dual-thermostat**
systems, out of scope for Path 1.

**Phased implementation order (smallest correct increment first):**
0. T0 reduction + T1 joint reversibility on a frozen-then-mobilized single water
   (no alchemy, λ≡1).
1. One explicit water mobilized, λ≡1 joint HMC; pass T2.
2. 2ala in a small explicit box, λ≡1; pass T3/T4/T5 — demonstrate the
   acceptance-vs-dt improvement.
3. Layer the alchemical λ:1→0→1 softening on top of the now-mobilized solvent.
4. Generalize "solvent subset" to a co-mobilized local shell of a contacting
   molecule (docking site, membrane lipids) — the GPCR/membrane case.

## Open questions the implementer must resolve before coding

1. **Half-kick exposure:** can a `CustomIntegrator` expose solvent half-kick +
   drift as separate F-reversible ops sharing the SAME force eval as the solute
   step, or must `verletStep` be refactored into explicit kick/drift
   sub-operations? (Determines whether the change is contained to `World.cpp:1912`
   or touches `RobotIntegrator.hpp`.)
2. **Periodic box / minimum image:** which coords (wrapped vs unwrapped) go to the
   robot-side `bridge.evaluate` for the solute–solvent force, and whether
   `enforcePeriodicBox` must be disabled on substep getState calls (hazard #5).
3. **SETTLE vs HMR-flexible water:** rigid (SETTLE, larger dt) or flexible (needs
   HMR)? Sets the solvent dt ceiling and whether the joint substep shares one dt
   with the solute step or needs an inner RESPA sub-cycle (codebase has r-RESPA
   `MTSIntegrator`, `OpenMMContext.cpp:30`).
4. **dt matching:** solvent NVE dt == solute internal dt, or RESPA-nest `k` solvent
   substeps per solute step? Nesting must itself be a symmetric palindrome.
5. **Frozen-solute context lifecycle:** one persistent second context with zeroed
   solute masses vs reusing the main context — memory/perf on 1M-atom systems, and
   whether virtual sites (TIP4P-Ew M-site, `ForceBridge.hpp:80-89`) are correctly
   recomputed each substep.
6. **Drude:** confirm deferred to Path 2 (dual thermostat); Path 1's single NVE
   stream cannot maintain the cold-Drude constraint.

## References

- NCMC framework — Nilmeier, Crooks, Minh & Chodera, *PNAS* 108(45):E1009 (2011),
  doi:10.1073/pnas.1106094108 (cited inline at `World.cpp:1833`).
- Fixman/coordinate-Jacobian — Spiridon & Minh, *JCTC* 13:4649 (2017), cited at
  `World.cpp:836` and `references/index.yaml`.
- Internal-coordinate EOM reversibility/conservation — Mazur, Dorofeev & Abagyan
  (1991), doi:10.1016/0021-9991(91)90210-C (`references/index.yaml` line 79).
- Ensemble-validation slope test — `references/index.yaml` ensemble-validation
  entry (~240-251).
