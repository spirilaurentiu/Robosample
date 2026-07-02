# Spec: Robotics oracle — mobilizer reaction forces (Scope-B pin + live-field check)

Status: ready for review
Owner: coder (after review)
Parent: `docs/specs/robotics-oracle-differential.md` (the campaign; sections cited below as "differential spec §N"), `docs/specs/robotics-oracle-data-provenance.md` (provenance conventions), `docs/specs/singular-dof-fixman.md` (the singular/phantom carve-out this reuses).

This spec closes one confirmed coverage gap: **mobilizer reaction forces on real molecules are unpinned by the oracle, and no fixture exercises a reaction under a live molecular force field.** It adds two independent legs. It does **not** re-derive the ABA recursion, the reaction operator, or the Scope-B build/correspondence machinery — those are binding from the parent specs and are cited, not repeated.

## 1. Problem

### 1.1 In the user's words
"Reaction forces on real molecules are never checked; the only reaction test that runs under a live force field is the thing we ship (`World::calcSpatialForces`), and nothing validates it."

### 1.2 In the codebase's vocabulary
A **mobilizer reaction force** is the spatial force `[angular=torque; linear=force]` that a body's inboard mobilizer transmits to it, expressed in Ground, reported either at the body-frame origin `Bo` (`reactionBo`) or at the outboard mobilizer-frame origin `Mo` (`reactionMo`). Both engines report the force on the **child** body; the force on the parent at `Mo` is its negative (`SimbodyMatterSubsystem.h:2393-2396`). Generalized (mobility) forces applied at the joint are **included** in the reaction, not subtracted (`SimbodyMatterSubsystemRep.cpp:6060-6061`; port note `RobotEngine.cpp:1243-1253`).

Three states of coverage today:

1. **Port self-consistency (already covered — do not duplicate).** `TestReactionForces.cpp` checks `RobotEngine::calcMobilizerReactionForces` / `findMobilizerReactionOnBodyAtMInGround` (`src/RobotEngine.cpp:1257-1315`) against its own Newton-Euler recursion, the free-mobilizer law, and a whole-tree balance, driven by the analytic harmonic bridge (no OpenMM). All port-internal.
2. **Scope-A synthetic cross-engine pin (already covered — do not duplicate).** The Scope-A oracle pins `reactionBo`/`reactionMo` vs live Simbody at `1e-8` (`tests/TestRoboticsOracle.cpp:429-441`, generator `Robosample/tools/gen_robotics_oracle.cpp:380-385`), on hand-built `MobilizedBody` systems, **including** applied-force discriminator states (differential spec §7.1 S3, §8.2 #7). So the reaction *math* and the `[torque;force]` slot / `Bo`↔`Mo` shift conventions are already externally anchored on synthetic bodies with nonzero applied force.
3. **THE GAP.** The Scope-B molecule dump (`Robosample/src/RoboticsOracleMoleculeDump.cpp`) dumps `q,u,udot,X_GB,V_GB,A_GB,minEigD,KE,logDetM` but **no reactions**; `tests/TestRoboticsOracleMolecule.cpp` asserts none. Reactions on a **real molecule's** articulated tree — where the per-atom→per-body force reduction, the branched inward accumulation, and the real mass dynamic range live — are never compared to Simbody. And `World::calcSpatialForces` (`Robosample/src/World.cpp:4355`), which stores exactly `reactionForces[mbx][0]` (torque), `[mbx][1]` (force), `u`, `uDot` per interesting body under the **live** OpenMM/DuMM field, is never validated against an independent method.

Two legs close the gap, each pinning a different quantity:

- **Leg A** (numerical differential, in the existing Scope-B fixture path): pin the **zero-applied-force** reaction on the three §4.4 molecules cross-engine, reusing the differential's staged comparison and its singular/phantom carve-outs.
- **Leg B** (live-force-field check, clone-side): validate the **live-field** reaction (exactly what `calcSpatialForces` stores) on a real molecule against Simbody's *independent* `calcMobilizerReactionForcesUsingFreebodyMethod` and a whole-tree Newton-Euler balance.

## 2. Binding conventions (from the parent specs and the two source functions)

- **Spatial layout** `[angular; linear]` on both sides; port force order `[torque; force]` (differential spec §7.2; `RobotEngine.cpp:1172`). `reactionBo`/`reactionMo` are `SpatialVec` = `[torque; force]` at `Bo`/`Mo`, expressed in Ground.
- **The reaction identity (Jain 2011, Eq. 7.34; `SimbodyMatterSubsystem.h:2421-2429`; `SimbodyMatterSubsystemRep.cpp:6045-6094`):** `F_reaction@Bo = zPlus[mbx] + PPlus[mbx]·APlus`, then shift `Bo→Mo`.
- **`shiftForceBy` sign (binding; `SpatialAlgebra.h:625-626`):** `shiftForceBy(F, r) = SpatialVec(F.angular − r % F.linear, F.linear)`. The linear (force) part is unchanged by the shift; only the torque changes by `−r × f`.
- **`findMobilizerReactionOnBodyAtMInGround(s)` == `calcMobilizerReactionForces(s)[mbx]`** (proven bit-identical from clone source; `MobilizedBody.cpp:406-414` vs `SimbodyMatterSubsystemRep.cpp:6076-6094`; confirmed at runtime `TestMobilizerReactionForces.cpp:166-169`).
- **`Bo` is a valid cross-engine anchor; `X_GB.R()` is not.** The differential spec §6 NOTE establishes that across two independently-built molecule models the body-frame *orientation* `X_GB.R()` legitimately differs by convention, while `X_GB.p()` (the body origin `Bo`) agrees (asserted at `kMolStage1Tol`). This governs §4 below: `reactionBo` is convention-free, `reactionMo` is not.

## 3. Claims

- **C1.** The Scope-B reference reaction at zero applied force can be reconstructed *inside the existing dump* from buffers `calcAccelerationForOracle` already returns (`zPlusAll`, `A_GB_All`) plus the `PPlus` getter, reproducing `calcMobilizerReactionForces`'s math on the zero-force `zPlus` — **without** calling the public reaction accessor (which would realize acceleration under the live DuMM/OpenMM field and contaminate the zero-force comparison).
- **C2.** `reactionBo` (torque+force at the body origin, expressed in Ground) is a **convention-free** cross-engine anchor; `reactionMo`'s **angular** part is convention-contaminated by `X_GB.R()` and is comparable only after the physical point `Mo` is shown to agree cross-engine. `reactionMo.linear == reactionBo.linear` is always comparable (the shift leaves the force unchanged).
- **C3.** At zero applied force the reaction is driven purely by the velocity-dependent (gyroscopic/Coriolis) term; the meaningful teeth of Leg A is the **`random` (u≠0)** state of the well-conditioned **regular** class (`10ala` flexible). The `rest` (u=0) state and the 0-DOF rigid class give identically-zero reactions on both sides (a weak but real guard); `1APQ`'s `random` state is singular and its reactions are skipped by the existing guard.
- **C4.** Leg A's cross-engine reaction diff inherits the differential's **exact** tolerance and guard structure (Stage-4 tier + FINDING#1 widening for the regular class + the whole-tree singular skip + the per-body phantom-angular skip). It introduces no new physics tolerance, only a magnitude-relative scaling appropriate to a force (§4.3).
- **C5.** Leg B validates the live-field reaction with **zero new Simbody patch**: `calcMobilizerReactionForcesUsingFreebodyMethod` is already public (`SimbodyMatterSubsystem.h:2964`) and assembles force-element + constraint + gyroscopic contributions independently of the PPlus method — a genuine cross-method check on the real molecule's live DuMM force reduction.
- **C6.** `World::calcSpatialForces` has a physically-correct guard with an inverted message (`World.cpp:4356-4358`); the fix is message-only.

## 4. Leg A — Scope-B zero-force reaction pin

### 4.1 Reconstruction (clone side, inside `dumpMoleculeOracleState`)

For each mobilized body `b` (mbx `b`, parent mbx `p`), reusing the buffers already produced by `matter.calcAccelerationForOracle(state, 0, 0, udotAll, A_GB_All, zAll, zPlusAll, epsilonAll)` (`RoboticsOracleMoleculeDump.cpp:196-201`; `zPlusAll` is currently captured-but-unused):

```
p_PB_G = X_GB.p() − X_GP.p()                                   // Ground, parent→body origin
A_GP   = A_GB_All[p]                                          // ZERO-force parent acceleration
APlus  = SpatialVec(A_GP[0], A_GP[1] + A_GP[0] % p_PB_G)      // parent accel shifted to Bo
FB_G   = zPlusAll[b] + PPlus[b] · APlus                       // reactionBo (Ground)
p_BM_G = X_GB.R() * mobod.getOutboardFrame(state).p()
FM_G   = shiftForceBy(FB_G, p_BM_G)                           // reactionMo (Ground)
```

with `PPlus[b] = matter.getArticulatedBodyInertiaPPlusForOracle(state, b)` and, when the parent is Ground, `APlus = 0` so `FB_G = zPlusAll[b]`. This is `SimbodyMatterSubsystemRep::calcMobilizerReactionForces` (`.cpp:6076-6094`) verbatim, with the two force-carrying inputs (`getArticulatedBodyForcesPlus`, `parent.getBodyAcceleration`) replaced by the zero-force `zPlusAll`/`A_GB_All`. `PPlus` is force-independent (a function of configuration and mass only, valid at `Stage::Position`+; `SimbodyMatterSubsystem.h:2841-2843`), so the persisted getter is consistent with the zero-force reconstruction.

- **NOTE (the load-bearing subtlety).** `A_GP` **SHALL** be `A_GB_All[p]` (the zero-force buffer), **not** `parent.getBodyAcceleration(state)` (the live-field cache). Using the live accessor here is the single most likely reconstruction bug and would silently agree with the port only when the force field is absent.
- **NOTE.** The public `findMobilizerReactionOnBodyAtMInGround` / `calcMobilizerReactionForces` accessors **SHALL NOT** be used to produce the Scope-B reference: they trigger `getArticulatedBodyForcesPlus(s)`, realizing acceleration under the registered DuMM/OpenMM forces — inconsistent with the §4.5 zero-force carve-out the whole molecule dump is built on.

### 4.2 Dump / npz / manifest additions

- The dump struct (`Robosample/include/RoboticsOracleMoleculeDump.hpp::BodyDump`) **SHALL** gain `double reactionBo_ang[3]`, `reactionBo_lin[3]`, `reactionMo_ang[3]`, `reactionMo_lin[3]` (Ground-expressed). The pybind marshaller (`Robosample/src/PyBind11.cpp:84-128` `dump_robotics_oracle_molecule`) **SHALL** emit them as `reaction_bo_ang`/`_lin`, `reaction_mo_ang`/`_lin`.
- The generator (`tests/fixtures/robotics_oracle_molecules/_generate_molecule_oracle.py`) **SHALL** write, per state/body, `state{s}_body{b}_reactionBo_ang`, `_reactionBo_lin`, `_reactionMo_ang`, `_reactionMo_lin` (each length-3 float64), following the existing per-body naming (`_generate_molecule_oracle.py:188-208`).
- The loader (`tests/RoboticsOracleMoleculeLoader.hpp::MoleculeOracleBody`) **SHALL** gain the four fields and read them in `loadMoleculeOracleCase`.
- `manifest.json` `schema_version` **SHALL** bump `1 → 2`; the generator writes 2 and the arrays list gains the four names. Old fixtures then fail loud (the loader's `need()` throws on a missing array), forcing regeneration.

### 4.3 Port-side comparison (in `runMoleculeNumericDifferential`, `tests/TestRoboticsOracleMolecule.cpp`)

The port already sets `bodyForceG = 0` and `mobilityForce = 0` and runs `realizePosition/Velocity/ArticulatedBodyInertias/calcUDot` (`.cpp:571-598`) — the exact zero-force state the reconstruction assumes. After Stage 4, add a reaction stage:

- Compute the port reference once per state: `RobotEngine::calcMobilizerReactionForces(m, s, reacBo.data(), reacMo.data())` (both buffers; `src/RobotEngine.cpp:1257`). This is the port's Newton-Euler freebody recursion at zero applied force — an independent method from the clone's PPlus reconstruction, so Leg A is cross-**engine** *and* cross-**method**.
- **Primary anchor (Normative): `reactionBo`.** Per matched body (atom-set map, §5 of the differential spec), compare `s`-side `reacBo[pb]` vs clone `state{s}_body{cb}_reactionBo_*` (angular and linear), Ground-expressed.
- **`reactionMo` (Recommended, gated).** `reactionMo.linear` equals `reactionBo.linear` and MAY be compared directly. `reactionMo.angular` **SHOULD NOT** be diffed cross-engine unless the physical point `Mo` is first shown to agree: assert `‖p_Mo_G(port) − p_Mo_G(clone)‖ ≤ kMolStage1Tol` where `p_Mo_G = X_GB.p + X_GB.R·X_BM.p`, and only then compare `reactionMo.angular`. If `Mo` disagrees, report it (it is itself a finding about the mobilizer-frame convention) and skip the `reactionMo.angular` diff. The reconstruction stores `reactionMo` regardless, for this gate and for diagnostics.

### 4.4 Guards and tolerances (reuse, do not re-invent)

- **Whole-tree singular skip.** The reaction diff **SHALL** sit under the existing `st.minEigD <= kMolMinEigDEps` guard (`.cpp:681-687`): `udot`/`A_GB`/reactions all depend on the same (near-)singular `D`-solve, so a singular tree makes the cross-engine reaction comparison ill-posed (this is what skips `1APQ`'s `random` reactions).
- **Per-body phantom-angular skip.** For a body with `bd.minEigD <= kMolMinEigDEps`, skip the **angular** reaction comparison (it inherits the gauge-unobservable hinge-axis direction via `A_GB.angular`), keep the **linear** comparison. Same structure as Stage 2/4 (`.cpp:651-666`, `710-724`).
- **Tolerance.** Reactions are Stage-4/5-class end products. They **SHALL** use the same `stage4Tol` parameter already threaded into `runMoleculeNumericDifferential` (default `kMolStage4Tol = 1e-8`; `kMolFlexibleStage4Tol = 1e-3` passed only by `RegularFlexibleNumeric`, per the FINDING#1 mass-property widening `.cpp:317-339`) — the regular class's clone-side mass-property reconstruction inaccuracy propagates into `reaction ≈ Mk·A_GB` and accumulates through the inward `Σ Φ·reac_child`, exactly as it does for `A_GB`/`udot`.
- **Force scaling (Recommended).** A reaction is a force, so compare with a magnitude-relative bound `|a−b| ≤ stage4Tol · max(1, ‖ref‖)` (as `logDetM`/`KE`/`‖udot‖` already do, `.cpp:729-741`), **not** the plain absolute `NearVec3` used for `A_GB` — reaction magnitudes scale with body mass (`10ala` masses up to ~16 Da), so an absolute `A_GB` tolerance would be miscalibrated. Calibrate the concrete constant empirically with ≥2× margin over the largest legitimate residual, as the file's other tiers were (`.cpp:285-356`), keeping it ≥3 orders below an O(1)-relative transcription bug.
- **NOTE.** Per-case teeth: `10ala_regular` `random` (well-conditioned, u≠0) is the only state that exercises a nonzero reaction non-trivially; `10ala_rigid` and every `rest` state give zero reactions both sides; `1APQ_cyclic` `random` is skipped (singular). State this in the test so a reviewer does not mistake the weak cases for strong ones.

### 4.5 §4.5 comparability (cyclic / `1APQ`)

No special-casing beyond what the differential spec §4.5 already guarantees at the dump/engine-call level: the clone reference comes from `calcAccelerationForOracle` (never the dead Rod-constraint path) and the port reference comes from `calcMobilizerReactionForces` run after `calcUDot` (pre-projection; SHAKE/RATTLE runs strictly later, `RobotIntegrator.hpp`). Both are the **unconstrained** tree reaction. The reaction is computed on the spanning tree only; the shared-tree precondition (`assertRingClosingSetsMatch`, `.cpp:489-501`) already gates the `1APQ` state comparison.

### 4.6 Regeneration

- Rebuild+install the clone bindings, then run (differential spec §11 step 1, generator header): `python3 tests/fixtures/robotics_oracle_molecules/_generate_molecule_oracle.py`.
- Fixtures regenerated: `10ala_rigid`, `10ala_regular`, `1APQ_cyclic` (`.moldyn.npz` + `.moldyn.manifest.json`). The `*.systopo.npz` (port topology) fixtures are unchanged.

## 5. Leg B — live-force-field reaction check (clone side)

### 5.1 What it validates and why it is not Leg A
Leg B validates the quantity `World::calcSpatialForces` actually stores: the reaction under the **live** DuMM/OpenMM field on a real molecule, realized to `Stage::Acceleration` with forces present. Leg A deliberately runs at zero applied force (§4.5); it never exercises the per-atom→per-body force reduction that feeds the live reaction. Leg B's independent oracle is Simbody's `calcMobilizerReactionForcesUsingFreebodyMethod` (`SimbodyMatterSubsystem.h:2964`, public), which assembles applied body forces (`getRigidBodyForces`, i.e. the DuMM-reduced per-body force), constraint forces, and gyroscopic effects from scratch — a tip-to-base freebody sweep, ~3× slower and structurally independent of the PPlus method (`SimbodyMatterSubsystemRep.cpp:6105-6168`).

### 5.2 Placement and wiring (Recommended architecture)
Because the AMBER-prmtop entry point is Python-only (differential spec §4.2), a real molecule under the live field is only reachable through the clone's `World`/`Context` pybind path — the same path `_generate_molecule_oracle.py` uses. Therefore:

- Add an **additive, test-fixture-only** method on the clone `World` (e.g. `World::checkLiveFieldReactionResiduals()` returning a small POD/`py::dict` of residuals), which:
  1. realizes `worldState` through `Stage::Acceleration` under the live force field (`compoundSystem->realize(worldState, Stage::Acceleration)`), **independent of the sampler type** — it does **not** call `calcSpatialForces` and so is not subject to its OMMVV guard;
  2. computes `FM_standard` via `matter->calcMobilizerReactionForces(worldState, .)` and `FM_freebody` via `matter->calcMobilizerReactionForcesUsingFreebodyMethod(worldState, .)`, and returns `max_mbx ‖FM_standard[mbx] − FM_freebody[mbx]‖`;
  3. returns a whole-tree Newton-Euler balance residual (§5.3);
  4. echoes, for the `interestingMobodIndices` set, the exact tuple `calcSpatialForces` stores — `torque_G = FM_standard[mbx][0]`, `force_G = FM_standard[mbx][1]`, `u = getOneU(worldState,0)`, `uDot = getOneUDot(worldState,0)` — so the check also pins the storage slot mapping, not only the physics.
- Expose it via one pybind `.def` next to `dump_robotics_oracle_molecule` (`PyBind11.cpp:943`).
- Drive it from a clone-side script/pytest under `tests/fixtures/robotics_oracle_molecules/` (or the clone's own test tree) that builds `10ala` flexible (and optionally `1APQ`) exactly as the generator does, then asserts the residuals are below tolerance.

- **NOTE (gate placement).** This is a **clone-side** check: it needs the clone's `robo_bindings` built+installed (`cuda-release`), which the disasm `nox -s tests` gate does **not** build. It therefore lives with the fixture-generation tooling and runs in the clone's validation, **not** the disasm project gate. Give it a `test_`-collectable form if the clone has a pytest gate; otherwise a manual dev tool alongside `_generate_molecule_oracle.py`. Do not wire it into the disasm gate.
- **NOTE (rejected alternative).** A self-contained C++ test in the Simbody clone test tree (the `TestMobilizerReactionForces.cpp` pattern) built from hand-authored `MobilizedBody`s + a `Force` element would validate the reaction *math* under a generic force, but **not** the DuMM per-atom→per-body reduction on a real molecule — and that math is already covered by `TestMobilizerReactionForces.cpp` and by Scope-A §7.1 S3. Leg B's value is precisely the real-molecule live-DuMM path, so the pybind-driven form is preferred.

### 5.3 Whole-tree Newton-Euler balance (independent of both reaction methods)
With all reactions reported on the child at `Mo`, the reaction on the parent at `Mo` is the negative (`SimbodyMatterSubsystem.h:2393-2396`). Shifting each child's reaction to its parent's origin and summing, the net force system on every non-Ground body must equal its `Mk·A_GB − (appliedBodyForce + gyroscopicForce)` residual to ~`1e-9` relative — the same balance `TestReactionForces.cpp` (port) and the freebody sweep (clone) already embody. Compute it from `FM_standard`, `getBodyAcceleration`, `getRigidBodyForces`, and `getGyroscopicForce`; it is a third, method-independent oracle.

### 5.4 The `calcSpatialForces` message-bug (trivial)
`World::calcSpatialForces` (`World.cpp:4356-4358`) throws **when** the integrator **is** `OpenMMVelocityVerlet`, but the message reads "...is only implemented for OpenMMVelocityVerlet". The **condition is physically correct** — under OMMVV OpenMM integrates in Cartesian space, so the Simbody `u`/`uDot`/reactions are not the dynamics — but the message is inverted. Fix the message text only (e.g. "...is not supported for the OpenMMVelocityVerlet sampler; it requires an internal-coordinate/Simbody sampler"). No behavior change. Classify: Implementation note, not a correctness item.

### 5.5 `u`/`uDot` mapping (documentation)
For the 1-DOF Torsion mobilizers that are Scope-B's only flexibility, `getOneU(state,0) == bd.u[0]` and `getOneUDot(state,0) == bd.udot[0]` (both index `getFirstUIndex + 0`; dump `RoboticsOracleMoleculeDump.cpp:222-225`). `u` is a state variable (force-independent), so it matches between the zero-force dump and the live `calcSpatialForces`; `uDot` does **not** (the dump's `udot` is zero-force, `calcSpatialForces` reports the live-field `uDot`). This is a mapping note only — Leg A and `calcSpatialForces` are not expected to agree on `uDot`.

## 6. Correctness conditions

- **PRECONDITION (runtime guard, clone dump).** `PPlus` and `zPlusAll`/`A_GB_All` must be read after `matter.calcAccelerationForOracle(...)` at `Stage::Dynamics`+; the parent acceleration used in `APlus` is `A_GB_All[p]`, never `parent.getBodyAcceleration(state)`.
- **PRECONDITION (Leg B).** `worldState` is realized through `Stage::Acceleration` under the live field before either reaction method is called (both require Acceleration stage; `SimbodyMatterSubsystem.h:2437, 2948`).
- **INVARIANT (clone dump, leaf self-check — LEMMA-tagged tolerance).** For every **leaf** body (no children) at zero applied force, the reconstructed `FB_G(Bo)` equals the freebody value `getBodySpatialInertiaInGround(s)·A_GB_All[b] + getGyroscopicForce(s, b)` to `1e-9` relative. This cross-validates the PPlus reconstruction (right `APlus`, right parent, right `zPlus`) against the Newton-Euler formula **within the clone at zero force**, so it catches a reconstruction coding error without needing the force field switched off. It must fail if `A_GP` is sourced from the live accessor.
- **INVARIANT (Leg A, cross-engine).** On `10ala_regular` `random`, port `reactionBo` matches the clone reconstruction within the §4.4 tolerance; must fail on a `[torque;force]` slot swap or a dropped inward `Σ Φ·reac_child` term. Discriminating structure: it requires a **flexible tree at u≠0** (a rigid or u=0 case has zero reactions and cannot distinguish these bugs).
- **INVARIANT (Leg B).** `max_mbx ‖FM_standard − FM_freebody‖ ≤ 1e-9·max(1,‖·‖)` under the live field, and the §5.3 whole-tree balance residual ≤ `1e-9` relative; must fail if the per-atom→per-body DuMM force reduction is inconsistent with the realized accelerations.
- **LEMMA (Leg A teeth, expected value).** The `rest` state and the entire `10ala_rigid` case yield `reactionBo = reactionMo = 0` on both sides (u=0 ⇒ zPlus=0, udot=0, A_GB=0). This is the analytically-derivable expected value; it is a weak guard, and the test must label it as such so its green is not read as strong coverage.
- **LEMMA (`reactionMo` gate).** `reactionMo.linear == reactionBo.linear` exactly (the shift changes only the torque). `reactionMo.angular` is comparable cross-engine iff `Mo` agrees to `kMolStage1Tol`; otherwise it carries the same `X_GB.R()` convention ambiguity that made `X_GB.R()` inadmissible (differential spec §6 NOTE).

## 7. Touch list (components changed / conventions at risk)

- **Clone dump** (`Robosample/src/RoboticsOracleMoleculeDump.cpp`, `include/…hpp`): add the §4.1 reconstruction + leaf self-check; add four `SpatialVec`-worth of fields.
- **Clone pybind** (`Robosample/src/PyBind11.cpp`): marshal the four reaction fields; add the Leg-B `.def`.
- **Clone World** (`Robosample/include/World.hpp`, `src/World.cpp`): add the additive Leg-B method; fix the `calcSpatialForces` message (§5.4).
- **Generator** (`tests/fixtures/robotics_oracle_molecules/_generate_molecule_oracle.py`): emit reaction arrays; bump `schema_version`.
- **Loader** (`tests/RoboticsOracleMoleculeLoader.hpp`): read the four fields.
- **Port test** (`tests/TestRoboticsOracleMolecule.cpp`): add the reaction stage in `runMoleculeNumericDifferential` under the existing guards.
- **Clone-side driver** (new, alongside the generator): Leg-B pytest/dev tool.
- **Conventions at risk:** (1) `[torque; force]` slot order vs `[force; torque]`; (2) `Bo` (convention-free) vs `Mo`/`X_GB.R()` (convention-contaminated); (3) `shiftForceBy` sign `angular − r%linear`; (4) live vs zero-force parent `A_GB` in `APlus`; (5) `Φ` child-shift direction in the inward accumulation (parent→child origin vs child→parent).

## 8. Verification plan

- **Leg A regression discriminators (Normative):**
  - Slot swap: transpose `reactionBo.angular`↔`.linear` in the loader → `10ala_regular` `random` must go red. (Guards against a `[torque;force]` mismap.)
  - Live-`A_GP` bug: source `APlus`'s parent acceleration from `parent.getBodyAcceleration(state)` in the dump → the leaf self-check (§6 INVARIANT) must fail. (Guards the zero-force reconstruction contract.)
  - Sign of `shiftForceBy`: flip to `angular + r%linear` → the gated `reactionMo.angular` diff (where `Mo` agrees) must go red while `reactionBo` stays green. (Isolates the `Bo→Mo` shift.)
- **Leg B discriminators (Normative):** perturb the DuMM per-body force reduction (or realize at the wrong stage) → the freebody-vs-standard residual and the whole-tree balance must go red. The two-method agreement must hold at the clone's own `1e-9`-class tolerance on a real molecule, not only on the synthetic `TestMobilizerReactionForces.cpp` systems.
- **Non-redundancy check:** confirm the added assertions fail on a bug that `TestReactionForces.cpp` (port self-consistency) and `TestRoboticsOracle.cpp` Scope-A (synthetic cross-engine) both pass — the target bug is a **real-molecule** reduction/inward-accumulation error at u≠0, which neither existing suite exercises.
- **Determinism:** regeneration with the pinned `RNG_SEED` reproduces byte-identical reaction arrays (same contract as the existing dump, provenance §4/§5).

## 9. Open questions / blocking unknowns (none blocking)

- **`Mo` physical-point agreement across engines.** Whether the mobilizer outboard-frame origin `Mo` coincides cross-engine for matched Torsion bodies is not established here. It is **not blocking**: `reactionBo` is the primary anchor and is convention-free; the `reactionMo.angular` diff is explicitly gated on a runtime `Mo`-agreement assertion (§4.3), and a disagreement is reported as a finding rather than silently tolerated or silently compared. If a reviewer wants `reactionMo.angular` promoted to a primary anchor, the needed input is a proof (or an asserted fixture check) that the port's `X_BM` and Molmodel's outboard frame place `Mo` at the same physical atom for the Torsion mobilizer.
- **Leg B gate status.** The recommendation is clone-side, out of the disasm gate, because the disasm gate does not build clone bindings. If the project wants live-field reaction coverage *inside* `nox -s tests`, that requires either building the clone bindings in that session or a disasm-side reaction path under a live field — a larger scope decision for the orchestrator, not resolved here.
