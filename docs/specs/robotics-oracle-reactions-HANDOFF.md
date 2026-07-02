# HANDOFF: robotics-oracle reaction-force validation

Date: 2026-07-02. Branch: `disasm`. Status: **implementation complete + verified green in isolation; full `nox -s tests` gate NOT confirmed (pytest phase hung on unrelated pre-existing sim tests; killed).**

## What the task was

Verify whether the values `World::calcSpatialForces()` stores (`Robosample/src/World.cpp:4355`) are part of the robotics-oracle campaign, and rigorously test them (analytical + numerical).

**Finding (verified from clone source):** `calcSpatialForces` stores, per interesting mobilized body, `torque_G=reactionForces[mbx][0]`, `force_G=reactionForces[mbx][1]` from `matter->calcMobilizerReactionForces(worldState,·)`, plus `u`/`uDot`. `SimbodyMatterSubsystemRep::calcMobilizerReactionForces` (`SimbodyMatterSubsystemRep.cpp:6076-6094`) is **bit-identical** to `findMobilizerReactionOnBodyAtMInGround` (`MobilizedBody.cpp:406-414`) = the oracle's `reactionMo` (SpatialVec `[angular=torque; linear=force]` at outboard frame origin Mo, in Ground).

**The gap that was closed:** the Scope-A synthetic oracle already pinned `reactionBo`/`reactionMo` (1e-8), but (a) the Scope-B **molecule** dump omitted reactions entirely, and (b) no fixture ran reactions under a **live** force field — the exact conditions `calcSpatialForces` runs in.

Authoritative spec: `docs/specs/robotics-oracle-reactions.md`. Parent: `docs/specs/robotics-oracle-differential.md` (§4.5 cyclic carve-out, §6 frame-convention NOTE), `docs/specs/singular-dof-fixman.md`.

## What was implemented (two legs, per user choice "both")

### Leg A — Scope-B cross-engine reaction pin (in the disasm gate)
Reconstructs the **zero-applied-force** reaction inside the molecule dump from buffers `calcAccelerationForOracle` already returns (`zPlus`, `A_GB`) + the `PPlus` getter — no new Simbody patch. `FB_G = zPlusAll[mbx] + PPlus[mbx]·APlus`, `APlus=(A_GP[0], A_GP[1]+A_GP[0]%p_PB_G)` with `A_GP=A_GB_All[parentMbx]` (zero-force, NOT the live accessor), then shift Bo→Mo `[F.ang − p_BM_G % F.lin ; F.lin]`.

Files changed (all under the **clone** `Robosample/` and the **port** `tests/`):
- `Robosample/include/RoboticsOracleMoleculeDump.hpp` — added `reactionBo_ang/lin`, `reactionMo_ang/lin` to `BodyDump`.
- `Robosample/src/RoboticsOracleMoleculeDump.cpp` — the reconstruction in the per-body loop (~line 240, after the A_GB block).
- `Robosample/src/PyBind11.cpp` — marshals the 4 reaction fields into the dump dict (~line 113).
- `tests/fixtures/robotics_oracle_molecules/_generate_molecule_oracle.py` — emits the 4 arrays per state/body; **schema_version bumped 1→2**.
- `tests/RoboticsOracleMoleculeLoader.hpp` — reads the 4 fields into `MoleculeOracleBody`.
- `tests/TestRoboticsOracleMolecule.cpp` — new "Stage 4b" in `runMoleculeNumericDifferential`: port reaction via `RobotEngine::calcMobilizerReactionForces` (independent Newton-Euler method) compared to the clone reconstruction. Compares `reactionBo` (ang+lin) AND `reactionMo` (ang+lin, N1) at magnitude-relative tol `stage4Tol*max(1,‖ref‖)`; per-body phantom-angular skip; under the existing whole-tree singular-skip guard. N3 non-vacuousness guard `EXPECT_GT(maxRefReactionNorm,1)` on the random state.

### Leg B — live-force-field check (CLONE-side, OUT of the disasm gate)
`Robosample/src/PyBind11.cpp` — free fn `check_live_field_reactions(World&)` (~line 130): realizes `worldState` to `Stage::Acceleration` under the LIVE field, cross-checks `matter.calcMobilizerReactionForces` (the `calcSpatialForces` path) vs the independent `matter.calcMobilizerReactionForcesUsingFreebodyMethod`; returns max abs/rel diff. Registered as `.def("check_live_field_reactions",…)`.
Driver: `tests/fixtures/robotics_oracle_molecules/_check_live_field_reactions.py` (manual, `_`-prefixed → not pytest-collected; needs clone bindings built). Asserts rel diff ≤ 1e-9 and reaction norm > 1.

### Incidental fix
`Robosample/src/World.cpp:4356` — inverted throw message in `calcSpatialForces` corrected (condition unchanged; it correctly throws for OpenMMVelocityVerlet).

## Build / fixture state (already done this session)
- Clone rebuilt + installed: `cmake --build /home/victor/Robosample/Robosample/build/cuda-release --parallel && cmake --install /home/victor/Robosample/Robosample/build/cuda-release` (installs `robo_bindings.so` into `Robosample/python/robosample/`). NOTE: `-j0` does NOT work with `cmake --build <dir>` (only with `--preset`); use `--parallel`.
- Fixtures regenerated (schema v2, with reactions): `10ala_rigid`, `10ala_regular`, `1APQ_cyclic` `.moldyn.npz`+`.manifest.json` via `python3 tests/fixtures/robotics_oracle_molecules/_generate_molecule_oracle.py`.
- Port `cuda-tests` built (`cmake --build --preset cuda-tests --target TestRoboticsOracleMolecule`).

## What is VERIFIED green
- `ctest --test-dir build/cuda-tests -R RoboticsOracleMolecule` → **9/9 pass** (including the new reaction stage with N1 `reactionMo.angular` + N3 guard).
- Non-vacuous: regular/random max ‖reactionBo_lin‖≈209.5; regular/rest=0; minEigD=0.032 ≫ 1e-12 (not singular-skipped).
- Teeth confirmed: negating one reference component turned `RegularFlexibleNumeric` RED, then reverted GREEN.
- N1 confirmed: `reactionMo.angular` passes cross-engine → the Mo frame origin coincides cross-engine (the shipped `calcSpatialForces` slot is now pinned).
- Leg B: live-field two-method rel diff **1.1e-14** (10ala) / **5.9e-14** (1APQ), reaction norms 2320/1601 (nonzero).
- Hostile reviewer: **no blocking findings**; reconstruction verified bit-faithful. N1/N3 were its follow-ups and are now DONE. S1 (spec §6 leaf self-check) NOT implemented — see below.

## WHAT NEEDS TO BE DONE NEXT

1. **Full C++ ctest suite regression run** — was interrupted (killed) partway. Run `ctest --test-dir build/cuda-tests --output-on-failure -j0` and confirm the WHOLE suite is green (only `TestRoboticsOracleMolecule.cpp` + its loader header changed on the port side, so a regression elsewhere is unlikely, but confirm). **Do NOT pipe through `tail` — it buffers and hides interim output;** stream it or write to a file and poll counts.
2. **The `nox -s tests` gate** — its **pytest phase hung ~2h+** on pre-existing MD-simulation tests (`test_fixman_idealized_chains`, `test_fixman_proline_realism`, `test_ensemble_pe_ladder`, `test_torsion_conformational`, `test_installation`). These are UNRELATED to this change (no collected Python test reads the reaction fixtures; the only Python touched is the `_`-prefixed manual generator/driver). Determine whether this pytest hang is environmental (GPU contention with concurrent builds) or a genuine pre-existing hang, independent of this work. The C++ half of the gate is the relevant check for this change and is green.
3. **S1 (reviewer should-fix): spec §6 leaf self-check not implemented.** The spec asks for a clone-local Newton-Euler self-check in the dump (`FB_G == getBodySpatialInertiaInGround(s)·A_GB_All[b] + getGyroscopicForce(s,b)` for leaf bodies, 1e-9). It is subsumed by Leg A's cross-engine `reactionBo` diff (an independent method that catches the live-`A_GP` bug even at rest), so it was left out. **Either** implement it in `RoboticsOracleMoleculeDump.cpp` (needs a clone rebuild + fixture regen), **or** add a NOTE to spec §6 documenting it as redundant with Leg A. (Recommend: document, unless a third independent check is wanted.)
4. **N2 (reviewer nit, low priority):** Leg B implements only the two-method agreement; spec §5.3 whole-tree Newton-Euler balance and §5.2 storage-slot echo (`torque_G/force_G/u/uDot` tuple) were dropped. Both Simbody methods share `getRigidBodyForces`, so Leg B validates reaction-operator self-consistency but not independently the DuMM per-atom→per-body reduction. Optional to strengthen; note the coverage nuance in the spec if not.
5. **Spec sync:** update `docs/specs/robotics-oracle-reactions.md` §4.3 to reflect that `reactionMo.angular` IS now compared cross-engine (Mo coincidence verified empirically), superseding the earlier "gated/not-diffed" recommendation. Mark spec Status → implemented once the above close.
6. **Commit** (user has not been asked yet; only commit when the user asks). The `disasm` branch has substantial pre-existing uncommitted work; scope any commit to THIS change's files (listed above + the regenerated fixtures + the two new docs).

## Key commands (CLAUDE.md-compliant, no cd / no -j0-on-bare-build)
- Clone build+install: `cmake --build /home/victor/Robosample/Robosample/build/cuda-release --parallel` then `cmake --install /home/victor/Robosample/Robosample/build/cuda-release`
- Regen fixtures: `python3 tests/fixtures/robotics_oracle_molecules/_generate_molecule_oracle.py`
- Port test build: `cmake --build --preset cuda-tests --target TestRoboticsOracleMolecule`
- Run: `ctest --test-dir build/cuda-tests -R RoboticsOracleMolecule --output-on-failure`
- Leg B driver: `python3 tests/fixtures/robotics_oracle_molecules/_check_live_field_reactions.py`
