# Robosample C++ Test Suite - Triage Record

Status: recovery draft, 2026-07-12. Companion to [`ARCHITECTURE.md`](ARCHITECTURE.md)
and [`MODULES.md`](MODULES.md). Scope: the 47 `Test*.cpp` and 13 `.hpp` helpers
under `tests/` (~20.2k lines). Python tests (`test_*.py`) are out of scope.

The test tree is evidence (it reveals intended contracts), oracle (it defines
behavior preservation for the source split), and refactor target (it is
AI-generated and bloated too). Under the two-phase freeze rule, tests move only
**after** the source splits are verified - this record is written now but its
`TEST-###` tickets execute in Phase B.

---

## 1. The cost model (decisive finding)

**Only one test - `TestAlchemy` - constructs a real OpenMM `Context`.** Every
other engine test drives forces through the pure-C++ `AnalyticForceBridge`
(harmonic wells, no OpenMM, no GPU). Consequences:

- The CMake GPU `RESOURCE_LOCK "gpu"` is nearly irrelevant to wall-clock; it
  matters only for `TestAlchemy`.
- Real cost is **CPU HMC/NCMC sampling-chain length**, not GPU work.
- A slow/fast CI split needs **no CMake label system and no relinking** - the
  boundary already exists as the `ROBOSAMPLE_SLOW_TESTS` environment variable,
  read by `StatTest.hpp::slowEnabled()` (line 35). A bare `ctest` `GTEST_SKIP()`s
  the million-sample cases; `nox -s tests` sets the variable.

The gate is applied **inconsistently**: 8 files honor it, but 5 heavy files run
their full chain unconditionally (see section 4).

---

## 2. FAST / SLOW classification - 34 FAST, 13 SLOW

Signal: *algebra* = math/helpers only, no sampling; *single-shot* = build + one
eval, no loop; *chain* = HMC/NCMC/integrator loop of N draws.

### SLOW (13) - the CI-split and gating target

| File | Lines | Heaviest signal | Currently gated? |
|---|---|---|---|
| TestEnsembleValidation | 260 | N = 400k / 1.2M / 1.5M draws | yes |
| TestEnsembleOrientation | 299 | 600k moves x3 | yes |
| TestNcmcTeleport | 533 | N = 400k / 2M / 4M | yes |
| TestFixmanBoltzmann | 237 | `drv.move()` chain | yes |
| TestMassScaleInvariance | 248 | chain + reversibility | yes |
| TestNcmcExplicitSolvent | 730 | slow gate; **3/8 permanent `GTEST_SKIP` stubs, 6/8 skip at default tier** | yes |
| TestTwoRobotContact | 643 | nDraws = 300k Fixman contact | **no** |
| TestCyclicBoltzmann | 242 | 3000 `drv.move()` | **no** |
| TestIntegrator | 781 | 200-1500-step chains x10 | **no** |
| TestFreeJointKEPump | 731 | nSteps = 4000 chains | **no** |
| TestAlchemy | 539 | only real OpenMM `Context`; GPU lock | **no** |
| TestEquipartition | 276 | HmcDriver + StatTest KE sampling | partial |
| TestNCMCWork | 915 | **mixed**: tests 1-8 algebra (FAST), 9-16 work chains (SLOW) | partial |

### FAST (34) - run on every push

Pure algebra / single-shot, sub-second-ish: TestTransform, TestGeometry,
TestOrientation, TestStability, TestGibbsWorlds, TestVectorMath,
TestBatAnchorInvolution, TestSpatialAlgebra, TestRotationConstruction,
TestLinearAlgebra, TestBuilders, TestConstraintSolver, TestAnalyticForce,
TestMobilizerKinematics, TestReverseMobilizer, TestTransfer, TestInertia,
TestNMALinearAlgebra, TestLinearAlgebraOracle, TestKineticEnergy,
TestJointKernels, TestBatScalingJacobian, TestPeriodicBoundary, TestMobilizer,
TestQuaternion, TestReactionForces, TestRexAcceptanceAlgebra, TestAtomTransfer,
TestBiasForces, TestConstraints, TestMassMatrix, TestFixmanIdealizedChains,
and the two oracle files below.

**FAST\* (heavy fixture volume, still no sampling):** TestRoboticsOracle (1283
lines, 44 tests, 85 `loadCase`, ABA per case) and TestRoboticsOracleMolecule
(1171 lines, builds a full `World`/`Context` per molecule fixture). Single-shot
per case but bulk fixture I/O; candidates for the case-family split in section 5.

---

## 3. Triage by role (behavior-preservation classes)

- **Contract tests** (feed the invariant list, survive untouched): the robotics
  oracles (differential vs. Simbody fixtures -> INV-8-adjacent kinematic
  contracts), TestMassMatrix / TestKineticEnergy (mass-metric, INV-5),
  TestConstraints / TestConstraintSolver (INV-6), TestBatScalingJacobian /
  TestBatAnchorInvolution (INV-7), TestRexAcceptanceAlgebra (INV-8),
  TestReactionForces (INV-1), the ensemble/Boltzmann suite (detailed balance).
- **Characterization tests** (pin current behavior, label as such): the
  integrator drift/reversibility chains (TestIntegrator, TestFreeJointKEPump),
  TestStability, TestCyclicBoltzmann.
- **Implementation-coupled** (will break under pure code motion; each needs a
  predicted resolution in the SPLIT ticket that causes it):
  `ConstraintTestAccess.hpp` friends into the private constraint solver
  (`solveSmallSpd`/`solveCoupling`) - retarget when Constraints splits.
  (CORRECTED 2026-07-12, from ticket authoring - two earlier claims here were
  wrong: `TestRexAcceptanceAlgebra` does **not** call `Context::attemptREXSwap`;
  it reimplements the acceptance formulas in-file and never instantiates a
  `Context`, so `C4` has no live breakage - README CX-7. And
  `RobotLinearAlgebra.hpp` is a **deliberately independent** oracle, not a
  duplicate to retarget: `TestLinearAlgebraOracle` diffs the engine against it, so
  pointing it at the extracted `hinge_linalg` would make the comparison
  tautological - it SHALL stay independent - README CX-5.)
- **Stub / vacuous** (deletion candidates, human-gated): `TestNcmcExplicitSolvent`
  has 3 permanently-deferred `GTEST_SKIP` stubs (INV1/INV2/INV3) awaiting
  OpenMM-backed fixtures, and 6/8 skip at the default tier; `tests/ForceBridge.hpp`
  is a 10-line near-dead no-op stub colliding in name with `AnalyticForceBridge.hpp`.

---

## 4. Duplication (the AI-generated tax)

- **`slowEnabled()` copied 5x.** TestFixmanBoltzmann:63, TestEnsembleOrientation:82,
  TestMassScaleInvariance:48, TestEnsembleValidation:48, TestNcmcTeleport:46 each
  redeclare a local copy of the shared `StatTest.hpp::slowEnabled()`. Only
  TestNcmcExplicitSolvent imports the shared one.
- **`kT300` / `kB` literal (`0.0083144626 * 300.0`) copied 10x** across the
  statistical tests.
- **Six ad-hoc force-bridge subclasses** reinvent `AnalyticForceBridge`:
  `ZeroBridge`, `NcmcLambdaBridge`, `TwoRobotBridge`, `InfForceBridge`,
  `LambdaWellBridge`, `LambdaAnalyticBridge` - all lambda/scaled harmonic variants
  with no shared base.
- **NCMC step-loop harness duplicated** (`runNcmcLoop`/`NcmcRun`) between
  TestNCMCWork and TestNcmcExplicitSolvent; the HMC `move()` chain is re-inlined
  in every ensemble test.
- **Zero gtest fixtures across all 47 files** - every test rebuilds
  model/state/bridge in its body. This is the primary copy-paste source.

---

## 5. Refactoring seams (Phase B `TEST-###`)

1. **Consolidate the slow gate.** Delete the 5 local `slowEnabled()` copies; add a
   `SlowStatTest` gtest fixture whose `SetUp()` `GTEST_SKIP()`s unless
   `ROBOSAMPLE_SLOW_TESTS`. Collapses ~15 repeated skip blocks to inheritance.
2. **Gate the 5 ungated slow files** (Cyclic, Integrator, FreeJointKEPump,
   TwoRobotContact, Alchemy) so every-push CI runs FAST + smoke, and the full slow
   tier runs nightly. No CMake changes required.
3. **`TestPhysConstants.hpp`** for `kB` / `kT300` / standard temperatures - removes
   the 10x duplication.
4. **`HarmonicBridge<Policy>` base** consolidating the six ad-hoc `*Bridge`
   subclasses; lambda-scaling and inter/intra split as policy parameters over
   `AnalyticForceBridge`. Retire `tests/ForceBridge.hpp` into it as the explicit
   "no force" policy.
5. **`SamplingHarness` helper** - `runHmcChain(drv, nMoves, stride) -> Marginals`
   and `runNcmcLoop(...) -> NcmcRun` - so ensemble/NCMC tests differ only in the
   assertion.
6. **`RobotFixture` base** holding a built `RobotModel` + `RobotState` +
   `AnalyticForceBridge` (parameterized by a builder lambda) - removes the
   ~20-file build prologue.
7. **Split the two oracle monsters** along the case families the loader already
   partitions (single-state / body-output / aggregate / fuzz / multi-system) into
   ~300-line binaries that parallelize under ctest.
8. **Split `TestNCMCWork`** into a FAST `protocolLambda` algebra file and a SLOW
   work-chain file (it currently mixes both, defeating file-granularity CI split).
9. **Target structure:** mirror the module map one-to-one
   (`src/world/HmcMove.cpp <-> tests/world/HmcMoveTest.cpp`), one fixture per class,
   behavior-sentence test names, shared helpers under `tests/support/` held to
   production header discipline.

---

## 6. Decisions (recorded 2026-07-12)

- **D-T1 - SUPERSEDED 2026-07-12: slow tests are impracticable; optimize them.**
  The earlier "document only" is replaced by an explicit directive to make the
  slow tier practical. Optimization plan in section 8. This becomes its own Phase-B test
  track, gated by the mutation-power criterion (section 8.5).
- **D-T2** `TestNcmcExplicitSolvent` has 3 permanently-deferred stubs (6/8 skip at
  default tier) - write the OpenMM
  fixtures, keep as a documented placeholder, or delete? (Open; human-gated;
  deletion weakens the oracle.)
- **D-T3** Retire `tests/ForceBridge.hpp` dead stub now, or after the
  `HarmonicBridge` consolidation? (Open.)
- **D-T4** Split the two 1.2k-line oracle files, or leave them as monolithic
  contract fixtures? (Open; they are the crown-jewel differential oracle;
  splitting is navigability-only, no behavior change.)

## 7. REX oracle coverage (gap, recorded 2026-07-12)

The existing replica-exchange test coverage is **algebraic, not behavioral**:

- `tests/TestRexAcceptanceAlgebra.cpp` (7 tests) - the swap *acceptance formula*
  in isolation via `Context::attemptREXSwap`: detailed balance of the acceptance
  ratio, Jacobian sign-flip blindness, RENEMC `ETerm_nonequil`, domain-error
  sentinel. Theory-derived, correct - but tests the function, not the driver.
- `tests/test_rex_label_swap_equivalence.py` - the label-swap `RunREX` vs
  coordinate-swap `runREX` INVARIANT-EQUIV equivalence (Python).
- `tests/test_rex_swap_acceptance_algebra.py` - Python mirror of the acceptance
  algebra (spec V8).

**Missing:** an end-to-end **stationary-distribution / detailed-balance oracle**
for the full REX chain - that the driver, wired to a correct acceptance formula,
actually samples each replica's Boltzmann distribution at its temperature (spec
INV-1/INV-2) and leaves the product distribution invariant. Acceptance-algebra
correctness does not prove the driver wires it correctly. This gap is the
behavioral oracle a Context-REX split needs before it moves. Scope note: the
runnable surface today is REMC + Default; RENE/REBASONTOP driven-REX is
uncompiled (see ARCHITECTURE OQ-5) and RENEMC's round-loop throws, so a behavioral
oracle can only cover REMC/Default until that code is built. A researcher-derived
spec for this oracle is a separate behavioral track, sequenced *before* the
Context REX code motion.

---

## 8. Slow-tier optimization plan (2026-07-12)

Why the 13 slow tests are slow: they are **statistical-convergence oracles**.
Each draws 10^5-10^6 HMC/NCMC samples, histograms an observable, and chi-squares
it against a closed-form Boltzmann/Gamma target at `alpha = 1e-4`. `N` is large for
two reasons: binned goodness-of-fit at small `alpha` needs many samples per bin, and
HMC samples are **autocorrelated** so the effective sample size is far below `N`
(the Fixman campaign already traced slow-tier trouble to autocorrelation, not the
engine). `TestIntegrator` / `TestFreeJointKEPump` are different: they run long
single trajectories to measure energy drift / reversibility - length is the
assertion, not a convergence artifact.

Levers, ordered by value:

### 8.1 Parallelize wall-clock (zero statistical change)
- **Drop the GPU lock from the 12 CPU-only slow tests.** The CMake wiring stamps
  `RESOURCE_LOCK "gpu"` on *every* robosample test under `USE_CUDA`, but only
  `TestAlchemy` builds a real OpenMM context (see section 1). The other 12 are pure-CPU
  analytic-bridge tests being needlessly serialized. Keep the lock on
  `TestAlchemy` alone; the rest run concurrently under `ctest -j`.
- **Intra-test OpenMP-parallel independent chains.** The analytic force bridge is
  stateless, so `K` independent HMC chains (each its own `RobotState` + RNG
  stream) run on `K` threads and pool into one histogram - same total samples,
  ~`Kx` less wall-clock. Highest-leverage per-test change.

### 8.2 Reduce total samples for equal power
- **Size `N` to a target effective sample size**, not a hardcoded raw count:
  estimate the integrated autocorrelation time and draw to a fixed ESS. Tuning
  HMC step size / `mdSteps` to maximize independent-samples-per-move then cuts
  wall-clock directly.
- **Use a more powerful statistic than binned chi-square.** For 1-D marginals,
  Anderson-Darling / KS need far fewer samples for equal power; moment checks
  (equipartition `<U> = (n/2)RT`, `<U^2>`) converge as `1/sqrt(N)` with small constants
  and catch most bias cheaply. Reserve the fine-histogram chi-square for nightly.

### 8.3 Three-tier gate (fixes the shared-function point)
Replace the binary `ROBOSAMPLE_SLOW_TESTS` with **one** shared function returning
a tier, read by all tests; delete the 5 duplicated `slowEnabled()` copies (section 4):

- **Smoke** (always on, seconds): moment checks only.
- **Standard** (every-push CI, ~1 min, threaded): coarse histogram + AD.
- **Exhaustive** (nightly): the current millions-sample fine chi-square.

Each test derives `N` from the tier instead of hardcoding it. This is the
concrete answer to "they all share the same function that checks the environment":
make it genuinely one function returning `{tier, sampleBudget}`, not five copies
of a boolean, and thread `N` through it.

### 8.4 Smaller systems where the distribution is still exact
Some ensemble tests can use fewer DOF (the closed-form target holds for any
positive-definite system) to cut per-move cost without weakening the check.

### 8.5 The binding constraint: mutation-power gates every shortening
These tests are the **oracle** for the source refactor; changing `N` or the
statistic changes what they verify. So every shortening SHALL be gated by a
**mutation test**: inject a known bias into the engine (a wrong Jacobian sign, a
swapped beta, a dropped Fixman term) and confirm the shortened oracle still fails
reliably at the same false-fail rate. Shorten `N` only as far as the injected
bias is still caught. If a 10x smaller `N` still fails on the known-biased
engine, the shortening is free; if not, it is over-shortened. This keeps
"impracticable but rigorous" from turning into "fast but blind."
