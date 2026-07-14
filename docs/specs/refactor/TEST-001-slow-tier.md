# TEST-001: Make the slow statistical tier practical

Status: draft, 2026-07-12. Phase-B test ticket (source frozen). Executor: `coder`. Style: `styles/spec.md`. Companion to
[`../../architecture/TESTS.md`](../../architecture/TESTS.md) section 8 (the plan this
ticket executes) and [`VERIFY.md`](VERIFY.md) section 5 (exit criteria). This is a
test-only change; it touches no engine source. It changes what a subset of tests
*cost*, and - where it changes what they *verify* - is gated by the mutation-power
criterion in section 6.

## Context (why)

The 13 slow tests (`TESTS.md section 2`) are statistical-convergence oracles: each draws
10^5-10^6 autocorrelated HMC/NCMC samples, histograms an observable, and
chi-squares it against a closed-form Boltzmann/Gamma target at `alpha = 1e-4`. They are
"basically impracticable" (user, 2026-07-12). Root causes: binned goodness-of-fit
at small `alpha` is sample-hungry; HMC samples are autocorrelated so effective sample
size << N; and all robosample tests are serialized under a GPU resource lock even
though only `TestAlchemy` uses the GPU (`TESTS.md section 1`).

Ownership boundary: these tests are the **oracle** for the source refactor, so any
change to `N` or the statistic changes what they verify. This ticket separates the
free wins (no statistical change) from the power-affecting wins (mutation-gated).

## Scope

- **In:** `tests/StatTest.hpp` (the shared gate + statistics), the 13 slow
  `Test*.cpp`, the CMake `RESOURCE_LOCK` wiring for the CPU-only slow tests, one new
  `tests/support/` helper for parallel chains, and one new mutation-fixture header.
- **Out:** any engine source; the FAST tests; `TestAlchemy`'s GPU lock (kept).

## Deliverables (ordered; the first two are free, the rest are gated)

### D1 - One shared tier gate (fixes the duplicated function)
Replace the binary `ROBOSAMPLE_SLOW_TESTS` check with a single function in
`StatTest.hpp`:

```
enum class TestTier { Smoke, Standard, Exhaustive };
TestTier testTier();          // reads ROBOSAMPLE_TEST_TIER (0/1/2); default Smoke.
                              // ROBOSAMPLE_SLOW_TESTS=1 maps to Exhaustive (back-compat).
long tierSampleBudget(long smoke, long standard, long exhaustive);  // pick N by tier
```

Delete the **5 local `slowEnabled()` copies** (`TestFixmanBoltzmann:63`,
`TestEnsembleOrientation:82`, `TestMassScaleInvariance:48`, `TestEnsembleValidation:48`,
`TestNcmcTeleport:46`); route every gate through `testTier()`. Add a `SlowStatTest`
gtest fixture whose `SetUp()` `GTEST_SKIP()`s the case if its required tier exceeds
`testTier()`, collapsing the ~15 repeated skip blocks (`TESTS.md section 4`) to
inheritance. Back-compat: `ROBOSAMPLE_SLOW_TESTS=1` keeps meaning "run the heavy
path" (-> Exhaustive), so `nox -s tests` is unaffected until its call site is updated.

### D2 - Parallelize wall-clock (no statistical change)
- **CMake:** remove `RESOURCE_LOCK "gpu"` from the 12 CPU-only slow tests; keep it
  on `TestAlchemy` alone. They then run concurrently under `ctest -j`. (The
  `robosample_add_test` helper currently stamps the lock on every test when
  `USE_CUDA`; make the lock opt-in per test, `TestAlchemy` being the only opt-in.)
- **Intra-test:** add `tests/support/ChainPool.hpp` - run `K` independent HMC
  chains, each with its own `RobotState` + independent RNG stream (seed offset),
  on OpenMP threads, and pool their samples into one histogram/accumulator. The
  analytic force bridge is stateless, so chains share the read-only `RobotModel`
  and do not race. Same total samples, ~`Kx` less wall-clock. This is a **pure
  refactor of how samples are gathered** (same seeds, same total draws, same
  statistic) and SHALL leave Exhaustive-tier results statistically identical.

### D3 - Cheap Standard-tier statistics (power-affecting -> gated by section 6)
For each slow test define three tier bodies:
- **Smoke** (always on, seconds): moment checks only - equipartition `<U> = (n/2)RT`
  and `<U^2>` within `4*stderr` (reuse `MeanAccumulator`). Catches gross bias.
- **Standard** (every-push CI, target <= ~1 min each, threaded via D2): a coarse
  histogram + Anderson-Darling / KS on the 1-D marginal (fewer samples than a fine
  binned chi-square for equal power), at a target *effective* sample size (section D4).
- **Exhaustive** (nightly): the current millions-sample fine chi-square at
  `alpha = 1e-4`, unchanged - this is the reference the shorter tiers are validated
  against.

`N` per tier comes from `tierSampleBudget(...)`, not a hardcoded literal.

### D4 - Size N to effective sample size (power-affecting -> gated by section 6)
Where a test currently hardcodes a raw draw count, replace it with a target
*effective* sample size: estimate the integrated autocorrelation time of the
observable and draw until the ESS target is met (cap at the Exhaustive budget).
This removes the "guess a big N" pattern the Fixman campaign flagged and lets the
Standard tier hit its power target with far fewer wall-clock draws. NOTE: tests
that measure trajectory-length physics - `TestIntegrator`, `TestFreeJointKEPump`
(energy drift / reversibility over a fixed number of steps) - are NOT
sample-count oracles; their step counts are the assertion and SHALL NOT be reduced.
They get D2 parallelism (independent seeds concurrently) only.

## Predicted test breakage
None beyond the intended edits: this is a Phase-B ticket, source is frozen, and the
FAST tests and `TestAlchemy` are untouched except `TestAlchemy` keeps its lock.

## Exit criteria (machine-checked; see VERIFY.md section 5)

- Build clean; **no test includes a `.cpp`**; three consecutive full-suite runs
  identical at each tier.
- **Exhaustive-tier pass set + assertion counts identical to the B1/B2 baseline**
  (the D1/D2 changes are behavior-neutral at Exhaustive). Any Standard/Smoke-tier
  count delta is recorded as an explicit, human-approved triage delta in this
  ticket.
- **Coverage map (B5) not reduced.**
- The 5 duplicated `slowEnabled()` copies are gone; one `testTier()` remains.
- `ctest -j` shows the 12 CPU-only slow tests running concurrently (GPU lock only
  on `TestAlchemy`), with a measured wall-clock reduction reported in the ticket.
- **section 6 mutation-power gate passes for every shortened oracle.**

## 6. The mutation-power gate (the binding constraint on D3/D4)

A shortened oracle is accepted only if it still catches a known bias. For each test
whose `N` or statistic changes (D3/D4), add a `tests/support/BiasedBridge.hpp`
mutation fixture that injects one documented bias into the analytic force / swap /
correction path - e.g. a wrong BAT-Jacobian sign, a swapped beta in the swap
acceptance, a dropped Fixman term (mirroring the corruptions already in
`TestRexAcceptanceAlgebra` and the REX oracle spec's O8). The Standard-tier oracle
run against the biased engine SHALL **fail** at the same false-fail budget the
Exhaustive tier has. `N` is shortened only as far as the injected bias is still
caught reliably (SHOULD: bias detected in >= 19/20 fixed-seed repeats). A shortening
that lets a known bias pass is rejected - "fast but blind" is not acceptable for an
oracle. This gate is what lets us make the tier practical without weakening the
refactor's behavior-preservation guarantee.
