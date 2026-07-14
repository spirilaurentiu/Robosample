# TEST-004: Promote the copy-pasted sampling loops into a `SamplingHarness`

Status: draft, 2026-07-12. Phase-B test ticket (source frozen). Executor: `coder`. Style: `styles/spec.md`. Companion to
[`../../architecture/TESTS.md`](../../architecture/TESTS.md) section 4 (the duplicated
`runNcmcLoop`/`NcmcRun` and the re-inlined HMC `move()` chain) and section 5.5 (the seam
this executes), and [`VERIFY.md`](VERIFY.md) section 5. Test-only motion.

## Context (why)

`TESTS.md section 4` records two duplicated sampling loops:

- **The NCMC step loop.** `NcmcRun` + `runNcmcLoop(...)` is defined in
  `TestNCMCWork.cpp:348` (struct) / `:361` (function) and re-implemented in
  `TestNcmcExplicitSolvent.cpp` (its `ncmcMoveConstructionI/II` +
  `innerGhmcStep` + `seedMomenta` machinery, `TestNcmcExplicitSolvent.cpp:466-694`,
  which mirror the same PERTURB/PROPAGATE structure). `runNcmcLoop` is called at
  seven sites inside `TestNCMCWork.cpp` alone (`:433,:478,:511,:769,:774,:818,:823`).
- **The HMC `move()` chain.** The `for (i<N) drv.move(); accumulate` loop is
  re-inlined in every ensemble test: `TestFixmanBoltzmann.cpp`,
  `TestMassScaleInvariance.cpp`, `TestEnsembleOrientation.cpp`,
  `TestCyclicBoltzmann.cpp`, `TestEnsembleValidation.cpp` (all call
  `drv.move()`), each histogramming or accumulating an observable with
  `StatTest.hpp`'s `Histogram`/`MeanAccumulator`. The move itself already lives in
  the shared `tests/HmcDriver.hpp` (`HmcDriver<Bridge>::move()`,
  `HmcDriver.hpp:180`); only the *loop around it* is duplicated.

This ticket lifts the two loops into one `tests/support/SamplingHarness` so the
ensemble/NCMC tests differ only in the assertion, per `TESTS.md section 5.5`.

## Scope

- **In:** one new `tests/support/SamplingHarness.hpp` exposing
  `runHmcChain(...) -> Marginals` and `runNcmcLoop(...) -> NcmcRun`; conversion of
  the call sites that inline these loops.
- **Out:** any engine source; `HmcDriver::move()` itself (unchanged - the harness
  calls it); the statistics in `StatTest.hpp`; the tier gate (TEST-001); the
  fixture prologue (TEST-002); the bridges (TEST-003).

## Moves (exactly what)

### New: `tests/support/SamplingHarness.hpp`

- `struct Marginals` - the pooled observable container the ensemble tests build:
  a `stat::MeanAccumulator` and one or more `stat::Histogram`s
  (`StatTest.hpp:149,190`), plus the driver's `attempted()`/`accepted()` counts
  (`HmcDriver.hpp:84-90`). It holds sample summaries, not per-sample state.
- `template <class Bridge> Marginals runHmcChain(HmcDriver<Bridge>& drv, long
  nMoves, int stride, ObserverFn observe)` - the shared chain loop: call
  `drv.move()` `nMoves` times, invoke `observe` every `stride` accepted/attempted
  step to push the test's scalar observable into the `Marginals`. The loop body
  SHALL be the exact `for`-loop the ensemble tests inline today (same move count,
  same stride, same accept/attempt bookkeeping) so the drawn stream is
  bit-identical.
- `struct NcmcRun` and `runNcmcLoop(...)` - **moved verbatim** from
  `TestNCMCWork.cpp:348-408` into the harness namespace. The struct fields
  (`work`, `Hstart`, `Hend`, `heat`, `ok`) and the PERTURB-at-fixed-q /
  PROPAGATE-one-Verlet-step body (`TestNCMCWork.cpp:377-407`, including the
  `ncmc::protocolLambda` schedule call) SHALL be reproduced unchanged. Templatize
  the bridge parameter (currently hard-typed `LambdaAnalyticBridge&`,
  `TestNCMCWork.cpp:363`) over the TEST-003 `HarmonicBridge<IntraLambdaInterPolicy>`
  so `TestNcmcExplicitSolvent.cpp`'s construction path can share it - that
  templatization is the one shape change, and it SHALL leave `TestNCMCWork.cpp`'s
  instantiation identical.

### Call-site conversions
- The five `runNcmcLoop` sites in `TestNCMCWork.cpp` include the harness and drop
  the local definitions; the numeric arguments are unchanged.
- Each ensemble test replaces its inline `drv.move()` loop with a `runHmcChain`
  call carrying its observer lambda; the assertion on the returned `Marginals`
  stays in the test body.

## Target structure (after this ticket)

```
tests/support/
  SamplingHarness.hpp   # runHmcChain(...)->Marginals, runNcmcLoop(...)->NcmcRun, NcmcRun struct
```

Header-only; not matched by `file(GLOB tests/Test*.cpp)` (`CMakeLists.txt`),
so **no CMake change**.

## Constraints

- Pure test motion. `runNcmcLoop` moves verbatim (modulo the bridge-template
  parameter); `runHmcChain` reproduces the ensemble loop call-for-call. Same
  seeds, same move/step counts, same stride, same observable - the sample stream
  is bit-identical. A converted test whose statistic shifts is not a pure move.
- The harness owns no engine state and no RNG; it drives the caller's `HmcDriver`
  and bridge, which already carry the seeded streams (`HmcDriver.hpp:76`). It
  SHALL NOT reseed.
- No test includes a `.cpp`. Header held to production discipline: `#pragma
  once`, IWYU-clean, `@file` block, Doxygen on `runHmcChain`/`runNcmcLoop`/
  `Marginals`/`NcmcRun` stating the loop contract and the telescope identity
  `runNcmcLoop` preserves (`TestNCMCWork.cpp:398-404`).
- One commit for the harness + conversions; formatting separate.

## Predicted breakage

None beyond the intended conversions. Source is frozen; the harness reproduces
both loops exactly. `TestNcmcExplicitSolvent.cpp` has 3 permanently-deferred
`GTEST_SKIP` stubs and 6/8 skip at the default tier
(`TESTS.md section 2`, `D-T2`) - converting its one live construction path to the shared
`runNcmcLoop` SHALL not change which tests skip. Any assertion-count delta means
the loop was not reproduced faithfully - revert.

## Exit criteria (machine-checked; see VERIFY.md section 5)

- Build clean; **pass set identical to B1**; **assertion counts identical to B2**
  in both modes; three consecutive full-suite runs identical.
- **Coverage map (B5) not reduced.**
- `NcmcRun`/`runNcmcLoop` exist once (in the harness); the local copies in
  `TestNCMCWork.cpp` and the mirrored machinery in `TestNcmcExplicitSolvent.cpp`
  are gone.
- The inline `drv.move()` loops in the five ensemble tests are replaced by
  `runHmcChain`; each test body retains only its assertion.
- **No test includes a `.cpp`.**
- `tests/support/SamplingHarness.hpp` <= 300 LOC.
