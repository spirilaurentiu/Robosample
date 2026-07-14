# TEST-006: Split `TestNCMCWork` into a FAST algebra file and a SLOW work-chain file

Status: draft, 2026-07-12. Phase-B test ticket (source frozen). Executor: `coder`. Style: `styles/spec.md`. Companion to
[`../../architecture/TESTS.md`](../../architecture/TESTS.md) section 2 (the mixed
FAST/SLOW file) and section 5.8 (the seam this executes), and [`VERIFY.md`](VERIFY.md)
section 5. Test-only motion.

## Context (why)

`tests/TestNCMCWork.cpp` (915 lines) is classified **mixed** in `TESTS.md section 2`:
"tests 1-8 algebra (FAST), 9-16 work chains (SLOW)". The two halves share one
translation unit, so the file-granularity CI split (FAST every push, SLOW
nightly) cannot separate them - a single file is either scheduled or skipped
whole. Splitting the file along the algebra/work-chain boundary lets the FAST
half run on every push while the work-chain half moves to the slow tier.

### The boundary (from the `TEST` macros)
- **FAST - pure protocol-schedule algebra** (`ncmc::protocolLambda`, no
  sampling): `NcmcProtocol.*` - `StartsAtOne` (`:79`),
  `BothEndpointsPinnedToOne` (`:90`), `IsPalindrome` (`:105`),
  `RangeIsUnitInterval` (`:118`), `TroughIsMinimumAndReachesZero` (`:133`),
  `MonotoneDownThenUp` (`:146`), `HoldZeroBlockIsCentered` (`:161`),
  `DegenerateBudgetsAreSafe` (`:178`). These call only `ncmc::protocolLambda`
  and assert schedule properties - no `runNcmcLoop`, no engine stepping.
- **SLOW - work-chain / integration identities** (drive `runNcmcLoop`, step the
  Verlet integrator): `NcmcWork.AccumulatorIsTheFixedQTelescope` (`:416`),
  `WorkPlusHeatEqualsDeltaH` (`:466`), `WorkDeltaHGapShrinksUnderStepRefinement`
  (`:495`), `FlatLambdaOneGivesZeroWorkAndPlainHmcDeltaH` (`:535`),
  `NcmcJacobian.FixedQPerturbationLeavesMassMetricAndPitchInvariant` (`:689`),
  `NcmcReversibility.PalindromeMapIsMomentumFlipReversible` (`:751`),
  `NcmcCrooks.SelfReverseProtocolFlipsWorkSign` (`:801`),
  `NcmcFixman.WeldedSolventLeavesMassMetricInvariant` (`:889`). These all
  construct the `LambdaAnalyticBridge` and call `runNcmcLoop`
  (`:433,:478,:511,:769,:774,:818,:823`).

The shared support lives above the `TEST`s: `LambdaAnalyticBridge` (`:205`),
`freeTorsionChain` (`:286`), `seedDerivatives` (`:308`), `kineticEnergy`
(`:322`), `primePositions` (`:342`), `NcmcRun` + `runNcmcLoop` (`:348`,`:361`).

## Scope

- **In:** split `TestNCMCWork.cpp` into two `Test*.cpp` binaries along the
  FAST/SLOW boundary above; route the shared support to the pieces that need it.
- **Out:** any engine source; the `robo::ncmc::protocolLambda` implementation
  under test; the tier-gate mechanism (TEST-001 owns the `SlowStatTest` fixture
  and `testTier()`); the bridge consolidation (TEST-003) and the
  `runNcmcLoop`/`NcmcRun` promotion (TEST-004). If TEST-003/TEST-004 have already
  landed, the SLOW half includes `support/HarmonicBridge.hpp` /
  `support/SamplingHarness.hpp` rather than the local copies; otherwise the local
  support moves with the SLOW half. State the sequencing assumption in the
  execution note.

## Moves (exactly what)

### New: `tests/TestNcmcProtocol.cpp` (FAST)
- The 8 `NcmcProtocol.*` `TEST`s (`:79-190`), moved verbatim. Includes only what
  the schedule algebra needs (`ncmc` header + gtest). No bridge, no engine, no
  `runNcmcLoop`. This binary is pure algebra and stays in the every-push FAST
  set (`TESTS.md section 2` FAST list).

### New: `tests/TestNcmcWorkChain.cpp` (SLOW)
- The 8 work-chain / Jacobian / reversibility / Crooks / Fixman `TEST`s
  (`:416-908`), moved verbatim, plus the shared support they consume
  (`LambdaAnalyticBridge`, `freeTorsionChain`, `seedDerivatives`, `kineticEnergy`,
  `primePositions`, `NcmcRun`, `runNcmcLoop`) - either moved with this file or
  pulled from `tests/support/` if TEST-003/TEST-004 landed first.
- This binary carries the sampling loops; it joins the slow tier
  (`ROBOSAMPLE_SLOW_TESTS` / `testTier()`, per TEST-001). Wiring the gate is
  TEST-001's job; this ticket only ensures the SLOW tests live in a file that
  *can* be gated as a unit.

`TestNCMCWork.cpp` is deleted once both successors exist.

## Target structure (after this ticket)

```
tests/TestNcmcProtocol.cpp    # FAST: NcmcProtocol.* schedule algebra (8 tests)
tests/TestNcmcWorkChain.cpp   # SLOW: work/Jacobian/reversibility/Crooks/Fixman (8 tests) + support
```

Both matched by `file(GLOB tests/Test*.cpp)` (`CMakeLists.txt`), so **no
CMake change** - the auto-glob registers them and drops the deleted monolith on
the next configure (`CONFIGURE_DEPENDS`).

## Constraints

- Pure test motion. Every moved `TEST` body is byte-identical to its origin; the
  shared support functions move unchanged. The set of `Ncmc*` tests across the
  two files SHALL equal the pre-split set exactly - same suite names, same test
  names, same count (16).
- No `TEST` is added, dropped, or renamed. The FAST/SLOW partition follows the
  existing suite boundary (`NcmcProtocol` vs the four work suites); it does not
  reclassify any individual test.
- No test includes a `.cpp`. Any support moved with the SLOW half stays
  header-discipline if promoted to `tests/support/`, otherwise remains
  file-local to `TestNcmcWorkChain.cpp`.
- One commit for the split + monolith deletion; formatting separate.

## Predicted breakage

None beyond the intended split. Source is frozen. The FAST half is pure algebra
and always runs; the SLOW half's gating is applied by TEST-001, not here, so at
this ticket's exit the two files together produce the identical pass set and
assertion counts as the pre-split monolith in both gated modes.

## Exit criteria (machine-checked; see VERIFY.md section 5)

- Build clean; **pass set identical to B1**; **assertion counts identical to B2**
  in both modes; three consecutive full-suite runs identical.
- **Coverage map (B5) not reduced.**
- Comment-stripped diff of every moved `TEST` body against its origin is empty.
- `TestNCMCWork.cpp` no longer exists; the 8 `NcmcProtocol.*` tests are in
  `TestNcmcProtocol.cpp` and the 8 work-chain tests in `TestNcmcWorkChain.cpp`.
- `TestNcmcProtocol.cpp` builds and runs without any engine-stepping symbol
  (verifiable: it links no `runNcmcLoop`/`RobotEngine::stepTo` reference) -
  proving the FAST half is genuinely sampling-free and CI can schedule it every
  push.
- **No test includes a `.cpp`.** Each successor <= 600 LOC.
