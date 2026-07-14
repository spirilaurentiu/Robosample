# TEST-003: Consolidate the six ad-hoc force bridges into `HarmonicBridge<Policy>`

Status: draft, 2026-07-12. Phase-B test ticket (source frozen). Executor: `coder`. Style: `styles/spec.md`. Companion to
[`../../architecture/TESTS.md`](../../architecture/TESTS.md) section 4 (the six
duplicated bridges), section 3 (the `ForceBridge.hpp` name-collision stub) and section 5.4 (the
seam this executes), and [`VERIFY.md`](VERIFY.md) section 5. Test-only motion.

## Context (why)

Six test files each reinvent the harmonic analytic bridge - all are
lambda/scaled variants of `tests/AnalyticForceBridge.hpp`'s per-atom isotropic
harmonic anchor with the **identical** per-atom->per-body reduction (clear
`bodyForceG`/`mobilityForce`, skip `atomMass==0`, moment about the body origin,
`AnalyticForceBridge.hpp:95-116`). `TESTS.md section 4` names them; located exactly:

| Subclass | File:line | Variation over AnalyticForceBridge |
|---|---|---|
| `ZeroBridge` | `TestBiasForces.cpp:144` | k = 0 (clears forces, applies none) |
| `NcmcLambdaBridge` | `TestNcmcExplicitSolvent.cpp:160` | two anchors, `setLambda` scales the inter term |
| `TwoRobotBridge` | `TestTwoRobotContact.cpp:214` | single anchor + caches raw force into `atomForceG()` when `wantsAtomForces()` |
| `InfForceBridge` | `TestIntegrator.cpp:708` | poisons one body force with `inf` (reject-path probe) |
| `LambdaWellBridge` | `TestNcmcTeleport.cpp:53` | single anchor, `setLambda` scales the whole well |
| `LambdaAnalyticBridge` | `TestNCMCWork.cpp:205` | two anchors (intra + lambda*inter), `setLambda` |

Separately, `tests/ForceBridge.hpp` is a 10-line no-op stub
(`ForceBridge.hpp:7-10`, `evaluate(RobotState&){}`) whose class name **collides**
with the real production `ForceBridge` and sits one letter from
`AnalyticForceBridge` - `TESTS.md section 3` flags it as near-dead, `D-T3` asks when to
retire it. This ticket answers `D-T3`: fold it in as the explicit "no force"
policy so the collision disappears.

This ticket introduces one `tests/support/HarmonicBridge<Policy>` base over the
`AnalyticForceBridge` contract and expresses each of the six (plus the stub) as a
policy.

## Scope

- **In:** one new `tests/support/HarmonicBridge.hpp`; conversion of the six
  subclass definitions and their call sites; deletion of `tests/ForceBridge.hpp`
  after its single use is retargeted.
- **Out:** any engine source; the real `include/.../ForceBridge.hpp` production
  header; the sampling loops that *drive* these bridges (TEST-004); the
  `RobotFixture` prologue (TEST-002).

## Moves (exactly what)

### New: `tests/support/HarmonicBridge.hpp`
- `template <class Policy> class HarmonicBridge` over the `AnalyticForceBridge`
  interface: constructor captures per-atom anchors from a position-realized start
  state `s0` (as every subclass already does, e.g.
  `TestNCMCWork.cpp:207-218`), and `evaluate(RobotState&)` runs the **single**
  shared reduction (`AnalyticForceBridge.hpp:95-116`), delegating the per-atom
  force to `Policy::atomForce(...)` and the energy to `Policy::energy(...)`.
- Policies expressing the observed variants, each SHALL reproduce the exact force
  and energy expressions of the subclass it replaces:
  - `NoForcePolicy` - k = 0; replaces `ZeroBridge` **and** the retired
    `tests/ForceBridge.hpp` stub (the explicit "no force" policy).
  - `SingleWellPolicy` - one anchor, stiffness `k`; replaces the
    `AnalyticForceBridge` usage where a plain well suffices and `TwoRobotBridge`'s
    force term.
  - `LambdaWellPolicy` - one anchor, `setLambda` scales the whole well; replaces
    `LambdaWellBridge` (`TestNcmcTeleport.cpp:53`).
  - `IntraLambdaInterPolicy` - intra anchor + lambda*inter anchor; replaces
    `LambdaAnalyticBridge` (`TestNCMCWork.cpp:205`) and `NcmcLambdaBridge`
    (`TestNcmcExplicitSolvent.cpp:160`) - they are the same two-anchor form
    (`intraAnchor_ = p`, `interAnchor_ = p + (0.05,-0.03,0.04)`).
  - `PoisonForcePolicy` - writes `inf` into one body force; replaces
    `InfForceBridge` (`TestIntegrator.cpp:708`).
- The `TwoRobotBridge` atom-force caching into `atomForceG()` when
  `s.wantsAtomForces()` (`TestTwoRobotContact.cpp:250-260`) is a reduction-side
  behavior, not a per-atom-force variation. Preserve it as a
  `HarmonicBridge` construction flag (`cacheAtomForces`), off by default, on for
  the two-robot suite. This keeps `evaluate`'s reduction the single source of
  truth while reproducing the Cartesian-solvent cache contract exactly.

### Retire: `tests/ForceBridge.hpp`
- Its lone consumer (the kinematics/dynamics TU that compiles OpenMM-free - the
  file's banner names `verletStep` as the never-exercised caller) retargets to
  `HarmonicBridge<NoForcePolicy>`. Delete `tests/ForceBridge.hpp`. The
  name-collision with production `ForceBridge` and `AnalyticForceBridge` is the
  motivation (`TESTS.md section 3`, `D-T3`).

## Target structure (after this ticket)

```
tests/support/
  HarmonicBridge.hpp   # HarmonicBridge<Policy> + the six policies (incl. NoForce)
```

Header-only; not matched by `file(GLOB tests/Test*.cpp)` (`CMakeLists.txt`),
so **no CMake change**. `tests/ForceBridge.hpp` is deleted; confirm no
`CMakeLists.txt` reference names it (it is not in `robosample_add_test`
EXTRA_SOURCES).

## Constraints

- Pure test motion. Each policy's `atomForce`/`energy` SHALL be **algebraically
  identical** to the subclass expression it replaces, including sign, stiffness,
  and anchor offset. The shared `evaluate` reduction SHALL match
  `AnalyticForceBridge.hpp:95-116` term-for-term. A converted test that draws a
  different force is not a pure move - revert.
- Retiring `ForceBridge.hpp` is a deletion of a near-dead stub with a single
  consumer, retargeted in the same commit; it is not a behavioral change to any
  assertion. If any test's pass/assertion count moves, stop - the stub was still
  relied upon and its removal needs a human-gated triage note.
- No test includes a `.cpp`. Header held to production discipline: `#pragma
  once`, IWYU-clean, `@file` block, Doxygen on each policy stating the potential
  it encodes.
- One commit for the consolidation + retarget + stub deletion; formatting
  separate.

## Predicted breakage

- The single TU that included `tests/ForceBridge.hpp` fails to compile until its
  include is switched to `support/HarmonicBridge.hpp` and its `ForceBridge` type
  to `HarmonicBridge<NoForcePolicy>`. This is the one predicted edit; it is
  behavior-preserving (the stub's `evaluate` was a no-op and `NoForcePolicy`
  clears-and-applies-nothing identically).
- No other breakage: source is frozen and each converted bridge is
  expression-identical.

## Exit criteria (machine-checked; see VERIFY.md section 5)

- Build clean; **pass set identical to B1**; **assertion counts identical to B2**
  in both modes; three consecutive full-suite runs identical.
- **Coverage map (B5) not reduced.**
- The six subclass definitions (`ZeroBridge`, `NcmcLambdaBridge`,
  `TwoRobotBridge`, `InfForceBridge`, `LambdaWellBridge`, `LambdaAnalyticBridge`)
  are gone; each call site names a `HarmonicBridge<...>`.
- `tests/ForceBridge.hpp` no longer exists; no source or CMake references it; the
  `ForceBridge` name no longer collides inside `tests/`.
- **No test includes a `.cpp`.**
- `tests/support/HarmonicBridge.hpp` <= 300 LOC (claim an exception with reason if
  the six policies push it over).
