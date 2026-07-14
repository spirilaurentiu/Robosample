# SPLIT-C5: extract the DrivenRexDriver from Context.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/Context.cpp` (1335 LOC) contains the driven (nonequilibrium-work) replica
exchange round-loop for `RUN_TYPE::RENE`/`REBASONTOP`: the equilibrium segment
plus pairing, the per-replica BAT-scaling drive, the main WTerm swap attempt, and
REBASONTOP's interleaved REMC sub-rounds. This ticket extracts that round-loop as
a cohesive translation unit, downstream of **SPLIT-C4** (which extracted the
label-swap driver and the acceptance algebra `attemptREXSwap`, whose RENE/
REBASONTOP/RENEMC branches this code drives).

Per ARCHITECTURE.md section 1, replica exchange is the intended center, so this is a
core REX module - but this specific code carries a status the ticket SHALL
preserve. Per **OQ-5**, the driven-REX code (`runDrivenRound`, `driveReplica`,
`runInterleavedRemcSubround`) is **documented as reviewed-on-paper but never
compiled or run** (coordinator directive, 2026-07-12 "drop compiling and running
entirely"; `Context.hpp:104-109`, `ReplicaExchange.hpp:10-15`, `Context.cpp:818
-822`). This ticket moves it **unchanged, behind that same status**: it is a pure
relocation, adds no test, compiles no new behavior, and the "reviewed-on-paper,
not build-confirmed" annotations move verbatim with the code.

The `RUN_TYPE::RENEMC` round-loop is not extracted here because it does not exist
yet - `RunREX` throws `std::logic_error` for RENEMC (`src/Context.cpp:1060-1066`),
a Stage 2c TODO; that throw stays in `RunREX` (extracted in C4). Only RENEMC's
*acceptance formula* exists, and it lives in `attemptREXSwap` (C4).

Applicable invariants, restated in full (as referenced by the moved code):

- **INV-3 Stateless Worlds.** A World carries no persistent per-replica state
  between sweeps; per-atom Ground coordinates (nm) are the sole inter-world
  currency. (`runDrivenRound`'s equilibrium sweep pushes/pulls coordinates through
  `replicas_[r].atomsLocations`; the accept decision swaps **labels** via
  `attemptREXSwap`, never coordinates.)
- **INV-5 Mass metric.** HMC momentum seeding and kinetic energy use the same
  mass-matrix operators; `calcLogDetM` is the Fixman kinetic term. (`driveReplica`
  resets `WORK`/`WORK_Jacobian` and reseeds the trial endpoint from the committed
  configuration once per driven range, `src/Context.cpp:890-892`.)
- **INV-7 BAT map/Jacobian agreement.** `applyBatScaling` and its Cartesian
  log-Jacobian read identical `(r,theta)` geometry; one anchor snapshot per round is
  shared by both swap partners so the paired map is an exact involution.
- **INV-9 (BAT-spec, the frozen-anchor involution referenced throughout this
  code; distinct from ARCHITECTURE.md section 5's INV-9 "quaternion double cover").**
  The Context-owned running-mean BAT anchor is state-independent by construction;
  `runDrivenRound` takes exactly **one** frozen snapshot per round
  (`src/Context.cpp:1031`, `batAnchorSnapshot()`) and passes the **same** snapshot
  to every drive of both partners, so the paired scaling map is an exact
  involution and `attemptREXSwap`'s `logCorrectionTerm == 0` holds.
- **Spec INV-10 (drive/run-type pairing; not in ARCHITECTURE.md section 5's INV-1...9
  list).** RENE/REBASONTOP SHALL drive with a volume-changing BAT-scaling world
  (`DistortOption::ScaleBendStretch`); the guard is enforced by
  `checkInv7AndInv10Guards` (extracted in C4) before any round runs
  (`src/Context.cpp:1070-1073`); `driveReplica` relies on it (`src/Context.cpp:899
  -901`).

## Moves (exactly what)
Cohesive translation-unit split (R3 pattern). All three methods stay `Context`
members; `include/Context.hpp` is unchanged. Only **definitions** relocate, and
they move **verbatim** - including every "NOT compiled or run / reviewed-on-paper"
comment.

To `src/workflow/rex/DrivenRexDriver.cpp` [Context member definitions, verbatim]:
- The section banner comment (`src/Context.cpp:818-822`). NOTE: this banner also
  names the INV-7/INV-10 precondition guard, whose function
  (`checkInv7AndInv10Guards`, `823-880`) C4 relocates to `SwapAcceptance.cpp`. The
  banner moves verbatim with the driven round here; splitting it to follow the guard
  would be a comment edit, not pure motion. C4 restates the guard's contract in
  `SwapAcceptance`, so no contract is lost - the split attribution is documented, not
  silent.
- `Context::driveReplica(int,int,double,const robo::BatAnchorStats::Snapshot&)`
  incl. its `std::domain_error` catch -> forced-reject path
  (`src/Context.cpp:882-955`).
- `Context::runInterleavedRemcSubround()` (`src/Context.cpp:957-968`).
- `Context::runDrivenRound(int,bool)` (`src/Context.cpp:970-1057`).

These read `openmmPotential` (stays in residual `Context.cpp`, C4), the shared
`kBoltzmann_kJ` is **not** used here (only `attemptREXSwap`/`runREX` use it, C4).
`runDrivenRound` calls `accumulateBatAnchorStats`/`batAnchorSnapshot` (inline in
`Context.hpp:233-245`, unchanged), `prepareExchangePairs`/`attemptREXSwap`
(extracted in C4, still `Context` members), and `World::applyBatScalingDrive`/
`getDistortJacobianDetLog`/`getDistortOption` (unchanged World API). REX state
members stay declared in `Context.hpp`; the moved definitions access them across
TU. Declarations stay at `include/Context.hpp:339` (`runDrivenRound`), `:355`
(`driveReplica`), `:365` (`runInterleavedRemcSubround`), all private.

`RunREX`'s `if (driven) runDrivenRound(round, verbose);` dispatch
(`src/Context.cpp:1085-1090`) stays inside `RunREX` (moved by C4); it resolves
against the unchanged declaration.

## Public API after this ticket
No new exported symbol; `include/Context.hpp` byte-unchanged. All three methods
are private, so the public-symbol nm set is untouched. A thin
`src/workflow/rex/DrivenRexDriver.hpp` MAY be omitted; the `.cpp` compiles against
`#include "Context.hpp"`, `#include "ReplicaExchange.hpp"`, and carries the
includes the moved bodies use - `<algorithm>` (`std::copy`), `<cmath>`
(`std::sqrt`), `<cstdio>`, `<limits>` (`std::numeric_limits<double>::infinity`),
`<stdexcept>` (`std::domain_error`), `<vector>`. `#pragma once` if a header is
added.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames. The RENE/REBASONTOP
  code moves **unchanged and behind its existing "reviewed-on-paper, uncompiled"
  status (OQ-5)** - this ticket neither wires it into a test nor compiles new
  behavior; the status comments travel with the code.
- `src/Context.cpp` loses lines 818-822, 882-955, 957-968, 970-1057; nothing else
  changes. Includes now only used by the moved code (`<limits>` if unused
  elsewhere) migrate to the new TU.
- `Context`'s class interface is unchanged.
- Invariants that SHALL remain true: INV-3, INV-5, INV-7, the BAT-spec INV-9
  (frozen-anchor involution), and spec INV-10 (all restated above).
- Build: compile `src/workflow/rex/DrivenRexDriver.cpp` into the object library;
  the non-recursive `src/*.cpp` glob (`CMakeLists.txt`) SHALL reach it (flat
  placement or `GLOB_RECURSE`). NOTE: since this code was previously "not
  compiled," adding its TU to the build is the first time the compiler sees it as
  part of the library; a compile error here is a **pre-existing latent defect
  surfaced by the move**, not a regression introduced by it - report it to the
  Architect (OQ-5) rather than "fixing forward," and do not let it force a logic
  edit inside a pure-motion ticket. One commit for the move; formatting separate.

## Predicted test breakage (Phase-A include fixes only)
- None. `runDrivenRound`/`driveReplica`/`runInterleavedRemcSubround` are private,
  named by no test, and reachable only via `RunREX(RENE/REBASONTOP, ...)`, which no
  test invokes (the REX suite exercises REMC and the acceptance algebra only).
  `tests/TestRexAcceptanceAlgebra.cpp` and `tests/test_rex_swap_acceptance_
  algebra.py` reimplement the WTerm/Jacobian formulas in-file and do not call the
  engine, so they are unaffected. `tests/TestBatAnchorInvolution.cpp` exercises
  `BatAnchorStats` directly, not this driver.

## Exit criteria (machine-checked)
- Full build passes (the new TU now compiles as part of the library - see the
  OQ-5 build note); test suite / baseline outputs identical to Stage-0 (B1/B2/B3).
- nm diff on public symbols vs B0: empty.
- include-cycle script clean vs B4; clang-tidy / clang-format / IWYU clean on
  touched files.
- `DrivenRexDriver.cpp` <= 600 LOC (actual ~190). `Context.cpp` shrinks by ~180
  LOC and, after C1-C5, is a thin orchestrator.
- Comment-stripped before/after diff: only the moved blocks (status comments
  included) at their new location and the added `#include`s.
