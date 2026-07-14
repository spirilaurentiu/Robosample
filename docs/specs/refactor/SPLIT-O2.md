# SPLIT-O2: fold the box math into PeriodicBox; remove the duplicate computePeriodicBoxVectors_Context

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/OpenMMContext.cpp` carries `OpenMMContext::computePeriodicBoxVectors_Context`
(lines 717-753): the triclinic reduced-box-vector construction from cell lengths
and angles. This is a byte-for-byte duplicate of the pure math already hoisted
into `include/PeriodicBox.hpp` as `robo::pbc::reducedBoxVectors` (lines 71-110).
The PeriodicBox header itself documents the duplication: its file comment states
`reducedBoxVectors` is "the body of `OpenMMContext::computePeriodicBoxVectors_Context`,
minus the OpenMM types" (`PeriodicBox.hpp:9`, `:68`). Two copies of the same
lattice-reduction logic is a divergence risk: a fix to one silently skips the
other.

`robo::pbc::reducedBoxVectors` is already the single source the unit tests bind
(`tests/TestPeriodicBoundary.cpp:47,131-132,191,232,248`) and the production
clash scan uses `robo::pbc::minimumImage` (`src/Context.cpp:304-305`). Only the
OpenMMContext method remains un-deduplicated.

This is a **dependency-defect resolution ticket (VERIFY section 3)**, not a plain code
motion: it removes a member function. Grep confirms
`computePeriodicBoxVectors_Context` has **zero production callers** - the only in-
scope references are its declaration (`include/OpenMMContext.hpp:201-207`), its
definition (`src/OpenMMContext.cpp:717-753`), and a comment in
`tests/TestPeriodicBoundary.cpp:15`. (The `Robosample/` nested SimTK copy at
`Robosample/{include,src}/OpenMM.*` is out of scope per ARCHITECTURE section Scope.)

No ARCHITECTURE.md invariant names this method. The convention it SHALL NOT
violate is the reduced lower-triangular box convention already stated in
`PeriodicBox.hpp:15-20`: `a=(ax,0,0)`, `b=(bx,by,0)`, `c=(cx,cy,cz)`; diagonal
positive; off-diagonals reduced. `reducedBoxVectors` already encodes it.

## Moves (exactly what)
Because the method is dead code that duplicates `robo::pbc::reducedBoxVectors`,
the recommended resolution is **removal**, not relocation:

- `OpenMMContext::computePeriodicBoxVectors_Context` definition
  (`src/OpenMMContext.cpp` lines 717-753) -> **deleted**. Its logic already lives
  verbatim in `robo::pbc::reducedBoxVectors` (`include/PeriodicBox.hpp:71-110`);
  no new code is written.
- Its declaration (`include/OpenMMContext.hpp` lines 201-207) -> **deleted**.

`include/PeriodicBox.hpp` is unchanged: it is already the fold target and the
single source of truth. `OpenMMContext.hpp` may drop `<tuple>` (line 8) if no
other member returns a tuple after removal (verify by grep before dropping).

AMBIGUITY TO CONFIRM (see human-approval note below): if the public symbol
`OpenMMContext::computePeriodicBoxVectors_Context` MUST be retained for an
external caller not visible in-tree, the alternative is to keep the method as a
thin forwarder whose body calls `robo::pbc::reducedBoxVectors(...)` and repacks
the result into the returned `std::tuple<OpenMM::Vec3,OpenMM::Vec3,OpenMM::Vec3>`.
That preserves the symbol and still kills the duplicate math. Removal is preferred
because the symbol is provably unused in scope.

## Public API after this ticket
`include/PeriodicBox.hpp` exports (unchanged): `robo::pbc::BoxVectors`,
`robo::pbc::minimumImage`, `robo::pbc::minimumImageDistance`,
`robo::pbc::reducedBoxVectors`, `robo::pbc::orthorhombic`. Already self-contained
(`<array>`, `<cmath>`, `#pragma once`).

`OpenMMContext` loses one public method:
`computePeriodicBoxVectors_Context`. This is the one public-symbol delta in this
ticket.

## Constraints
- Logic-preserving dedup, not blind motion: the deleted body is a proven
  duplicate of `reducedBoxVectors`; nothing new is authored.
- **Public-symbol delta is non-empty** (one method removed). Per VERIFY section 3 this
  ticket is flagged **human-approved** and carries its own nm-delta rationale:
  the removed symbol has zero in-scope callers (grep evidence above), and the
  surviving `robo::pbc::reducedBoxVectors` is behavior-identical. VERIFY criterion
  4 is relaxed only under this flag.
- Before/after include edge removed (checked vs B4): no include edge changes -
  `OpenMMContext.cpp` already includes neither `PeriodicBox.hpp` (it did not use
  it) nor gains it; it simply loses dead code. Confirm `<tuple>` removal does not
  break other decls.
- Invariants that SHALL remain true: the reduced lower-triangular box convention
  (`PeriodicBox.hpp:15-20`); INV-1...INV-9 are untouched (no force/wrench path here).
- One commit for the removal; any `#include` pruning in a separate commit.

## Predicted test breakage (Phase-A include fixes only)
- None functional. `tests/TestPeriodicBoundary.cpp` and `tests/TestAlchemy.cpp`
  already bind `robo::pbc::*` from `PeriodicBox.hpp`
  (`TestPeriodicBoundary.cpp:40,44-47`; `TestAlchemy.cpp:45,126`), never the
  removed OpenMMContext method. `TestPeriodicBoundary.cpp:15` mentions the method
  only inside a comment; a comment is not a compile dependency, so no include fix
  is required. If the reviewer wants the stale comment updated, that is a Nit, not
  a Phase-A include fix, and is out of scope for this pure-dedup commit.

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0 (B1/B2/B3).
  `TestPeriodicBoundary` and `TestAlchemy` pass unchanged.
- nm diff on public symbols vs B0: exactly one removal
  (`OpenMMContext::computePeriodicBoxVectors_Context`), matching the
  human-approved rationale. No other delta.
- include-cycle script clean vs B4; clang-tidy / clang-format / IWYU clean.
- `OpenMMContext.cpp` shrinks by ~37 LOC; `PeriodicBox.hpp` unchanged and <= 300
  LOC (actual ~119).
