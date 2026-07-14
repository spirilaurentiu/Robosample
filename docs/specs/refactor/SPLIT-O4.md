# SPLIT-O4: extract the alchemy/decoupling forces and state (AlchemyForceFactory) from OpenMMContext.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/OpenMMContext.cpp` carries the NCMC alchemical decoupling machinery: two
force builders, the enable/lambda control surface, and the alchemy member state.
Unlike the standard builders (SPLIT-O3), these read and write `OpenMMContext`
member state (`alchemyEnabled`, `alchemyAtoms`, `alchemyAtomSet`, `alchemyForce`)
and one (`setAlchemicalLambda`) drives the live `OpenMM::Context`. This ticket
groups the alchemy concern - builders + control + state - into an
`AlchemyForceFactory` so intermolecular decoupling lives in one unit.

Two alchemy paths exist, both driven by the single `lambda_inter` global
parameter:
- vacuum/implicit: `createAlchemyCorrectionForce` - exact linear Axrest scaling
  via a `(lambda_inter-1)*standard` interaction group (lines 918-952);
- explicit-solvent (PME/Ewald/CutoffPeriodic): `createAlchemyDecouplingForces` -
  PME-exact charge scaling on the main `NonbondedForce` plus soft-core Axrest and
  hard intra-A custom forces (lines 1090-1165).

The protocol is guidance only: acceptance uses the full Hamiltonian at
`lambda=1`, where both paths reduce exactly to the unmodified field (the comments
at lines 1096-1101, 1141-1147 state this). That exactness at `lambda=1` is the
contract to preserve.

No INV-1/INV-2 force->wrench convention is touched (these define the potential,
not the per-body reduction). The relevant invariant is behavioral: `enableAlchemy`
SHALL yield an ascending, deduplicated Region-A set (lines 966-972), and every
excluded pair stays intramolecular so the added exclusions are energy-neutral
(the CPU shared-neighbor-list rule, lines 947-950, 1141-1147).

## Moves (exactly what)
Target `bridge/AlchemyForceFactory.{hpp,cpp}`. The factory OWNS the alchemy state
so the builders keep reading it; `OpenMMContext` holds one `AlchemyForceFactory`
member and forwards its public alchemy calls to it.

- Alchemy state members
  (`include/OpenMMContext.hpp` lines 276-279:
  `bool alchemyEnabled`, `std::vector<int> alchemyAtoms`,
  `std::set<int> alchemyAtomSet`, `OpenMM::CustomNonbondedForce* alchemyForce`)
  -> `AlchemyForceFactory` members [private].
- `createAlchemyCorrectionForce` (`src/OpenMMContext.cpp` lines 918-952; decl
  `include/OpenMMContext.hpp:228-229`) -> `AlchemyForceFactory` [reads
  `alchemyAtomSet`; now its own member]
- `createAlchemyDecouplingForces` (lines 1090-1165; decl `:257-259`) ->
  `AlchemyForceFactory` [reads `alchemyAtoms`, `alchemyAtomSet`; mutates the
  passed-in `main` NonbondedForce]
- `enableAlchemy(const std::vector<int>&)` (lines 966-972; decl `:47`) ->
  `AlchemyForceFactory` [writes the state]
- `enableAlchemy(int, int)` (lines 974-980; decl `:49`) -> `AlchemyForceFactory`
- `setAlchemicalLambda(double)` (lines 982-988; decl `:50`) -> see AMBIGUITY: its
  body calls `context->setParameter("lambda_inter", ...)` and reads
  `alchemyEnabled`/`ensureInitialized()`, all `OpenMMContext`-owned. It moves as a
  factory method taking the `OpenMM::Context*` (or stays a thin `OpenMMContext`
  forwarder delegating `alchemyEnabled` + `setParameter` to the factory).

`OpenMMContext::initialize()` alchemy dispatch (`src/OpenMMContext.cpp` lines
564-580) stays in `initialize()` (extracted later under SPLIT-O5) but calls the
factory: `alchemyEnabled` -> `alchemyFactory_.enabled()`;
`createAlchemyDecouplingForces(...)` / `createAlchemyCorrectionForce(...)` ->
factory methods; the `alchemyForce = nullptr` / `alchemyForce = ...` assignments
follow the state into the factory.

## Public API after this ticket
`OpenMMContext`'s public alchemy surface is PRESERVED as forwarders:
`enableAlchemy(const std::vector<int>&)`, `enableAlchemy(int,int)`,
`setAlchemicalLambda(double)` keep their exact signatures and delegate to the
member `AlchemyForceFactory`. `bridge/AlchemyForceFactory.hpp` exports
`class AlchemyForceFactory` with those three control methods, the two build
methods, and `enabled()`. Header self-contained: `#include "OpenMM.h"`,
`#include "TopologyElements.hpp"`, `<set>`, `<vector>`, `<algorithm>`.
`#pragma once`.

AMBIGUITY TO FLAG (correctness-adjacent): this is more than byte-for-byte motion
because the four builders/controls read `OpenMMContext` member state. Two shapes
preserve behavior:
1. Factory owns state; `OpenMMContext` forwards (recommended; public signatures
   unchanged, only *internal* state relocates - no public nm delta).
2. Factory stays stateless free functions taking the Region-A set as an argument;
   `OpenMMContext` keeps the state. Simpler motion but `setAlchemicalLambda`/
   `initialize` still touch `context`.
Shape 1 keeps the alchemy concern cohesive and is assumed below. Either way the
public `enableAlchemy`/`setAlchemicalLambda` symbols SHALL survive with identical
signatures (nm public delta empty). Confirm before implementing.

## Constraints
- Behavior-preserving relocation. Zero logic changes to the force expressions,
  the soft-core `alpha=0.5`, the charge-offset scaling, or the exclusion loops.
  Zero renames of the public control methods.
- `src/OpenMMContext.cpp` gains `#include "AlchemyForceFactory.hpp"`, holds an
  `AlchemyForceFactory` member, and loses the moved bodies + state; the
  `initialize()` dispatch requalifies to the factory.
- Invariants that SHALL remain true: `lambda=1` reproduces the unmodified field
  exactly (both paths); Region A stays ascending/deduplicated; added exclusions
  stay energy-neutral (CPU neighbor-list rule). No INV-1...INV-9 force-reduction
  invariant is touched.
- New files auto-globbed into `robosample_objects` (`CMakeLists.txt`); no
  CMake edit. One commit; formatting separate.

## Predicted test breakage (Phase-A include fixes only)
- `tests/TestAlchemy.cpp` and `tests/TestNcmcExplicitSolvent.cpp` /
  `tests/TestNCMCWork.cpp` / `tests/TestNcmcTeleport.cpp` drive alchemy through
  the public `OpenMMContext`/`ForceBridge` surface (`TestAlchemy.cpp:44,222-233`
  uses `OpenMMContext::get()`, `initialize`, `computePotentialEnergyByGroup`,
  `ForceGroupEnergy` - all preserved). They never name `createAlchemy*Force`
  (private). So under shape 1 there is **no** test include fix: the forwarders
  keep the public calls resolving. If shape 2 is chosen and any test reached a
  now-relocated helper, it would retarget to `AlchemyForceFactory.hpp` - none do
  today, so predicted breakage is empty.

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0. The
  alchemy energy tests (`TestAlchemy`, NCMC suite) pass with byte-identical
  `lambda=1` energies and unchanged force-group breakdown.
- nm diff on public symbols vs B0: empty (public `enableAlchemy`/
  `setAlchemicalLambda` retained as forwarders; alchemy state was private).
- include-cycle script clean vs B4; no new upward layer edge; clang-tidy /
  clang-format / IWYU clean.
- `AlchemyForceFactory.hpp` <= 300 LOC; `.cpp` <= 600 LOC (actual ~150).
  `OpenMMContext.cpp` shrinks by ~120 LOC.
