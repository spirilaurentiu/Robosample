# SPLIT-W7: extract the CartesianSolvent sub-integration helpers from World.cpp

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)
`src/World.cpp` (2619 LOC) is a god-object mixing ~12 responsibilities
(`ARCHITECTURE.md section 7`). This ticket extracts the **Cartesian-solvent** concern
(`docs/specs/ncmc_solvent_relax.md`, `docs/specs/two-robot-contact/`): marking a
set of atoms to be advanced in flat Cartesian space by velocity-Verlet inside the
proposal, plus the Maxwell-Boltzmann velocity draw and the flat-metric kinetic
energy for those atoms. The atoms stay welded as 0-DOF rigid bodies (no
Fixman/Jacobian contribution); only their per-atom `(x,v)` move. An empty set is
the welded engine, bit-for-bit.

`setCartesianSolvent` also carries PRECONDITION P1: a flagged atom SHALL sit in a
0-DOF Weld/Rigid body whose atom membership the mask covers exactly, else KE is
double-counted and momentum double-drawn - a silent Boltzmann corruption. Its
guard is the one correctness-critical branch here. Sequenced after
W1/W3/W5/W2/W6 (`README.md section 4`). Invariants below come from `ARCHITECTURE.md` and
are restated in full.

## Moves (exactly what)
- `void World::setCartesianSolvent(const std::vector<int>&)` (`src/World.cpp`
  lines 1833-1903) -> `src/world/sampler/CartesianSolvent.cpp` [public World
  member definition relocated; declaration stays `World.hpp:259`].
- `void World::drawSolventVelocities()` (lines 1905-1917) ->
  `src/world/sampler/CartesianSolvent.cpp` [private World member; decl
  `World.hpp:618`].
- `double World::calcSolventKE() const` (lines 1919-1931) ->
  `src/world/sampler/CartesianSolvent.cpp` [private World member; decl
  `World.hpp:619`].

No file-local anonymous helpers are needed by this trio. The Cartesian-solvent
runtime state lives on `RobotState` (`state_.setCartSolvent` / `cartSolventAtoms`
/ `cartSolventInvMass`, `RobotState.hpp`) - unchanged; this ticket does not touch
`RobotState`.

## Public API after this ticket
No new public header. `World::setCartesianSolvent` keeps its existing public
declaration in `World.hpp`. `src/world/sampler/CartesianSolvent.cpp` is a
`World`-class TU (`#include "World.hpp"`) defining the three relocated members.
IWYU: `<cmath>` (`std::sqrt`, line 1913), `<cstddef>` (`std::size_t`), `<cstdio>`
(`std::fprintf`, line 1897), `<stdexcept>` (`std::runtime_error`, lines 1871,
1884), `<string>` (`std::to_string`), `<vector>`, `World.hpp`.

## Constraints
- Pure code motion. Zero logic changes. Zero symbol renames.
- `src/World.cpp` loses lines 1833-1931. Nothing else changes.
- Files land at `MODULES.md section 1` path (`src/world/sampler/CartesianSolvent.cpp`);
  recursive `src/` glob per SPLIT-W1.
- Invariants that SHALL remain true:
  - **INV-5 Mass metric.** HMC momentum seeding and kinetic energy use the same
    mass-matrix operators; the kinetic term SHALL be assembled from the same masses
    it was drawn against. `drawSolventVelocities` samples `v_s ~ N(0, RT/m_s)`
    (sigma = sqrt(RT*invM), line 1913) and `calcSolventKE` returns `1/2 Sigma m_s |v_s|^2` (line
    1930) over the **same** `cartSolventAtoms`/`cartSolventInvMass` arrays - the
    flat-metric analogue of the generalized draw. Preserve both verbatim so the
    solvent KE the acceptance uses matches the marginal it was drawn from.
  - **INV-2 Virtual sites.** Virtual-site forces are already projected onto parent
    atoms by OpenMM; the reduction skips massless particles. `setCartesianSolvent`
    keeps only real (mass>0) atoms (lines 1840-1844) so massless virtual sites are
    never Verlet-integrated here; the relocated mask filter SHALL stay exact
    because it also drives the skip in `fillAtomPositionsFromBodies`.
  - **INV-3 Stateless Worlds.** Per-atom Ground coordinates are the sole
    inter-world currency; the Cartesian-solvent `(x,v)` remain part of that
    per-atom Cartesian state. PRECONDITION P1's throwing guard (a DOF>0 body, or
    a partially-covered 0-DOF body, throws - lines 1858-1894) SHALL move verbatim:
    it is the check that prevents a silent double-counted Boltzmann density, not a
    cosmetic assertion.
- One commit for the move; formatting fixes in a separate commit.

## Predicted test breakage (Phase-A include fixes only)
- None. No test `#include`s `World.cpp`; no `friend`s. `TestTwoRobotContact.cpp`
  (`#include "World.hpp"`) exercises exactly the `setCartesianSolvent` **public**
  guard (INV-WELD: `world.setCartesianSolvent({2,3})` no-throw / throw cases,
  lines 344, 362, 373) through the unchanged public declaration; it explicitly
  notes `drawSolventVelocities`/`calcSolventKE` are private and re-derives their
  formulas via the analytic bridge instead of calling them
  (`TestTwoRobotContact.cpp:15-24`). It keeps `#include "World.hpp"` verbatim. No
  include-path fix is required by any test.

## Exit criteria (machine-checked)
- Full build passes; test suite / baseline outputs identical to Stage-0.
- `nm` diff on public symbols: empty (`setCartesianSolvent` keeps its mangled
  name; the two draw/KE helpers are private).
- include-cycle script: clean; clang-tidy / clang-format: clean; IWYU clean.
- `CartesianSolvent.cpp` <= 600 LOC (~100 expected).
- Comment-stripped before/after diff of `World.cpp` + `CartesianSolvent.cpp` is
  empty apart from the moved blocks at their new location.
