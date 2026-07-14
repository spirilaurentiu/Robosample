# TEST-005: Split the two oracle monsters along their case families

Status: draft, 2026-07-12. Phase-B test ticket (source frozen). Executor: `coder`. Style: `styles/spec.md`. Companion to
[`../../architecture/TESTS.md`](../../architecture/TESTS.md) section 2 (FAST\* fixture
volume), section 3 (contract tests - the crown jewels) and section 5.7 (the seam this
executes), `D-T4`, and [`VERIFY.md`](VERIFY.md) section 5. Test-only motion, navigability
only.

## Context (why)

Two files carry the differential robotics oracle against Simbody fixtures - the
crown-jewel contract fixtures that feed the INV-8-adjacent kinematic contracts
(`TESTS.md section 3`):

- `tests/TestRoboticsOracle.cpp` - **1283 lines**, 44 `TEST`s, one per fixture
  case, ABA per case.
- `tests/TestRoboticsOracleMolecule.cpp` - **1171 lines**, builds a full
  `World`/`Context` per molecule fixture.

`TESTS.md section 5.7`/`D-T4` propose splitting them along the case families the loader
**already partitions** so each binary stays ~300 lines and the cases parallelize
under `ctest`. The partition is explicit in
`tests/RoboticsOracleLoader.hpp`: one loader + one runner per family.

### Family boundary (from the loader and the runners)
| Family | Loader entry | Runner (`TestRoboticsOracle.cpp`) | Cases |
|---|---|---|---|
| single-state | `loadCase` (`RoboticsOracleLoader.hpp:356`) | `runOracleCase` (`:225`) | `Torsion`,`Free`,`Rigid`,`Slider`,`Cylinder`,`Cartesian`,`Ball`,`BallNeedle`,`BendStretch`,`SphericalCoords`,`FreeLine` (`:969-1019`) |
| multi-system + body-output | `loadMultiCase` (`:388`) | `runOracleMultiCase` (`:524`) | `MixedChain`,`Forest`,`WideStarHub`,`RigidMidChain`,`DuplicateMolecules`,`ForceDiscriminators` (`:1024-1050`) |
| aggregate | `loadAggregateCase` (`:431`) | `runOracleAggregateCase` (`:750`) | `DepthChain`,`ConditioningStress` (`:1057-1065`) |
| fuzz | `loadFuzzCase` (`:463`) | `runOracleFuzzCase` (`:867`) | `FuzzStates_*`, `FuzzTopo_*` (`:1075-1134+`) |

The body-output `expectedBodyOutputFields` list (`RoboticsOracleLoader.hpp:65`) is
consumed **inside** the multi family (`checkNoSilentGap`'s `kind=="multi"` branch,
`:142-155`), so body-output travels with multi-system, not as a separate binary.
The Molecule file partitions by "symbolic differential" vs "numeric differential"
runners (`runMoleculeNumericDifferential`,
`TestRoboticsOracleMolecule.cpp:543`) and by fixture system (`10ala`, `1APQ`).

## Scope

- **In:** split each monolith into ~300-line `Test*.cpp` binaries, one per case
  family; extract the shared runners/loader glue into a header both halves
  include; register the new binaries in CMake.
- **Out:** any engine source; the fixture data under
  `tests/fixtures/robotics_oracle*`; **any change to what is compared.** The
  staged comparison bodies (`runOracleCase`/`runOracleMultiCase`/
  `runOracleAggregateCase`/`runOracleFuzzCase`) and the loader are the migration's
  hard constraint - they move unchanged (`RoboticsOracleLoader.hpp:11-14`: "this
  file only changes where the struct's data comes from, never what is compared").

## Moves (exactly what)

### New shared header: `tests/support/RoboticsOracleRunners.hpp`
- The four runner functions (`runOracleCase` `:225`, `runOracleMultiCase` `:524`,
  `runOracleAggregateCase` `:750`, `runOracleFuzzCase` `:867`) and the
  `kFixtureDir` binding move here **verbatim** so every split binary shares one
  copy. They are currently file-local to `TestRoboticsOracle.cpp`; promotion to a
  shared header changes their linkage location, not their bodies.

### Split `TestRoboticsOracle.cpp` (1283 -> ~5 files)
- `TestRoboticsOracleSingleState.cpp` - the 11 `runOracleCase` `TEST`s
  (`:969-1019`).
- `TestRoboticsOracleMultiSystem.cpp` - the 6 `runOracleMultiCase` `TEST`s
  (`:1024-1050`), body-output included.
- `TestRoboticsOracleAggregate.cpp` - the 2 `runOracleAggregateCase` `TEST`s
  (`:1057-1065`).
- `TestRoboticsOracleFuzz.cpp` - the `FuzzStates_*` / `FuzzTopo_*` `TEST`s
  (`:1075-1134+`).
- Each `#include "support/RoboticsOracleRunners.hpp"` +
  `RoboticsOracleLoader.hpp`. The `RoboticsOracleLoader.hpp` header comment
  "Included by exactly one TU" (`:16`) SHALL be updated to reflect the new
  inclusion set (comment-only; the `inline` functions already tolerate multiple
  TUs).

### Split `TestRoboticsOracleMolecule.cpp` (1171 -> ~4 files)
- Along the symbolic vs numeric differential runners and the fixture systems:
  e.g. `TestRoboticsOracleMoleculeSymbolic.cpp` (the `RigidWeldRoot`/`RigidFreeRoot`/
  `RegularFlexible`/`Cyclic1APQ` symbolic `TEST`s, `:840-953`) and
  `TestRoboticsOracleMoleculeNumeric.cpp` (the `*Numeric` +
  `runMoleculeNumericDifferential` `TEST`s, `:998-1100`), with the shared molecule
  loader glue in a sibling `tests/support/` header if it too is file-local.

### CMake (`CMakeLists.txt`)
- The oracle wiring `foreach(_oracle_target TestRoboticsOracle
  TestRoboticsOracleMolecule)` matches by **exact target name** and attaches
  cnpy/zlib + the `ROBOTICS_ORACLE*_FIXTURE_DIR` macros. Extend the loop list to
  the new target names (all `TestRoboticsOracle*` binaries need cnpy/zlib and the
  fixture-dir macros). The auto-glob picks up the new `Test*.cpp` as
  executables automatically; only this fixture-wiring list needs the new names.

## Target structure (after this ticket)

```
tests/support/RoboticsOracleRunners.hpp        # runOracleCase/Multi/Aggregate/Fuzz + kFixtureDir
tests/TestRoboticsOracleSingleState.cpp        # single-state family
tests/TestRoboticsOracleMultiSystem.cpp        # multi-system + body-output family
tests/TestRoboticsOracleAggregate.cpp          # aggregate family
tests/TestRoboticsOracleFuzz.cpp               # fuzz family
tests/TestRoboticsOracleMoleculeSymbolic.cpp   # molecule symbolic differential
tests/TestRoboticsOracleMoleculeNumeric.cpp    # molecule numeric differential
```

## Constraints

- **Navigability only. Zero behavior/assertion change.** These are crown-jewel
  contract fixtures; their comparison content SHALL survive untouched. Each moved
  `TEST` body and each moved runner is byte-identical to its origin (modulo the
  file it lives in). The loader (`RoboticsOracleLoader.hpp`) does not change.
- The set of `TEST(RoboticsOracle, *)` and `TEST(RoboticsOracleMolecule, *)`
  cases across the new binaries SHALL equal the pre-split set exactly - same
  suite names, same test names, same count. Splitting is a re-partition of the
  same tests into more files, never an add/drop.
- No test includes a `.cpp`. The runners move to a header, not a `.cpp`.
- Runners in `RoboticsOracleRunners.hpp` stay `inline`/header-only (ODR-safe
  across the new TUs). Header held to production discipline.
- One commit for the split + CMake target-list extension; formatting separate.

## Predicted breakage

- Until the CMake `foreach` target list (`CMakeLists.txt`) names the new
  `TestRoboticsOracle*` binaries, they link without cnpy/zlib and without the
  `ROBOTICS_ORACLE*_FIXTURE_DIR` macros, failing to build/run. Extending the list
  is the one predicted CMake edit.
- No other breakage: source is frozen, comparisons are byte-identical, fixtures
  unchanged.

## Exit criteria (machine-checked; see VERIFY.md section 5)

- Build clean; **pass set identical to B1** (same case names, now spread across
  the new binaries); **assertion counts identical to B2**; three consecutive
  full-suite runs identical.
- **Coverage map (B5) not reduced.**
- Comment-stripped diff of every moved `TEST` body and every moved runner against
  its origin is empty (mechanically proves zero comparison change).
- `TestRoboticsOracle.cpp` and `TestRoboticsOracleMolecule.cpp` no longer exist as
  monoliths; each successor binary <= ~300 LOC (claim an exception with reason for
  any that a family's case count pushes over).
- The `ctest` case list shows the oracle cases distributed across the new
  binaries (they now parallelize under `ctest -j`).
- **No test includes a `.cpp`.**
