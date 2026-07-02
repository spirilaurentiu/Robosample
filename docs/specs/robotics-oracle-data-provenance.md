# Robotics-oracle data provenance

Companion to `docs/specs/robotics-oracle-differential.md`. This document traces every
reference value in the oracle fixture back to the exact location in the original
(vendored Simbody) implementation that produces it, and to the generator line that reads
and serializes it. Purpose: auditability — a reviewer can confirm each datum comes from
the intended Simbody cache, not a re-derivation, without re-reading the whole generator.

Line numbers are as of Phase 1b (schema_version = 1). NOTE: they SHALL be re-verified
after any edit to the generator or the Simbody patch (a regeneration that shifts lines
does not invalidate the mapping, only the numbers).

**Phase 2 (storage rewrite, §4/§6 below) note:** the §1–3/5/7 `GEN:` line numbers below
are as of Phase 1b and are now off by a constant +5 (four `#include`s — `cnpy.hpp`,
`sha256.hpp`, `<array>`, `<cstring>`, `<sstream>` minus the pre-existing blank line —
were added at the top of `gen_robotics_oracle.cpp` for the npz/manifest rewrite, §4);
`fillState`/`fillBodyOutput`/`fillMultiState`/`fillAggregateState` and everything they
call are otherwise byte-identical to Phase 1b, so the READ side (§1–3/5/7's actual
subject) is unaffected in substance, only in line offset.

Files referenced:
- Generator: `Robosample/tools/gen_robotics_oracle.cpp` (`GEN`)
- Patch decl: `Robosample/Simbody01/Simbody/include/simbody/internal/SimbodyMatterSubsystem.h` (`SMS.h`)
- Patch impl: `Robosample/Simbody01/Simbody/src/SimbodyMatterSubsystem.cpp` (`SMS.cpp`)
- Patch rep:  `Robosample/Simbody01/Simbody/src/SimbodyMatterSubsystemRep.h` (`Rep.h`)
- Node dump:  `Robosample/Simbody01/Simbody/src/RigidBodyNodeSpec.h` (`RBNS.h`), `RigidBodyNode.h` (`RBN.h`)
- Cache defs: `Robosample/Simbody01/Simbody/src/SimbodyTreeState.h` (`STS.h`)
- Schema:     `tests/fixtures/robotics_oracle/RoboticsOracleTypes.hpp` (`TYPES`)

## 1. Extraction path per quantity

Each row: the quantity, the Simbody accessor used, where that accessor reads the value
(public API or the `[robotics-oracle patch]` `…ForOracle` surface, and the underlying
realized cache field), the generator line that reads it, and the generator line that
serializes it. `[P]` marks a patched (test-fixture-only) accessor; all others are stock
public Simbody API.

| Quantity | Simbody accessor | Reads from (cache field) | GEN read | GEN write |
|---|---|---|---|---|
| `X_GB` | `MobilizedBody::getBodyTransform` | — (public) | GEN:298 | GEN:299–300 |
| `X_FM` | `MobilizedBody::getMobilizerTransform` | — (public) | GEN:301 | GEN:302–303 |
| `V_GB` | `MobilizedBody::getBodyVelocity` | — (public) | GEN:306 | GEN:307–308 |
| `qdot` | `State::getQDot` + `getFirstQIndex` | — (public) | GEN:309 | GEN:310–312 |
| `P` | `SimbodyMatterSubsystem::getArticulatedBodyInertia` | `SBArticulatedBodyInertiaCache::articulatedBodyInertia[mbx]` (STS.h:870) | GEN:316 | GEN:318 |
| `PPlus` | `getArticulatedBodyInertiaPPlusForOracle` `[P]` (SMS.h:2845 / SMS.cpp:2304) | `SBArticulatedBodyInertiaCache::pPlus[mbx]` (STS.h:871) | GEN:317 | GEN:319 |
| `D` | `getHingeInertiaTermsForOracle` `[P]` (SMS.h:2848 / SMS.cpp:2308 → Rep.h:427 → RBNS.h:835 `dumpHingeInertiaForOracle`, `getD` RBNS.h:784) | `SBArticulatedBodyInertiaCache::storageForD` (STS.h:873) | GEN:327 | (→ `minEigD`) GEN:335 |
| `DI` | same call `[P]` (`getDI` RBNS.h:791) | `SBArticulatedBodyInertiaCache::storageForDI` (STS.h:874) | GEN:327 | GEN:330 |
| `G` | same call `[P]` (`getG` RBNS.h:798) | `SBArticulatedBodyInertiaCache::storageForG[2*uIndex]` (STS.h:875) | GEN:327 | GEN:332–333 |
| `minEigD` | local `minEigSym` on `D` (GEN:97–154 cyclic Jacobi) | derived from `D` above | GEN:335 | GEN:335 |
| `Z` | `getArticulatedBodyForceZForOracle` `[P]` (SMS.h:2860 / SMS.cpp:2315) | `SBTreeAccelerationCache::z[mbx]` (STS.h:1158) | GEN:338 | GEN:340–341 |
| `zPlus` | `getArticulatedBodyForceZPlusForOracle` `[P]` (SMS.h:2865 / SMS.cpp:2320) | `SBTreeAccelerationCache::zPlus[mbx]` (STS.h:1159) | GEN:339 | GEN:342–343 |
| `eps` | `getHingeInertiaTermsForOracle` `[P]` (RBNS.h:857, valid once `Stage::Acceleration`) | `SBTreeAccelerationCache::epsilon` (STS.h:1157) | GEN:327/345 | GEN:345 |
| `udot` | `State::getUDot` + `getFirstUIndex` | — (public) | GEN:347 | GEN:349 |
| `A_GB` | `MobilizedBody::getBodyAcceleration` | — (public) | GEN:351 | GEN:352–353 |
| `Mdense` | `SimbodyMatterSubsystem::calcM` | — (public) | GEN:357 | GEN:360 |
| `logDetM` | `SimbodyMatterSubsystem::calcDetM` (returns `ln|M|`) | — (public) | GEN:370 | GEN:371 |
| `reactionBo` | `MobilizedBody::findMobilizerReactionOnBodyAtOriginInGround` | — (public) | GEN:373 | GEN:375–376 |
| `reactionMo` | `MobilizedBody::findMobilizerReactionOnBodyAtMInGround` | — (public) | GEN:374 | GEN:377–378 |

The single-body read/serialize loop is `fillState` (GEN:249–379). Multi-body reuses it
per body via `fillBodyOutput` (GEN:911–974) / `fillMultiState` (GEN:986–1055). Stress
cases store only aggregates via `fillAggregateState` (GEN:1435–1495).

## 2. Patch surface (the only test-only code paths)

The differential's non-public reads go through five `[robotics-oracle patch]` symbols,
all banner-marked and citing §4.3 of the spec:

- `getArticulatedBodyInertiaPPlusForOracle` — decl SMS.h:2845, impl SMS.cpp:2304.
- `getHingeInertiaTermsForOracle` (D, DI, G, eps) — decl SMS.h:2848, impl SMS.cpp:2308 → Rep.h:427 → `RigidBodyNode::dumpHingeInertiaForOracle` (virtual, default no-op RBN.h:422) overridden once in the `RigidBodyNodeSpec<dof>` template (RBNS.h:835).
- `getArticulatedBodyForceZForOracle` — decl SMS.h:2860, impl SMS.cpp:2315.
- `getArticulatedBodyForceZPlusForOracle` — decl SMS.h:2865, impl SMS.cpp:2320.
- `calcAccelerationForOracle` — decl SMS.h:2881, impl SMS.cpp:2327. Capture-path fallback anticipated by the spec; NOT used in Phase 1a/1b (the persisted-cache reads above suffice because forces enter through the official realize pipeline, §4).

NOTE: `P`, `Z`, `zPlus`, `eps` are read from the *persisted* realize caches, which are
valid only after the state is realized to the relevant stage through the official
`System::realize()` pipeline. The patch adds no computation — it only exposes existing
cache fields (and, for D/DI/G, unpacks them dof-generically).

## 3. Force injection and its self-check

Applied body forces (`bodyForceG`) enter Simbody as real force elements, not
`Force::Custom` and not the capture path:

- Single body: `Force::ConstantTorque` + `Force::ConstantForce` (GEN:262–263), added at
  system-build time so they flow through `realize()` into `SBTreeAccelerationCache::z/zPlus`.
- Multi body: per-body plus a Ground force for the §8.2 #7 "force on Ground must not leak
  into udot" discriminator (GEN:1004–1008).

The convention SHALL be checked before any `Z`/reaction value is trusted: `checkReadback`
(GEN:386–400) asserts `Z(force state) − Z(force-free state) == −bodyForceG` to 1e-9,
which holds exactly because `Z` is affine in the applied force. The multi-body variant
`checkReadbackMulti` (GEN:1064–1101) restricts this identity to leaf bodies, since at an
interior body `Z` also accumulates children's contributions.

## 4. Serialization

Phase 2 (docs/specs/robotics-oracle-differential.md §9 storage rewrite): output is one
`<CaseName>.npz` + one `<CaseName>.manifest.json` per case under
`tests/fixtures/robotics_oracle/`, written by `main` (GEN:2051–2113) via
`writeCase`/`writeMultiCase`/`writeAggregateCase` (GEN:1905, 1986, 2033). The generated
C++ header (`RoboticsOracleData.hpp`) and its `emit*` helpers are gone — this replaces,
not layers on top of, Phase 1b's serialization (§9 of the differential spec has the
retrospective on why: a 1.1 MB/~22k-line header compiled into every build, produced
unreviewable diffs, and would not survive Scope B or the fuzz batch).

- **`.npz` (reference arrays, cnpy, `tests/third_party/cnpy`).** Every double-valued
  field of `OracleState`/`OracleBodyOutput`/`OracleAggregateState` is written as a named
  array at its TRUE nq/nu/numBodies extent (never `kMax`-padded), via
  `writeStateArrays`/`writeBodyOutputArrays`/`writeAggregateStateArrays`
  (GEN:1857, 1934, 2013). Naming convention: `state{i}_<field>` for case-wide quantities
  (`i` = 0-based state index, matching `manifest.json`'s `states` list order) and
  `state{i}_body{b}_<field>` for the per-body block in multi-body cases (`b` = 0-based
  body index). Precision is full IEEE-754 double (no truncation — the old header's "set
  precision to 17 significant digits" concern (previously GEN:1929) does not apply to a
  binary `.npy` payload); this is exactly the byte-reproducibility the differential spec's
  §9 storage rule requires. `NpzWriter` (`cnpy.hpp`) batches all of a case's arrays in
  memory and flushes once (a deliberate divergence from upstream cnpy's incremental
  append-and-rewrite `npz_save`, documented in `cnpy.hpp`'s banner) — the on-disk `.npz`
  is still a byte-for-byte standard ZIP (`unzip -l`/`numpy.load` both open it unmodified,
  verified against real `numpy.load()` during implementation).
- **`.manifest.json` (reviewable metadata, hand-rolled JSON emission, `writeManifest`
  GEN:1776).** Carries `schema_version`, `case_name`/`case_kind`
  (`"single"|"multi"|"aggregate"`), the model spec (`model.bodies[]`: `parent`, `joint`
  int + `jointName`, `X_PF`/`X_BM`, `mass`, `com_B`, `unitInertia_B` — the §4.1
  correspondence, small and human-reviewable, one entry per body), `nq`/`nu`,
  `reportBody` (aggregate cases only), the ordered `states` label list, the full `arrays`
  list (`{name, shape}` per npz member — this is what the disasm-side no-silent-gap guard
  checks against, `tests/RoboticsOracleLoader.hpp::checkNoSilentGap`), `rng_seed` (`null`
  — no RNG is used by the Phase 1a/1b battery; §8.3's future fuzz batch is required to
  populate this), and `sha256_npz` (SHA-256 over the case's `.npz` FILE BYTES, computed
  by re-opening the just-written file — `tests/third_party/sha256`, a from-scratch
  FIPS-180-4 implementation vendored for the same no-network-access reason as `cnpy`).
  Provenance git SHAs (`git rev-parse HEAD`, captured via `popen`, GEN:1669–1692):
  `disasm_git_sha` (this repo), `generator_git_sha` (`Robosample/Robosample`, the
  "refactor" clone `gen_robotics_oracle.cpp` itself lives in), and
  `simbody_patch_git_sha` (`Robosample/Robosample/Simbody01`, the tree carrying the §4.3
  clone-getter patch) — a superset of the differential spec's literal "generator and
  Simbody-patch SHA" ask, cheap and strictly more useful for provenance. Each SHA has a
  paired `*_dirty` boolean (`git status --porcelain` non-empty); `simbody_patch_git_dirty`
  is EXPECTED true (the §4.3 patch is deliberately uncommitted local state on the vendored
  Simbody01 tree, never shipped).
- **Determinism (byte-reproducibility, verified during implementation).** Two consecutive
  runs of the generator into different output directories produce byte-identical `.npz`
  files (checksummed) and manifests identical except for `regen_command` (which embeds the
  output path). The ZIP writer never embeds a wall-clock timestamp (fixed placeholder
  mod-time/date, `cnpy.cpp`), so no run-to-run entropy leaks into the archive bytes.
- A schema change (new/renamed `OracleState`/`OracleBodyOutput`/`OracleAggregateState`
  field) SHALL update `TYPES`, the matching `write*Arrays` helper here, the loader's
  `load*Case`/`load*Arrays` counterpart (`tests/RoboticsOracleLoader.hpp`), AND the
  no-silent-gap guard's `expected*Fields(schemaVersion)` lists — bump `kSchemaVersion` so
  stale fixtures fail loud instead of silently reading a missing field as zero.

## 5. Determinism (regeneration)

The generator uses no RNG. Every `(q, u, force)` input is a hard-coded constant or a
deterministic formula, so identical generator + Simbody build → bit-identical output.

- Per-case inputs: `buildTorsionCase` GEN:427, `buildFreeCase` GEN:464, `buildSliderCase`
  GEN:551, `buildCylinderCase` GEN:587, `buildCartesianCase` GEN:623, `buildBallCase`
  GEN:663, `buildBallNeedleCase` GEN:714, `buildBendStretchCase` GEN:764,
  `buildSphericalCoordsCase` GEN:806, `buildFreeLineCase` GEN:846, `buildRigidCase`
  GEN:515 (0-DOF).
- Stress inputs are deterministic index formulas: depth chain `q[i]=0.05 sin(0.7 i+1)`,
  `u[i]=0.08 cos(0.5 i+0.3)` (GEN:1556–1569); conditioning GEN:1614–1626.

The §8.3 fuzz batch (future) SHALL record its seed in the fixture metadata so its
"random" states are reproducible under this same no-live-RNG-at-gate-time contract.

## 6. On-disk schema

NOTE: `TYPES` (`RoboticsOracleTypes.hpp`) is now PURELY the in-memory schema both sides
fill at runtime — it no longer describes the on-disk layout (that was true only while
fixtures were a generated C++ header the disasm side `#include`d directly and the struct
literal *was* the file format). The struct is still the single source of truth for field
names/packing conventions (spatial vectors `[angular; linear]`; rotations row-major 3×3;
`ArticulatedInertia` as J (SymMat33 packed xx,xy,yy,xz,yz,zz), F (Mat33 row-major), M
(SymMat33 packed); `DI`/`Mdense` row-major `nu×nu`; `G` as `nu` spatial columns —
TYPES:22–31) and is engine-agnostic (no SimTK/robo dependency), so it still compiles
unchanged on both the generator and the disasm loader. Structs: `OracleState`/`OracleCase`
(single body, TYPES:42–117), `OracleBodyOutput`/`OracleMultiState`/`OracleMultiCase`
(multi body, TYPES:155–220), `OracleAggregateState`/`OracleAggregateCase` (stress,
aggregates only, TYPES:225–258).

The ACTUAL on-disk layout is the `.npz`/`.manifest.json` pair described in §4: the
generator's `write*Arrays` helpers flatten a `TYPES` struct into named npz arrays at
generation time, and the loader's `load*Arrays` helpers (`tests/RoboticsOracleLoader.hpp`)
reassemble the identical struct at test-load time from those same arrays plus the
manifest's `model` block. Neither side ever serializes the struct's raw memory layout —
only its named fields — so a `TYPES` field reordering is schema-inert (only a
rename/removal, guarded by the no-silent-gap check, is a breaking schema change).

## 7. Audit notes (known traps, each already handled)

1. `logDetM` is already in log space (`calcDetM` returns `ln|M|`); the generator does not
   re-log it (GEN:363–367).
2. Quaternion inputs are pre-normalized at generator build time (e.g. GEN:464, 663) so
   both engines consume the same unit quaternion.
3. `BendStretch`/`SphericalCoords` states stay off the `q=0` coordinate singularity
   (GEN:764, 806; §8.2 #1).
4. Frame-dependent stage-3/4 comparisons are trusted only after the frame-equality gate
   (`X_FM` == `getMobilizerTransform`) passes (TYPES:57–59).
5. Element-wise `DI/PPlus/G` diffs are gated on `minEigD` above the lock (§6.1); the
   conditioning-stress case compares aggregates instead.
