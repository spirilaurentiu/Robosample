# SPLIT-DEDUP-BATBENDSTRETCH: split applyBatScaling into applyBend/applyStretch

Status: draft (Phase A, source split - pure code motion). Executor: `coder`. Style: `styles/spec.md`. Exit criteria machine-checked (VERIFY section 2).

## Context (why)

`MODULES.md section 2` and `ARCHITECTURE.md section 7` flag `applyBatScaling` in
`src/BatScaling.cpp` as a 112-line god-function doing bend + stretch + Jacobian in
one loop (lines 124-236). This ticket splits the per-body work into two helpers -
`applyBend` (the theta rotation of the rigid subtree, lines 163-207) and
`applyStretch` (the radial shift, lines 210-229) - and routes the input geometry read
through the existing `readBatCoord` helper (`BatScaling.cpp` lines 54-94) so the
forward map and the Cartesian log-Jacobian read one `(r, theta)` definition.

Invariant from `ARCHITECTURE.md section 5`, restated in full:

- **INV-7 BAT map/Jacobian agreement.** `applyBatScaling` and its Cartesian
  log-Jacobian read identical `(r, theta)` geometry; one anchor snapshot per round is
  shared by both swap partners so the paired map is an exact involution (INV-9 in the
  BAT spec). Today `calcBatVolumeLogJac` (the Jacobian) reads `(r, theta)` via
  `readBatCoord`, while `applyBatScaling` (the map) re-derives `r0`/`theta0` inline
  (lines 153, 169-170) with a DIFFERENT floating-point expression. The dedup's intent
  is to collapse that to one read so agreement is structural.

## Moves (exactly what)

New file-local helpers in the `src/BatScaling.cpp` anonymous namespace (they need the
subtree, axis, and anchor state, so they take those by parameter):

- `applyStretch` - the radial-shift block (lines 210-229): `muR` lookup, `r1`
  computation, the physical-domain throwing guard (lines 216-222), `dr`, and the
  subtree shift. Returns whether it displaced a DOF (for `nScaled`).
- `applyBend` - the angle block (lines 163-207): `ez`/`theta0`/`axis` construction,
  `muTheta` lookup, `theta1`, the physical-domain throwing guard (lines 187-193),
  `dtheta`, the `rotateAboutAxis` of the subtree, and the `bondDir` update. Returns
  the possibly-updated `bondDir` and whether it displaced a DOF.

`applyBatScaling` (lines 124-236) becomes the driver: the J0 snapshot
(`calcBatVolumeLogJac`, line 132), the per-body loop with `batScaledDofs`/`bodyZRow`/
z-row reads, the `subtree` collection (line 160), then `applyBend` (if `sel.theta`)
and `applyStretch` (if `sel.r`), then the J1 snapshot and `outLnJac` (lines 233-234).
The loop structure, the throwing domain guards, and the `outNScaled`/`outLnJac`
arithmetic are unchanged.

Route-through-`readBatCoord` - SCOPED to what is bitwise-safe (see the FP note below):
- The STRETCH radial read is bitwise-identical between the two definitions:
  `readBatCoord` computes `r = (atomPos[zI]-atomPos[zJ]).norm()`; `applyBatScaling`
  computes `r0 = (atomPos[zI]-posJ).norm()` with `posJ = atomPos[zJ]`. Same operands,
  same operation. The map's `r0` MAY be sourced from `readBatCoord(...).r` with no
  numerical change. NOTE: the two functions use different domain guards -
  `readBatCoord` returns early on `r > 0.0` (`BatScaling.cpp:72`), `applyBatScaling`
  skips on `r0 > 1e-9` (`:154`). Sourcing only the `r0` value from `readBatCoord` is
  bitwise-safe only while `applyBatScaling` keeps its own `> 1e-9` guard; the guard
  SHALL NOT be replaced by `readBatCoord`'s.
- The BEND angle read is NOT bitwise-identical: `readBatCoord` computes
  `cosTheta = dot(v1,v2)/(r*r2)` (one division of the raw dot by the norm product),
  whereas `applyBatScaling` computes `cosTheta0 = dot(bondDir, ez)` on
  pre-normalized vectors. These agree mathematically but round differently. Routing
  `theta0` through `readBatCoord` would perturb `dtheta`, the rotated positions, and
  hence the deterministic HMC output - a B3 (bitwise) violation. Therefore `applyBend`
  SHALL keep its current inline `cosTheta0 = dot(bondDir, ez)` computation. The
  "route through readBatCoord" directive is satisfied for the stretch coordinate only;
  the bend coordinate stays inline and this divergence is documented in the code so
  the INV-7 agreement claim is precise (the map and Jacobian agree in value, and the
  stretch read is now shared; the bend read remains two expressions of the same
  quantity, as it is today).

## FP note (flag for the executing coder and reviewer)

The one-line summary "route the input read through readBatCoord" cannot be taken
literally for theta without breaking B3. This ticket therefore does a PARTIAL route:
stretch shares `readBatCoord`, bend does not. If a future change wants full sharing
(bend included), it is a deliberate numerical change requiring a re-derivation of the
involution FD gate (`tests/TestBatScalingJacobian.cpp` V5) and a tolerance-relaxed B3
baseline - that is a separate, human-approved ticket, not this behavior-preserving
dedup.

## Public API after this ticket

`include/BatScaling.hpp` is UNCHANGED - `applyBatScaling`, `calcBatVolumeLogJac`,
`countScaledDofsStatic`, `batScaledDofs`, and `BatAnchorStats` keep their signatures.
`applyBend`/`applyStretch` are file-local (anonymous namespace), not exported. No
`robo_bindings...so` symbol delta.

## Constraints

- Behavior-preserving. Zero numerical change: the moved bend/stretch blocks keep their
  exact expressions; only the stretch `r0` read is re-sourced from `readBatCoord`,
  which is bitwise-identical.
- No include-structure change; everything stays in `BatScaling.{hpp,cpp}`.
- Invariant that SHALL remain true: **INV-7** (restated above). The map and Jacobian
  read the same `r`; the bend `theta` remains two byte-for-byte-preserved expressions
  of the same geometric angle (documented), exactly as today.
- The two throwing domain guards (`theta1 in (0, pi)`, `r1 > 0`) move verbatim with
  their blocks; their `throw std::domain_error` messages are unchanged.
- One commit for the split; formatting fixes separate.

## Predicted test breakage (Phase-A include fixes only)

**None.** `applyBatScaling` keeps its signature; `applyBend`/`applyStretch` are
file-local and unnamed by any test. The BAT tests call `applyBatScaling` /
`calcBatVolumeLogJac` through the public header:
- `tests/TestBatScalingJacobian.cpp` (the V5 finite-difference map/Jacobian gate),
- `tests/TestBatAnchorInvolution.cpp` (the paired-involution INV-7 check),
- `tests/TestRexAcceptanceAlgebra.cpp`,
- `tests/test_bat_scaling_drive_engine.py`.
All see identical outputs (stretch read is bitwise-safe; bend read unchanged). No
include fix, no assertion-count change.

## Exit criteria (machine-checked)

- Full build passes; test suite / baseline outputs identical to Stage-0 (B1, B2, B3),
  both modes. `TestBatScalingJacobian.cpp`'s FD gate and `TestBatAnchorInvolution.cpp`
  are the sharpest detectors that the map's `(r, theta)` and its Jacobian still agree
  and the involution still holds bitwise.
- Public-symbol nm diff vs B0: empty (only file-local helpers added).
- include-cycle script clean (no new edge).
- clang-tidy / clang-format / IWYU clean on `BatScaling.cpp`.
- `BatScaling.cpp` <= 600 LOC (shrinks); `applyBatScaling` drops from ~112 lines to a
  short driver.
- No test file changed.
