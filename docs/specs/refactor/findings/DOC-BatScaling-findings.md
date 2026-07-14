# DOC-BatScaling findings

## Invariant numbering note
- ARCHITECTURE calls the BatScaling map/Jacobian-agreement invariant INV-7; the
  BAT spec and the in-code comments call the involution property "INV-9". Per the
  ticket ("INV-7 ... the involution property (INV-9 in the BAT spec)"), the docs
  use INV-7 consistently for both the shared-geometry agreement and the
  anchor-shared involution to avoid clashing with robot_math's INV-9 (quaternion
  double cover). Recorded so the two INV-9 usages are not conflated.

## INV-7 (documented as contract)
- Shared geometry: `applyBatScaling` computes the map and `outLnJac` together, so
  the returned log-Jacobian is consistent with the geometry actually produced
  (`outNScaled`). Documented as `@post`.
- Involution: the paired map is an exact involution only when both swap partners
  share the SAME frozen anchor `mu`; documented as a `@warning` on
  `applyBatScaling` and on `BatAnchorStats::snapshot`. Verified against the
  anchor-shared call path in the REX round: `include/Context.hpp`,
  `src/workflow/rex/DrivenRexDriver.cpp`, and the oracle
  `tests/TestBatAnchorInvolution.cpp`.

## Flagged limitation (reported, not fixed)
- A scaled body whose parent is itself Ground-adjacent has no `zK` reference
  atom; its angle DOF is silently skipped (omitted from both `outNScaled` and
  `outLnJac`, self-consistently). Documented as a `@note` limitation on
  `applyBatScaling`. Severity: Medium (a real, valid molecular configuration the
  current z-matrix has no fallback reference for); production configs should
  avoid stacking a scaled joint directly on a root-adjacent body until a fallback
  reference atom is added.

## Notes
- No `@note Assumed:` used; the involution is verified at the REX-round anchor
  call site, not inferred from the scaling function alone (per ticket).
