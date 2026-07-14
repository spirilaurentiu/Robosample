# DOC-NMA findings

## Naming note (not a contradiction)
- The ticket names the public symbol `NMA::RouteBNMA`. The code has
  `robo::RouteBNMA` (a result struct) and the entry point `robo::computeRouteBNMA`
  (free function). Documented the actual symbols; no `NMA` namespace exists.

## Verified mode convention (documented as contract)
- `eigval` is `omega^2` ascending under the mass-weighted u-metric; `modeU[k]`
  is the u-space mode `v_k = N y_k` with `N = sqrt(M^-1)`; `uScaleFactors` is the
  softest non-rigid mode renormalized to unit length. Consumer:
  `src/world/sampler/VelocityDistortion.cpp:90-94` copies `uScaleFactors` into
  `uScaleFactors_`; `src/world/sampler/HmcMove.cpp:317-339` re-normalizes it to a
  unit direction and forms the momentum bias `mu = alpha * uhat`. The consumer
  re-normalizes, so it relies on the DIRECTION and the u-space/mass-weighting
  convention, not the exact magnitude. Documented accordingly.
- `computeRouteBNMA` restores `state` to the input `q0` before returning
  (verified `NMA.hpp:355-356`); documented as a `@post` because callers rely on
  the state being unchanged.

## Test coverage
- Only `tests/TestNMALinearAlgebra.cpp` exercises this module, and only the
  internal `jacobiEigh` (ascending order) against `jacobiSymEig`; it explicitly
  states `computeRouteBNMA` needs the OpenMM force singleton and is out of scope.
  The end-to-end NMA path is not unit-tested. Coverage gap noted; not a doc
  failure.

## Notes
- No `@note Assumed:` used.
