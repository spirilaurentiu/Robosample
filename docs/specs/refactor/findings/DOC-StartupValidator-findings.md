# DOC-StartupValidator findings (2026-07-13)

## Verified criterion + failure action (from StartupValidator.cpp + initialize caller)
- `Context::checkStartupGeometry` (called once inside `initialize`, before any
  sampling; never mutates sampled state) flags three conditions:
  1. any non-finite coordinate;
  2. hard steric clashes: a non-excluded, non-virtual atom pair closer than
     `kClashNm = 0.08 nm` (minimum-image distance when the box is periodic);
  3. non-finite or > 1e4 kJ/mol initial potential energy.
  Excluded from the clash scan: 1-2/1-3/1-4 bonded/exclusion/scaling14 pairs and
  massless virtual sites.
- Failure action: throws `std::runtime_error` with a descriptive message, UNLESS
  env `ROBO_ALLOW_BAD_START` is set to a non-empty non-"0" value, in which case
  the same message is logged as a warning and the run continues. Documented both.
- The O(N^2) cost is mechanism, not contract (ticket section 2 agreed); not documented as
  a guarantee.

## Coverage gap (recorded)
- No dedicated test exercises `checkStartupGeometry`. tests/TestStability.cpp is
  adjacent evidence for the clash criterion but does not drive this function. The
  criterion here was recovered from source and the `initialize` call site only.
