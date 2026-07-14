# DOC-robo_debug findings

## Contradicted hypothesis (Critical-for-doc, not a code change)
- Ticket section 1 hypothesis: "the scanner is compiled out entirely unless
  `ROBO_DEBUG` is defined; production sampling behavior is identical with it
  absent." CONTRADICTED by the code. `robo_debug.hpp:26` hard-defines
  `#define ROBO_DEBUG 1` (not `#ifndef`-guarded), and
  `src/RobotEngine_dynamics.cpp:18` includes the header unconditionally. So the
  probes and the D/DI non-finite dump at `RobotEngine_dynamics.cpp:151-170` are
  compiled INTO and active in the production build. Documented the actual
  contract: the probes are read-only and write only to std::cout, so they do not
  perturb the sampled distribution; "compiled out in production" is achieved by
  a TU not including the header (e.g. RobotIntegrator, which supplies its own
  no-op `ROBO_CHECK`), not by the guard being off.

## Suspected defects (reported, not fixed)
- `ROBO_DEBUG` and `ROBO_VERBOSE` are hard `#define`d, so a `-DROBO_DEBUG=0` /
  `-DROBO_VERBOSE=1` on the command line clashes with these definitions rather
  than overriding them. The header's own comment claiming `-DROBO_VERBOSE=1`
  works is stale. Severity: Medium (maintainability / accidental always-on
  diagnostics), not a sampling-correctness issue.
- `ROBO_VERBOSE` is defined but read by no probe (dead macro). The verbosity
  levels described in the original comment (1 = per-step summary, 2 = per-body
  transforms) are not implemented in this header.

## Coverage gap (per ticket)
- No dedicated test exercises this module; it is guarded, developer-only. The
  only invoked path in the shipped build is the D/DI dump in
  `realizeArticulatedBodyInertias`. Recorded here as a coverage gap, not a doc
  failure.

## Notes
- No `@note Assumed:` used; all documented behavior is backed by the header
  source and its single include site.
