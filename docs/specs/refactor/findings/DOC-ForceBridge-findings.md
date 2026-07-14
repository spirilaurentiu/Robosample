# DOC-ForceBridge findings

## ARCHITECTURE 6.1 include state (required record)
`include/ForceBridge.hpp` still includes `OpenMMContext.hpp` **concretely**
(`#include "OpenMMContext.hpp"`, and the class holds
`std::vector<OpenMM::Vec3> posCache_` plus OpenMM-typed forwarders). The
forward-declaration resolution described in ARCHITECTURE 6.1 is not applied to
this header. Documented the interface, not the include structure; the coupling
is recorded here as a finding, not written into the contract. The thin-forwarder
/ forward-decl resolution is a separate (non-doc) ticket.

## Bridge concept (recovered from all instantiations)
There is no named `Bridge` concept; the requirement is enforced by
`RobotIntegrator::verletStep<Bridge>` instantiation. The sole required member is
`evaluate(RobotState&)` (`include/RobotIntegrator.hpp:443`, `bridge.evaluate(s)`
is the only call). Confirmed across the three instantiations:
- `ForceBridge` (OpenMM host + fused CUDA in one class);
- `tests/AnalyticForceBridge.hpp` (analytic oracle; same `evaluate` shape, same
  INV-1/INV-2 reduction);
- `tests/ForceBridge.hpp` (OpenMM-free stub whose only member is an empty
  `evaluate(RobotState&)`), which confirms `evaluate` is the minimal surface.
The concept contract (fill `bodyForceG` under INV-1/INV-2, zero `mobilityForce`,
make PE retrievable) is documented from the union of these, not from one.

## Other
- Header-only: no `ForceBridge.cpp`; every method is in-class inline.
- The wrench reduction is delegated to `reduceAtomForcesToBodies` (ForceReducer),
  not re-implemented, consistent with the dedup.

No `@note Assumed:` entries.
