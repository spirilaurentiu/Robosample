# DOC-World findings

## Style decision (applies to every World-layer DOC ticket)
The approved DOC-CALIBRATION trio (bridge/ForceReducer.{hpp,cpp},
bridge/GpuKinematics.cpp, tests/TestReactionForces.cpp) documents in the house
`//` contract-block style with inline `verify:` annotations, NOT Javadoc
`/** @brief */` (documenter.md Sec.8). Per Sec.8 "match the approved calibration
exemplar", the exemplar was followed; the exemplar-vs-Sec.8 conflict is recorded
here as required. Doxygen would need `///`/`/**` to index these; the exemplars
show the human gate accepted `//` blocks as the house form. No files were
converted to Javadoc.

## State of the module
World.hpp declarations and the SPLIT-relocated definitions were already
documented to contract quality by the SPLIT tickets. Residual doc-phase edits
were three genuinely bare public declarations, filled in include/World.hpp:
`setAtomsLocationsInGround`, `getAtomsLocationsInGround`, `setTemperature`.

## Verified hypotheses
- INV-3 stateless Worlds: verified. `getAtomsLocationsInGround` /
  `setAtomsLocationsInGround` are the sole inter-world currency; the latter
  reconstructs all frames/q from coordinates alone (GeometryFitter.cpp:89-141).
  Documented at both public declarations.
- Owns RobotModel/RobotState/ForceBridge/ConstraintSet/rng_ by value
  (World.hpp:651-655, 696): verified.

## Contradicted hypothesis
- Ticket Sec.2 lists a "RNG accessor" among public symbols. CONTRADICTED: there
  is NO public RNG accessor. `rng_` (std::mt19937_64), `gaussian_`, `uniform_`
  are private (World.hpp:696-698) and reached only through member moves. No doc
  written for a non-existent accessor.

## Hidden-global coupling (ticket Sec.6.4, for the Architect)
- `World::nDof()` (World.hpp:318-323) reaches
  `OpenMMContext::get().getNumDegreesOfFreedom()` on the Cartesian branch -- a
  process-global singleton not visible in the signature. The observable return
  is documented; the hidden coupling is flagged here for the thin-forwarder
  resolution ticket. Not resolved in the doc phase.

## Assumed notes
None. No `@note Assumed:` emitted.
