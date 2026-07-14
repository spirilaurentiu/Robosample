# DOC-PeriodicBox findings

Source: include/PeriodicBox.hpp. No code changed (comment-stripped diff empty).

## O2 fold / pre-fold agreement — the check the ticket asked for

reducedBoxVectors() has NO production caller in the current tree. Production
consumes ParmEd's already-reduced box vectors directly:
- OpenMM setup: OpenMMSystemBuilder sets setDefaultPeriodicBoxVectors from
  SystemTopology::boxVectors (src/bridge/OpenMMSystemBuilder.cpp:73);
- DCD box record: OutputWriter::boxFromReducedVectors derives the CHARMM box from
  the SAME SystemTopology::boxVectors (src/workflow/OutputWriter.cpp:143).
Both consumers therefore read one identical reduced-vector source (agreement
confirmed by construction, not two divergent copies). TopologyElements documents
that ParmEd supplies the vectors already reduced, so no C++-side reduction runs
in production.

reducedBoxVectors() itself is exercised only by TestPeriodicBoundary (FAST
contract test), which pins the triclinic-construction + reduction convention.
Documented reducedBoxVectors as the reduction contract with an explicit note
that it has no production caller. The two historical copies the O2 ticket names
(the Context.cpp minimum-image lambda and
OpenMMContext::computePeriodicBoxVectors_Context) are already folded and their
predecessors removed, so a live pre-fold diff is not recoverable; the current
single source is internally consistent.

## minimumImage / minimumImageDistance

Production caller: StartupValidator clash scan
(src/workflow/StartupValidator.cpp:90). Test callers: TestPeriodicBoundary,
TestAlchemy. The reduced-cell precondition (single c->b->a pass is canonical only
for a lower-triangular cell) is documented as @pre.

## Assumed: notes
None.
