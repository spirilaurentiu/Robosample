# DOC-TopologyElements findings

Source: include/TopologyElements.hpp (SystemTopology). No code changed
(comment-stripped diff empty).

## Producer / consumer boundary — verified (sampled)

Producer: Python fills SystemTopology in place via Context; consumers:
World::buildModel/ModelBuilder and the OpenMM system builder. Verified sample
fields against ModelBuilder (src/world/ModelBuilder.cpp):
- atomsBegin[mol] is used as molecule mol's ROOT atom (line 140), in addition to
  being the molecule's atom-range begin. Documented both roles.
- rootMobilities[mol] is the per-molecule root joint (line 143), parallel to
  numMolecules.
- bondsI/bondsJ (parent->child), atomsMass, numMolecules, atomsX/Y/Z all consumed
  as the field-group comments state.
No producer/consumer disagreement found in the sampled fields.

## God-struct documented as-is

Per the ticket, the deferred sub-struct split was not anticipated. Added a
struct-level contract (producer/consumer, global-BFS atom index, per-section
parallel-array convention, CSR [begin,end) molecule ranges, INV-3 units) plus
brief per-field annotations for the previously bare fields (molecule ranges,
atom-name/atomic-number arrays, CMAP, NBFIX, GBSA, virtual-site count,
thermostat). Fields with pre-existing `///<` unit annotations were left as the
existing authors' evidence.

## Coverage gap

Full ~110-field producer verification against the Python side was not
exhaustive; only representative fields were traced to ModelBuilder and the
box-vector path. The remaining fields are documented at group granularity from
their inline annotations and the consumer's usage pattern; a full producer audit
is out of this pass's scope.

## Doxygen build config

No engine-level Doxyfile exists (only submodule Doxyfiles under
Robosample/openmm/Simbody01/pybind11). The doc-build configuration for the
engine is owned by the calibration ticket, not this one.

## Assumed: notes
None.
