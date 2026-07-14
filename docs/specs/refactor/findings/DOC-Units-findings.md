# DOC-Units findings

Source: include/Units.hpp. No code changed (comment-stripped diff empty).

## Units verified against a consuming call site (not the symbol name)

The engine's internal length unit is nm (INV-3). Confirmed at the DCD boundary:
OutputWriter multiplies engine coordinates by 10 to write Angstrom
(src/workflow/OutputWriter.cpp:82-84), matching nmToAngstrom's documented
direction (nm in, Angstrom out). The boundary policy (typed quantities only at
I/O edges; raw double nm/dalton/ps/kJ-per-mol in the core) is documented from the
file header and the conversion signatures.

## Context (not a claim about this module)

The kT300 / kB literal is duplicated ~10x in the statistical tests rather than
sourced from here; Units.hpp exposes only unit conversions, not physical
constants. Recorded as context per the ticket; no action in this pass.

## Coverage gap

No dedicated unit test; the conversions are used pervasively at I/O boundaries.

## Assumed: notes
None.
