# DOC-DCDWriter findings

Source: include/DCDWriter.hpp. No code changed (comment-stripped diff empty).
DCDWriter.cpp was not touched: it defines only the header-declared methods and
introduces no additional public symbol, and documentation lives at the
declaration.

## Coordinate / imaging boundary — resolved (the ticket's known gap)

The caller, workflow/OutputWriter, converts engine nm to Angstrom (x10) and
performs whole-molecule COM periodic imaging BEFORE calling Writer::append
(src/workflow/OutputWriter.cpp:82-84 scatter with *10; lines 102-139 COM
imaging). The Writer itself performs no unit conversion and no imaging: it stores
the coordinates verbatim as 32-bit float and treats the Box as already in the
CHARMM on-disk convention (Angstrom lengths, degree angles, produced by
OutputWriter::boxFromReducedVectors). Documented this boundary on the class and
on append(): callers supply Angstrom, imaged coordinates.

## Latent Doxygen warning fixed

The header had a `@param path/numAtoms/withBox` block attached to the paramless
default constructor (`Writer() noexcept = default;`) — the parameters belong to
initialize(). Moved the parameter documentation onto initialize() (comment-only
change) so Doxygen does not warn about parameters on a paramless function.

## Coverage gap

No dedicated C++ gtest for the DCD writer. Exercised end-to-end by the Level-1
representative run (.dcd output) and by OutputWriter. Recorded per the ticket.

## Assumed: notes
None.
