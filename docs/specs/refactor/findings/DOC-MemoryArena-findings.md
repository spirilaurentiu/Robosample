# DOC-MemoryArena findings

Source: include/MemoryArena.hpp. No code changed (comment-stripped diff empty).

## Ownership / lifetime contract — recovered from the sole consumer

The arena owns its slab and frees it on destruction/move-assign; handed-out
pointers are borrowed and valid only while the arena is alive and un-reassigned.
Recovered from RobotState's use (include/RobotState.hpp): allocateFull/
allocateCompact do `arena_ = MemoryArena{}` to reset before re-reserving, which
frees the previous slab and invalidates all previously handed-out cache
pointers. Documented reserve()-before-commit() protocol and the reserve-after-
commit throw and idempotent commit() from the code.

## Alignment guarantee is relied upon (documented as contract)

kAlign = 64 and every sub-array is 64-byte aligned. This is not merely
incidental: RobotEngine::realizePosition runs `#pragma omp simd` loops over arena
arrays (src/RobotEngine_kinematics.cpp:150), for which the aligned layout is a
correctness/legality assumption. Documented the alignment as a guarantee.

## Coverage gap

No dedicated MemoryArena unit test. Exercised indirectly by every
RobotState-building test (TestKineticEnergy, TestMassMatrix, the oracles).
Recorded per the ticket.

## Assumed: notes
None.
