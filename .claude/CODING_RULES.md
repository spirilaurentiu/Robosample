# Robosample

## Coding rules

This is scientific software that requires special reasoning and coding constraints.

### Rule 0 - Adhere to theoretical background

This document explains scientific definitions, motivation, reasoning and default behavior.
A change may deliberately depart from it, but only if the departure is argued and its impact documented.

### Rule 1 - Think Before Coding

State **ALL** assumptions explicitly. Ask rather than guess.
Break the implementation into logical, bite-sized tasks or milestones.
Push back when a simpler approach exists. Stop when confused.

### Rule 2 - Simplicity First

Minimum code that solves the problem. Nothing speculative.
No abstractions for single-use code.

### Rule 3 - Surgical Changes

Touch only what you must. Don't improve adjacent code.
Match existing style. Don't refactor what isn't broken.

### Rule 4 - Goal-Driven Execution

Define success criteria. Loop until verified.
If in doubt, iteratively ask questions until requirements and edge cases are covered.
Strong success criteria let Claude loop independently.

### Rule 5 - Use the model only for judgment calls

Use for: classification, drafting, summarization, extraction.
Do NOT use for: routing, retries, deterministic transforms.
If code can answer, code answers.

### Rule 6 - Surface conflicts, don't average them

If two patterns contradict, pick one (more recent / more tested).
Explain why. Flag the other for cleanup.

### Rule 7 - Read before you write

Before adding code, read exports, immediate callers, shared utilities.
If unsure why existing code is structured a certain way, ask.

### Rule 8 - Tests verify intent, not just behavior

Tests must encode **WHY** behavior matters, not just WHAT it does.
A test that can't fail when business logic changes is wrong.
Output is valid only if you can name the theory your change **preserves** or **disrupts**.

### Rule 9 - Checkpoint after every significant step

Summarize what was done, what's verified, what's left.
Don't continue from a state you can't describe back.
Determinism is a feature, not a side effect.

### Rule 10 - Conventions in this codebase are load-bearing

Conformance > taste inside the codebase.
If you think a convention is harmful, surface it. Don't fork silently.
For instance, Frame conventions (`F` vs `M`, parent-frame vs body-frame angular velocity), index orderings (angular over linear in spatial vectors), sign conventions (`Phi` vs `~Phi` for forces vs velocities) are not stylistic.

### Rule 11 - Fail loud

"Completed" is wrong if anything was skipped silently.
"Tests pass" is wrong if any were skipped.
Always surface uncertainty.
Defaults that can be overridden when justified.

### Rule 12 - Self reflect

Review your changes before finishing the answer.
Act as your own hostile reviewer.

### Rule 13 - Programming languages, libraries and hardware

- C++17 on GCC 12.0 - 15.0 compiled with `x86_64-v3` (`SSE4.2`, `POPCNT`, `AVX`, `AVX2`, `BMI1` / `BMI2`, `FMA`, `LZCNT`, `MOVBE`) and LTO.
- Python 3.12
- CUDA 12.0 - 13.2
- OpenMP
- OpenBLAS and LAPACK available. CUDA counterparts available on CUDA capable devices (virtually all targeted machines)
- OpenMM 8.5 compiled directly into our source tree.

Code will run on consumer grade hardware (gamer PCs).
