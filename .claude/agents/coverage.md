# Coverage model

Line coverage alone misreports scientific code: a covered line can be guarded by a test that cannot
fail when the science breaks (the Rule 8 failure the `reviewer` hunts). Coverage here is three axes.

## The three axes

1. **Structural (gcov).** Was the line/branch executed. Measured on the `coverage` build (the Coverage
   build type, `--coverage`, no sanitizers, registers the gtest suite and builds the `.so`). gcovr
   merges the C++ `.gcda` with pytest-cov's XML. This is the axis the current `nox -s tests` targets;
   it is necessary and not sufficient.
2. **Oracle-binding.** Is every spec `PRECONDITION` / `INVARIANT` / `LEMMA` bound to a passing,
   discriminating test. A static checklist generated from the researcher's `validation:` block, not a
   gcov number. An `INVARIANT` with no bound test is a Blocking gap regardless of line coverage.
3. **Mutation.** The discriminator that catches the tautological test. Inject a fault into a covered
   line and rerun the fast sanitized suite; if nothing goes red, the line is *covered but unguarded*.
   The mutation operators are the entries of `numerical_hazards.md` — so a surviving mutant lands
   exactly where a silent physics bug would.

Structural coverage says a line ran. Mutation coverage says a wrong version of that line would be
caught. For scientific code the second is the real target.

## Mutation operators = the hazard checklist

Domain-agnostic operators (off-by-one, negate-condition, drop-statement) are the baseline. The
scientific operators come from `numerical_hazards.md` and are the ones that matter:

- sign flip on a `Phi` / `~Phi` term,
- swap `F` / `M` frame,
- reorder a pinned reduction (reintroduce non-determinism),
- drop a term from the `dH` sum,
- widen a subtraction that was arranged to avoid cancellation,
- remove a non-finite guard (let `NaN` propagate),
- off-by-one / signedness flip on an atom or DOF index.

A surviving mutant from this set is a missing or tautological oracle on a scientific path — the highest
-priority coverage gap, above any uncovered line.

## Tooling

- Structural: `gcovr` / `lcov` on the `coverage` build (exists).
- Mutation: `mull` (LLVM-IR mutation) for C++, `cosmic-ray` for Python, run against the fast sanitized
  suite. ASan kills any mutant that corrupts memory at no extra cost. Scope mutation to the diff in
  change-scoped mode; whole-module in campaign mode.
- Oracle-binding: a static check that every `validation:` item in the governing spec resolves to a test
  id, and that the test's `fails_on` is non-trivial. No new runtime.

## `coverage()` output

The `coverage()` tool returns all three: the gcov map, the oracle-binding checklist (bound / unbound
per validation item), and the mutation survivors (per file:line, with the operator that survived).

## Iterative loop

1. `coverage()` produces the three-axis report.
2. Rank gaps: surviving scientific mutants first, then uncovered scientific branches, then uncovered
   mechanical code, then uncovered dead code.
3. `reviewer` triages each gap into one bucket and routes it:
   - **dead code** -> `architect` `OPEN-QUESTIONS` (is it reachable at all?),
   - **scientific path, no discriminating oracle** -> `researcher`, to add a `validation:` item, then
     `coder` implements it,
   - **mechanical code, no test** -> `coder` writes the test directly.
4. Re-run `coverage()`.

Exit criterion (not a line percentage):

- every `INVARIANT` / `LEMMA` bound to a discriminating test, AND
- no surviving mutant on a scientific path, AND
- mechanical code line-covered.

Scientific code is held to mutation kill; mechanical code is held to line coverage. A badge that reports
only structural coverage SHALL be labeled as such, so a green number is not mistaken for guarded
science.

## Where each axis runs

- Mutation: fast sanitized Tests build — quick, and ASan is a free second oracle.
- Structural gcov: Coverage build (never the same binary as ASan).
- Oracle-binding: static, against the spec.
