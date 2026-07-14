# Project hazard checklist (Robosample)

The generic `reviewer` and `coder` agents describe hazard *classes*. This file names the *instances*
for this codebase. Agents load it as data when a change is classified scientific; they do not embed
it. When the science changes, this file changes — the agent prompts do not.

Each entry gives four things: the hazard class it instantiates, the specific check, the **oracle** that
closes it (a `verifier`/Wolfram template — computation-tier, so the proof is executed, not read), and
the **mutation** operator that probes whether a test actually guards it (`coverage.md`, axis 3). The
oracle is what a static finding carries so it can be machine-closed; the mutation operator is what
proves the guard is real rather than tautological.

## Theory adherence (hazard class: undocumented departure from established definitions)

- The change SHALL name the theory it preserves or disrupts. A departure from the documented
  background MUST be argued and documented.
- Re-derive the risky step in the codebase's notation. Cite the passage in `references/index.yaml`
  or `references/papers/*`; pull only what you cite.
- Oracle: none mechanical — this is provenance, checked by reading against corpus.
- Mutation: n/a.

## Correctness of the sampling move (hazard class: broken invariant in a transform)

- If sampling is touched: reversibility + volume preservation (or its Jacobian) SHALL be intact.
- The exact-`dH` acceptance test SHALL use the full Hamiltonian.
- No term SHALL land on the wrong side of the guidance/acceptance split.
- The Fixman / PMF distinction SHALL be respected.
- Oracle (volume preservation): symbolic Jacobian determinant of the transform in Wolfram; confirm it
  reduces to 1 (or to the Fixman factor) via `FullSimplify[det - expected] -> 0`. Proof-grade.
- Oracle (reversibility / detailed balance): the two-row verifier rule — prove the kernel ratio
  symbolically AND run the numeric reversibility test. The gap between them is the
  guidance/acceptance-split bug.
- Mutation: sign flip on the Jacobian term; move one term across the guidance/acceptance split. A
  surviving mutant means the acceptance test is not actually exercised.

## Numerical (hazard class: plausible wrong number that passes the suite)

The compiler axes are pinned in the build (`-fno-fast-math`, explicit `-ffp-contract`). These are
the algorithmic hazards no flag reaches, plus confirming the pins were not overridden.

- **Catastrophic cancellation in `dH`.** The acceptance test subtracts large energy totals for a
  small difference. Lost precision there biases acceptance silently — no crash, no failing test.
  Highest priority: the bias reaches every published result and nothing else detects it.
  Oracle: Wolfram computes the exact difference symbolically and the condition number of the
  subtraction; compare against the FP path at a configuration where totals are large and the
  difference small. The condition number is a stronger statement than any single FP test.
  Mutation: rewrite the compensated/reordered subtraction as a naive one; a test that still passes is
  not guarding the cancellation.
- **Summation / reduction order.** Reordered FP reductions are non-associative; OpenMP `reduction`
  and CUDA atomics vary run to run even when the race is benign, breaking determinism (Rule 9).
  Oracle: bit-for-bit determinism harness across two runs with identical seed (mechanical, not CAS).
  Mutation: swap a pinned pairwise/tree reduction for an unordered one; a green suite means no
  determinism test guards it.
- **NaN / Inf handling.** A non-finite energy reaching an acceptance comparison passes silently
  (`NaN < x` is false), biasing sampling instead of halting. A non-finite value SHALL halt (Rule
  11), not propagate.
  Oracle: a `PRECONDITION` test that injects a non-finite energy and asserts the run aborts with a
  diagnostic.
  Mutation: delete the non-finite guard; the injection test SHALL then fail.
- **Build-pin override.** Flag any per-target `-Ofast` / `-ffast-math` or changed `-ffp-contract` in
  a `CMakeLists.txt` or NVCC flag that re-enables what the toolchain pinned off.
  Oracle: static grep/scan of the build files (mechanical).
  Mutation: n/a (build-level, not a code path).

## Conventions (hazard class: silent convention fork)

The code MUST match exactly; these are not stylistic:

- Frame `F` vs `M`.
- Parent- vs body-frame angular velocity.
- Angular-over-linear ordering in spatial vectors.
- Sign conventions `Phi` vs `~Phi`.
- No explicit assembly of `M` / `M^-1` / `det M`.
- Oracle (mass operator): Wolfram SPD check and symbolic identity confirming the operator form matches
  `~Phi` without assembling `M^-1`; `FullSimplify[lhs - rhs] -> 0`. Proof-grade for the identity.
- Mutation: swap `Phi`/`~Phi`; flip the angular/linear ordering; a case where the wrong convention
  gives a detectably different result SHALL turn a test red.

## Boundary hazards the sanitizer sees poorly (hazard class: correctness the tool is blind to)

- **pybind11 boundary** — object ownership and lifetime across the C++/Python edge, reference
  counting, GIL assumptions.
- **Uninitialized reads on unexercised paths** — the sanitizer sees only the path the test ran.
- **Integer / sign confusion** in atom or DOF indexing that silently produces a wrong result.
- Oracle: none symbolic — reading plus, where a runtime trigger exists, an ASan/UBSan reproducer.
- Mutation: signedness flip on an index; drop a refcount increment. Reserved for campaign mode.
