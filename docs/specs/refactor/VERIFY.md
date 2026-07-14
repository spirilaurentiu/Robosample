# VERIFY - Behavioral-Preservation Verification Plan

Status: draft, 2026-07-12. Defines "done" for the whole refactor and per ticket.
Companion to [`README.md`](README.md). The refactor preserves behavior; this plan
is how that claim is made machine-checkable. If a check here cannot be run, the
ticket it guards is not complete.

---

## 1. Phase 0 - capture the Stage-0 baselines (blocking, before any ticket)

Run once on the current `disasm` tree and store under `docs/architecture/baseline/`.
No `SPLIT-###` ticket starts until these exist. Build env: an active conda env
(`CONDA_PREFIX` set); reference config `cuda-release` on the RTX 3090 / CUDA 13.0
that already runs the Level-1 example clean, plus `cpu-release` for the
sanitizer-friendly path.

- **B0 Exported public symbols.** `nm -C --defined-only` on the built
  `robo_bindings...so` (and each object library) -> `baseline/symbols.txt`. This is
  the API contract; the target delta on public symbols is **zero**.
- **B1 Pass set x3.** Build `*-tests`; run the full gtest suite three times.
  Identical results = deterministic baseline; any divergence = the flaky list,
  recorded now so flakiness is never misattributed to a later split ->
  `baseline/passset.txt`, `baseline/flaky.txt`.
- **B2 Assertion counts.** From gtest XML (`--gtest_output=xml`), per-test
  assertion count -> `baseline/assertions.txt`.
- **B3 Representative-run outputs.** The Level-1 example
  (`run.py ... 2ala.implicit ... --prod_steps 5 --seed 6000 --validate true`) stdout +
  the `2ala.implicit.*.csv`/`.dcd` it emits, plus the 4-replica REMC proof run
  (fixed seed) swap matrices -> `baseline/run_outputs/`. These diff to bit- or
  tolerance-identical after every ticket.
- **B4 Include graph + cycle report.** `clang -MM` (or CMake `--graphviz`) ->
  `baseline/include_graph.dot`; a cycle-detection pass -> `baseline/cycles.txt`
  (expected: the one managed RobotEngine<->Constraints forward-decl, ARCHITECTURE section 6.2).
- **B5 Coverage map (SHOULD).** gcov/lcov linking each public symbol to the tests
  that exercise it -> `baseline/coverage/`. Uncovered public symbols are recorded
  gaps, not silently ignored.

NOTE: the suite gates its million-sample cases behind `ROBOSAMPLE_SLOW_TESTS`.
Capture B1/B2 in BOTH modes (gated-off smoke and gated-on full) so a later change
to gating (Phase B `TEST-001`) has both baselines to diff against.

---

## 2. Per-`SPLIT-###` ticket exit criteria (Phase A)

Every SPLIT ticket is complete only when all hold:

1. **Build clean** in the reference config; no new warnings vs B0's build log.
2. **Pass set identical to B1**; **assertion counts identical to B2** (both
   modes). A red suite after a pure-motion ticket means the motion was not pure -
   revert, do not "fix forward."
3. **Representative outputs identical to B3** (bitwise for the deterministic HMC
   log; the OpenMM force-group validation stays within its existing 1e-6 band).
4. **Public-symbol nm diff vs B0: empty.** Any non-empty delta on a public symbol
   requires a separate human-approved API-change spec; it never rides inside a
   code-motion ticket. Internal-symbol churn (now `detail::`/anonymous) is allowed
   and expected.
5. **No new include cycles** (diff vs B4); no upward layer edges introduced
   (Application->Workflow->Domain->Algorithms->Bridge->Infrastructure->Platform;
   `MODULES.md section 1`).
6. **clang-format, clang-tidy, IWYU clean** on touched files.
7. **Comment-stripped before/after diff of every touched file is empty** -
   mechanically proves zero code change. (Strip comments + normalize whitespace,
   then diff; the only allowed non-empty hunks are the moved blocks appearing at
   their new location and the added `#include`.)
8. **File-size caps:** header <= 300 LOC, `.cpp` <= 600 LOC, unless the ticket
   claims a justified exception (parser / dispatch table / generated code).
9. **Predicted test breakage matched exactly.** The only test edits are the
   include-path fixes the ticket predicted; no other test file changed.

## 3. Dependency-defect resolution tickets (may change include structure)

The `MODULES.md section 3` defect resolutions (force-bridge interface, singleton
forwarders, ForceReducer dedup, joint-taxonomy centralization) narrow includes
and may add a thin forwarding symbol. Each such ticket additionally:

- states the before/after include edge it removes, checked against B4;
- if it adds or removes any *public* symbol, is flagged human-approved and
  carries its own nm-delta rationale (criterion 4 is relaxed only with that flag);
- for the ForceReducer dedup: adds a test asserting the host reducer and the CUDA
  `reduceForces` kernel produce identical per-body wrenches on a fixed input
  (closes the INV-1 divergence risk), and that test enters the B1 baseline before
  the dedup lands.

## 4. Per-`DOC-###` ticket exit criteria (doc phase)

- Doxygen builds **warning-free** (warnings-as-errors).
- Coverage script shows **every public symbol in scope documented**.
- **Comment-stripped diff of every touched file is empty** - proves the DOC ticket
  changed only comments, zero code.
- A `findings/DOC-###-findings.md` file is present (even if empty); no
  `@note Assumed:` without a matching findings entry.
- Contracts state behavior/intent, never mechanism: no sentence becomes false
  under a behavior-preserving rewrite. Contracts require call-site unanimity;
  disagreement between callers, or with ARCHITECTURE.md, is a finding, not prose.

## 5. Per-`TEST-###` ticket exit criteria (Phase B)

- Build clean; **pass set identical to B1**; **test + assertion counts identical
  to B2**, or matching an explicit human-approved triage delta recorded in the
  ticket.
- **Coverage map (B5) not reduced.** B5 is a SHOULD baseline (section 1), but for
  the Phase-B test refactor it is a hard exit: a test move that drops coverage is a
  regression, not a reorganization.
- **No test includes a `.cpp`.**
- Three consecutive full-suite runs identical (re-establish determinism).
- Contract comments are the deliverable of `TEST-007` (mirror-structure) and its
  exit criterion there. `TEST-001`/`004`/`005`/`006` are moves and optimizations
  that SHALL NOT be blocked on contract-comment coverage; they defer it to
  `TEST-007`, which brings every remaining test to the `TESTS.md section 5.9`
  structure.
- For `TEST-001` (slow-tier) specifically: the **mutation-power gate** -
  each shortened statistical oracle is re-run against a deliberately biased engine
  (wrong Jacobian sign, swapped beta, dropped Fixman term) and SHALL still fail at the
  pre-change false-fail rate. A shortening that lets a known bias pass is rejected.

---

## 6. Whole-effort done

The refactor is complete when: every ticket meets its criteria above; the full
build is clean in all reference configs; B1/B2/B3 hold end-to-end; the
public-symbol set equals B0 (modulo human-approved API specs); no file exceeds its
cap without a claimed exception; the include graph has no cycles or upward edges;
Doxygen builds warning-free; and the CI enforces the file-length check, the
include-cycle script, and the Doxygen-warnings-as-errors build so the structure
cannot silently regress.
