# DOC-Units: unit constants and conversions

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the unit constants and conversion factors used across the engine
  (nm / kJ*mol^-^1 / ps conventions), wrapping the nholthaus units platform dependency
  (ARCHITECTURE section 2, section 7; MODULES.md section 5: single-concept, no split).
- **Layer:** Infrastructure / util (ARCHITECTURE section 2).
- **Ownership:** none - compile-time constants and pure conversion functions.
- **Invariants (HYPOTHESES):** the inter-world currency is per-atom Ground
  coordinates in **nm** (INV-3, ARCHITECTURE section 5); the constants here define that unit
  system. Document each constant's unit and reference frame so callers cannot
  misread a factor; verify the actual unit against a consuming call site, not the
  symbol name (names are not evidence, `documenter.md` section 3.5).

## 2. Scope

- **Files:** `util/Units.hpp`.
- **Public symbols:** the exported constants and conversion helpers.
- **Known gaps to close:** none structural. State the physical unit and frame for
  each constant; if two call sites disagree on a factor's meaning, that is a finding.

## 3. Evidence pointers (tests exercising the module)

- None dedicated. Used pervasively; the `kT300` / `kB` literal is duplicated 10x in
  the statistical tests (TESTS.md section 4) rather than sourced here - record that as
  context, not a claim about this module.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; comment-stripped diff empty
  (VERIFY section 4).
- `findings/DOC-Units-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
