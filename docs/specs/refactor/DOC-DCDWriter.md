# DOC-DCDWriter: CHARMM DCD trajectory writer

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** the self-contained CHARMM-format DCD trajectory writer - the cleanest
  module in the set, zero engine dependency (ARCHITECTURE section 7, MODULES.md section 5:
  move-as-is, no split).
- **Layer:** Infrastructure / io (ARCHITECTURE section 2).
- **Ownership:** owns its file handle for the writer's lifetime; borrows the
  per-frame coordinate buffer the caller passes. Recover header/frame ordering
  obligations (open -> header -> N frames -> close) from the call site in the output
  path.
- **Invariants (HYPOTHESES):** the on-disk layout is the CHARMM DCD binary contract
  (magic, header record, per-frame X/Y/Z blocks, box record when periodic). The
  writer's contract to callers is the frame cadence and the coordinate unit/frame it
  expects; document those, not the byte-packing mechanism.

## 2. Scope

- **Files:** `io/DCDWriter.{hpp,cpp}` (moved as-is by the io relocation).
- **Public symbols:** open/header/append-frame/close surface.
- **Known gaps to close:** confirm whether the caller (OutputWriter) supplies
  whole-molecule periodic-imaged coordinates or raw coordinates - the imaging lives
  in `OutputWriter`, not here (ARCHITECTURE section 3). State the boundary precisely.

## 3. Evidence pointers (tests exercising the module)

- No dedicated C++ gtest. Exercised end-to-end by the Level-1 representative run
  (`.dcd` output, VERIFY B3) and by `OutputWriter`. Record the unit-test gap in
  findings.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; contracts behavioral (frame
  cadence, coordinate convention, file-lifecycle order), not byte-layout narration.
- Comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-DCDWriter-findings.md` present; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
