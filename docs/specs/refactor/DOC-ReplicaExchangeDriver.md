# DOC-ReplicaExchangeDriver: label-swap RunREX driver

Status: draft (written pre-split, EXECUTED in the doc phase after all `SPLIT-###`
verified). Executor: `documenter`. Style: `styles/reference.md`. Exit criteria machine-checked (VERIFY section 4).

## 1. Context (HYPOTHESES - verify against call sites)

- **Purpose:** `RunREX` - the label-swap replica-exchange driver: `R` replicas on a
  temperature ladder, alternating-parity neighbour swaps via `SwapAcceptance`, the two
  inverse index maps, and the swap matrices. **The intended center of the engine** -
  the refactor makes REMC `RunREX` the default driver and `R = 1` its degenerate case
  (ARCHITECTURE section 1, section 3; MODULES.md C4).
- **Layer:** Workflow (ARCHITECTURE section 2).
- **Ownership:** owns the replicas (`Replica`/`ThermodynamicState`, DOC-Replica), the
  two index maps, the swap matrices, and `rexRng_` (ARCHITECTURE section 4). Drives each
  replica through `GibbsSweep` (DOC-GibbsSweep) between swap attempts.
- **Invariants (HYPOTHESES):**
  - **INV-8 REX detailed balance** (ARCHITECTURE section 5): `RunREX` and legacy `runREX`
    produce equivalent sampling; the two index maps stay mutually inverse; the swap
    matrix has the correct alternating-parity structure (validated: 4-replica REMC on
    CUDA, ARCHITECTURE section 1). Document the driver's contract: what a round does (sweep ->
    attempt swaps -> update maps/matrices), that labels swap while configurations stay
    in place, and that the product distribution is left invariant.
  - Scope note: the runnable surface is **REMC + Default**; RENE/REBASONTOP driven-REX
    is uncompiled (OQ-5) and RENEMC's round-loop throws (TESTS.md section 7). Document only
    the runnable REMC/Default path with confidence; mark the rest as not-runnable and
    route to findings.
  - The behavioral (stationary-distribution / detailed-balance) oracle for the full
    chain is **pending** (TESTS.md section 7, a separate track before the Context-REX code
    motion) - today's evidence is algebraic (SwapAcceptance) plus the label-swap
    equivalence. Do not claim the driver samples each replica's Boltzmann distribution
    as established; cite it as the equivalence-to-`runREX` contract (INV-8), not as a
    proven stationary-distribution guarantee.

## 2. Scope

- **Files:** `workflow/rex/ReplicaExchangeDriver.{hpp,cpp}` (`RunREX`, the maps, the
  swap matrices).
- **Public symbols:** `RunREX`, the map/matrix accessors, the swap-round entry.
- **Known gaps to close:** state the round structure and the inverse-map/parity
  invariants; keep the REMC/Default vs uncompiled-driven-REX boundary explicit; do not
  overstate the (pending) stationary-distribution guarantee.

## 3. Evidence pointers (tests exercising the module)

- TestRexAcceptanceAlgebra (acceptance, not the driver; TESTS.md section 7),
  `test_rex_label_swap_equivalence.py` (label-swap vs coordinate-swap INVARIANT-EQUIV,
  Python - out of C++ scope but the equivalence evidence). The end-to-end behavioral
  REX oracle is the recorded gap (TESTS.md section 7). Record it in findings.

## 4. Exit criteria

- Doxygen warning-free; every public symbol documented; INV-8 + map/parity contract
  stated; the runnable-vs-uncompiled boundary and pending-oracle status recorded;
  comment-stripped diff empty (VERIFY section 4).
- `findings/DOC-ReplicaExchangeDriver-findings.md` present; the pending behavioral
  oracle and OQ-5 boundary recorded; no unmatched `@note Assumed:`.

## 5. Calibration

Attach the approved DOC-CALIBRATION trio once human-reviewed; match their form; spec
wins any conflict; conflict is a finding.
