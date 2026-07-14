# Spec: Make label-swap RunREX REMC the default driver; R=1 as the degenerate case

Status: **draft** — behavioral change spec, ready for review. Author: orchestrator,
2026-07-12, from the RunREX validation run (4-replica REMC on CUDA) and the REX
port campaign. This is NOT pure code motion; it changes the runnable default and
is human-approved separately from the refactor SPLIT tickets.

Requirement keywords per RFC 2119 (`styles/spec.md`): SHALL / SHOULD / MAY / NOTE.

Scope: the Python driver layer (`run.py`, `roborun.py`, the `run_*.py` scripts)
and, if needed, a thin `Context` convenience entry point. No sampling-algorithm
change. Both drivers already exist and are bound; this rewires which one ships.

---

## Motivation

### Problem restatement (user's phrasing)

> unblock runrex. robosample SHALL be centered around replica exchange. current
> implementation describes only a degenerate case (1 replica only).

### In the codebase's vocabulary

Two REX drivers are bound to Python (`src/PyBind11.cpp:588,597`):

- `run_rex` → `Context::runREX` — **coordinate-swap** REMC. Swaps `replicaCoords_`
  between replicas; INV-3 non-compliant by construction; marked "do not extend"
  (`Context.hpp:84-90`). Retained as the differential oracle.
- `run_rex_label_swap` → `Context::RunREX` — **label-swap** REMC. Coordinates stay
  per replica; thermodynamic-state labels swap; carries the `Replica` /
  `ThermodynamicState` object model and all run types.

Every runnable driver script calls `run_rex` at `NOF_REPLICAS = 1`
(`roborun.py:210`, `run.py:171`, all `run_*.py`). So what ships is a single
replica with no exchange — a degenerate REX. The label-swap `RunREX` — the real
multi-replica REX — is exercised only by tests.

**Verified prerequisite (2026-07-12):** the CUDA build runs the Level-1 example
clean (RTX 3090, CUDA 13.0), and a 4-replica REMC ladder [300, 356, 422, 500] K
runs end to end via `run_rex_label_swap(RunType.REMC, …)` producing the correct
alternating-parity neighbour swap matrix. The engine is not the blocker; the
driver wiring is.

### Constraints

- **Behavior preservation for the R=1 path.** A single-replica run SHALL sample
  the same distribution it does today. `R = 1` is the degenerate case of REX, not a
  separate code path.
- **Equivalence oracle stays.** `run_rex` (coordinate-swap) is retained as the
  INVARIANT-EQUIV oracle (ARCHITECTURE OQ-3); this spec does not delete it.
- Precedence (`CLAUDE.md`): the sampling distributions are fixed by theory
  (parallel tempering); this spec changes only which correct driver is the default.

---

## Behavior

### B1. `run_rex_label_swap` REMC becomes the default driver

The driver scripts SHALL call `run_rex_label_swap(RunType.REMC, equil, prod,
write_freq, verbose)` in place of `run_rex(...)`. The temperature ladder SHALL come
from `context.initialize(temperatures)` with `len(temperatures) = R ≥ 1`.

### B2. R=1 reduces to single-replica HMC, exactly

For `R = 1`, `RunType.REMC` attempts no swaps (`prepareExchangePairs` yields no
neighbour pair when `T = 1`; `mixReplicas` with `T ≤ 1` is a no-op). The run is
then `equil + prod` Gibbs sweeps of the one replica at its one temperature — the
same sampler the current single-replica `run_rex` drives. `RunType.Default` (no
exchange attempts) is the explicit alias for "independent replicas / no swaps" and
MAY be used for `R = 1` to make the intent legible. Either SHALL reproduce the
current single-replica behavior.

### B3. Default ladder policy

A driver that requests `R > 1` SHALL build a **geometric** ladder
`T_i = T_0 · r^i`, `r = (T_max/T_0)^{1/(R−1)}`, so swap acceptance is roughly
uniform across rungs (Kofke geometric-ladder result; see the oracle spec
`rex-stationary-distribution-oracle.md`). `roborun.py` already computes exactly
this (`roborun.py:211`); the change is to raise `NOF_REPLICAS` from 1 and route to
the label-swap call. `T_0`, `T_max`, and `R` SHALL be driver parameters, not
hardcoded constants.

### B4. Output indexing unchanged

`RunREX` already indexes output CSV/DCD by **thermodynamic state**, matching
`run_rex`'s convention that a fixed output slot is a fixed temperature
(`Context.hpp:120-125`). No output-format change; the reaction-force CSV and
per-replica DCD cadence are preserved.

### B5. Run-type scope

Only `RunType.REMC` and `RunType.Default` are wired as default-driver options.
`RunType.RENE`/`REBASONTOP` require an `mdSteps = 0` distort-world setup and are
out of scope here (their driver wiring is a separate spec); `RunType.RENEMC`
throws (Stage 2c). A driver SHALL reject an out-of-scope run type with a clear
error rather than silently degrading.

### Derivation sketch (why the swap is safe to default)

Label-swap and coordinate-swap parallel tempering induce the same distribution
over (configuration, temperature) assignments: swapping two replicas' coordinates
at fixed temperatures and swapping their temperature labels at fixed coordinates
are the same move up to relabeling (Sugita–Okamoto; Chodera–Shirts). `runREX` is
retained precisely as the differential oracle that this equivalence holds
(INVARIANT-EQUIV, already asserted bit-identically for REMC in
`tests/test_rex_label_swap_equivalence.py`). Therefore defaulting to
`run_rex_label_swap` REMC preserves the sampled distribution while giving the
per-rung continuous trajectories, per-rung schedules, and object model the
coordinate-swap driver cannot.

---

## Invariants

- **INV-D1 (R=1 equivalence).** A single-replica `run_rex_label_swap` REMC/Default
  run SHALL produce the same sampled distribution as the current single-replica
  `run_rex` run. Same seed → statistically identical observables.
- **INV-D2 (multi-replica equivalence to the oracle).** For `R ≥ 2` REMC, the
  per-thermodynamic-state statistics SHALL match the retained `runREX`
  coordinate-swap oracle (INVARIANT-EQUIV), within the equivalence test's tolerance.
- **INV-D3 (no format drift).** Output CSV/DCD/reaction files SHALL keep their
  current schema and by-state indexing.

---

## Interface

### Touch list

- `python/robosample/run.py` — replace `initialize([300])` + `run_rex(...)` with a
  ladder + `run_rex_label_swap(RunType.REMC, ...)`; expose `--n_replicas`,
  `--t_min`, `--t_max` (defaults: 1, 300, 300 → degenerate, matching today).
- `python/robosample/roborun.py` — set `NOF_REPLICAS` from an argument (default 1);
  route line 528 to `run_rex_label_swap(RunType.REMC, ...)`.
- The remaining `run_*.py` GPCR/ligand scripts — same substitution; keep their
  current single-replica defaults so no scientific run silently changes until its
  owner opts into a ladder.
- `Context` (optional convenience, MAY): a `run_rex(...)`-signature-compatible
  wrapper `run_replica_exchange(run_type=REMC, ...)` so a driver need not name the
  internal method. No behavioral effect; thin forward.

### Externally observable change

The default REX driver becomes label-swap. With the default `R = 1` the observable
behavior is unchanged (INV-D1). A driver that sets `R > 1` gets true multi-replica
REMC. No output schema changes (INV-D3).

---

## Validation strategy

Oracles tagged PRECONDITION / INVARIANT / LEMMA.

- **VD1 (PRECONDITION) — build + smoke.** The CUDA build runs the Level-1 example
  clean (already verified 2026-07-12). Re-confirm after the driver edit.
- **VD2 (INVARIANT, INV-D1) — R=1 equivalence.** Run the current single-replica
  `run_rex` and the new single-replica `run_rex_label_swap` REMC on the same system
  and seed; the deterministic HMC accept/reject trace and the emitted CSV SHALL be
  bit- or tolerance-identical.
- **VD3 (INVARIANT, INV-D2) — multi-replica equivalence.** The existing
  `tests/test_rex_label_swap_equivalence.py` already asserts label-swap ≡
  coordinate-swap REMC bit-identically; extend/confirm it covers the default-driver
  configuration (`R = 4` geometric ladder).
- **VD4 (LEMMA) — swap sanity.** For `R ≥ 2`, the attempted-swap matrix is
  nonzero, nearest-neighbour, alternating-parity; accepted/attempted ∈ (0,1). (The
  4-replica proof run already exhibits this.)
- **VD5 (LEMMA) — behavioral oracle (when landed).** The
  `rex-stationary-distribution-oracle.md` oracle (per-state marginal vs closed-form
  Boltzmann) validates the driver samples the right distribution — the strongest
  check, sequenced with that spec.

NOTE: VD2/VD3 are the load-bearing checks; they hold the "R=1 unchanged, multi-R
equivalent to the trusted oracle" contract that makes defaulting to the label-swap
driver safe.

---

## Consequences and trade-offs

- **The scientific `run_*.py` scripts keep single-replica defaults** until each
  owner opts into a ladder, so no production run's distribution changes silently.
  Centering on REX is a default-and-capability change, not a forced migration.
- **Retiring `run_rex`** (OQ-3) is deliberately out of scope: it remains the
  equivalence oracle. Retirement is a later human-gated decision once the behavioral
  oracle (VD5) is trusted.
- **Driven run types (RENE/REBASONTOP/RENEMC)** are not defaulted; their driver
  wiring + `mdSteps=0` distort-world setup + Stage-2c RENEMC drive are separate
  specs.

## Open questions

- **OQ-D1** Should the convenience `Context::run_replica_exchange` wrapper be added
  now (legibility) or deferred to avoid widening the public API before the Context
  REX split (`MODULES.md C4`)? Recommendation: defer; bind after the split so the
  new entry point lands on the extracted `ReplicaExchangeDriver`.
- **OQ-D2** Default `T_max` and `R` for drivers that opt into a ladder are
  system-dependent; leave as required parameters (no silent default ladder) or ship
  a conservative default (e.g. `R = 4`, `T_max = 1.3·… `)? Recommendation: required
  parameters — a silent default ladder is a scientific choice the driver owner must
  make.
