# Tier 1 — configurational PE-distribution ladder (pytest, end-to-end)

Scope: does each sampler produce a canonical PE distribution, and do DOF-matched
samplers agree. Molecule: ala-dipeptide (examples/ala-dipeptide.*, implicit or
vacuum) — small, has real torsional metastability, fixtures exist.

## Claim
(i) Every Robosample sampler run at two temperatures obeys the Shirts PE-slope
−(β₂−β₁) (foundations §5.1). (ii) Native OpenMM Langevin MD, the Robosample
OpenMM-only world, and the Robosample OpenMM+GCHMC composite all span the SAME
(full, flexible) DOF, so their PE distributions are equal.

## Ground truth
Native OpenMM Langevin MD (openmm.app.AmberPrmtopFile/AmberInpcrdFile on the
Reference platform, matched exactly as openmm_validation._build_reference_system
does) = external anchor. This reuses the single-point PE machinery already trusted
by test_openmm_potential_energy.py, now extended to a *distribution*.

### Rung A — native OpenMM MD PE distribution (anchor)  [new file]
- Langevin/BAOAB, matched FF/T/nonbonded/constraints to the Robosample worlds.

### Rung B — Robosample OpenMM-only world PE distribution  [new file]
- add_cartesian_world + add_sampler(..., AlwaysAccept or MetropolisHastings).
- INVARIANT (B≡A): KS/χ² equality of PE histograms (DOF-matched, §4). This is
  framing #4 (native ↔ Robosample-OpenMM).

### Rung C — Robosample OpenMM+GCHMC composite PE distribution  [new file]
- add_cartesian_world + add_robotic_world (Fixman ON), Gibbs-alternated.
- INVARIANT (C≡B): KS/χ² equality of PE histograms. This is framing #3 (energy
  leak). DOF-matched because the Cartesian world relaxes bonds/angles.

### Slope check on every rung
- LEMMA: for each of A, B, C, run at T₁,T₂ (ΔT/T ≈ √(2k_B/C_V), §6) and fit the
  Shirts slope by ML logistic regression (§5.1 Estimator B) AND weighted-linear
  (Estimator A) on N_eff. Expect z = (α̂₁+Δβ)/se within 3σ, consistently across
  ≥3 seeds. This validates canonical-ness independent of B≡A≡C.

## Data source
PE time series from {base}.moves.csv column `PE`; subsample to N_eff via
batstat.py/autoblock.py before any statistic. Native side logs its own PE series.

## Correctness conditions
- PRECONDITION: A, B, C use identical FF, T, nonbonded method+cutoff+PME tol,
  GBSA model, and constraint set (§6). A mismatched constraint set changes the
  ensemble and is a setup bug, not a sampler bug.
- PRECONDITION: histogram-equality (B≡A, C≡B) is asserted ONLY here because these
  three are DOF-matched; it SHALL NOT be extended to a pure-internal world (§4.1,
  framing #5-total is rejected).
- INVARIANT: B≡A and C≡B at KS p > 1e-3 (or χ² at alpha=1e-4) on N_eff samples.
- LEMMA: per-rung Shirts slope = −(β₂−β₁) within 3·se(N_eff).
- Moments (⟨PE⟩, Var(PE)) equal within 4·stderr(N_eff) as a coarse pre-filter.

## Touch list
- New tests/test_ensemble_pe_ladder.py (pytest, gated like
  test_openmm_potential_energy.py: skip w/o openmm/extension, hard-fail under
  ROBOSAMPLE_REQUIRE_OPENMM).
- Reuse openmm_validation._build_reference_system for the matched native System.
- New shared helper tests/ensemble_stats.py (Shirts ML + linear slope, N_eff via
  batstat). Conventions at risk: Shirts sign (§5.1), N_eff vs N, moves.csv KE/PE
  column semantics (see Open questions).

## Verification plan
Discriminating structure: B≡A alone can pass while C is silently wrong; the
C-vs-B leak test is what exercises the GCHMC world's effect on the *configurational*
distribution. A plausible-but-biased Fixman (wrong sign, §3 trap) leaves B≡A intact
but breaks C≡B and shifts the per-rung C slope — so both the equality AND the slope
on rung C must be checked together. Cross-check: the C-world Shirts slope must hold
even where C≡B is marginal, isolating "canonical but shifted" from "non-canonical".
