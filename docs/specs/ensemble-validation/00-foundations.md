# Statistical-mechanics validation suite — Foundations

Shared vocabulary, conventions, estimators and correctness principles for the
tiered ensemble-validation tests. Tiers 0/1/2 (`10-*`, `20-*`, `30-*`) reference
this file and do not redefine any term here.

The suite answers one question in three layers: **does each Robosample sampler
draw configurations and momenta from the canonical (NVT) ensemble it claims, and
do two samplers that are supposed to target the *same* ensemble agree?** This is
"correctness against theory, not intuition" (CLAUDE.md).

> RESUMING? Read `90-implementation-status.md` first — it records what is written,
> what passes, and the open Tier-2 slow-tier failures (with diagnoses and the next
> steps to run). This file and the tier files (`10/20/30-*`) state *what is correct*;
> the status file tracks *what is done*.

---

## 1. Worlds, coordinate systems, and the two kinetic-energy paths

- **Cartesian / OpenMM world** — `Context.add_cartesian_world(...)`. Coordinates
  are atomic Cartesian positions `x ∈ R^{3N}`. Potential energy `U` and kinetic
  energy `KE` are both computed by OpenMM. The mass matrix is the constant
  diagonal Cartesian matrix `M_cart` (atomic masses). Fixman is OFF (flat metric).
- **Internal / GCHMC / torsional world** — `Context.add_robotic_world(...)` /
  `add_torsional_world(...)`. Coordinates are joint (mobilizer) generalized
  coordinates `φ`. `U` is still computed by OpenMM on the atomic positions
  reconstructed from `φ`; `KE` is computed by the internal articulated-body engine
  (`RobotEngine::calcKineticEnergy`) from the configuration-dependent generalized
  mass metric. Fixman is ON by default (`World.cpp:256`).

**Generalized mass metric** `M(φ)` (a.k.a. mass-metric tensor): the `n_dof × n_dof`
matrix `M(φ) = Jᵀ M_cart J`, `J = ∂x/∂φ` (Spiridon 2017 `eq:mass-metric-tensor`).
`RobotEngine::calcLogDetM` returns `ln det M(φ)` via the O(n) articulated-body
(LDLᵀ / innovations) factorization. In this suite `M(φ)` always denotes this true
(non-inverted) metric.

> NOTE (notation trap). Spiridon 2020 (`spiridon_2020_robosample/equations.md`
> eq:2,4,6) writes `M_tot` for the **inverse** metric and prints the marginal as
> `|M_tot|^{-1/2}`; the file's own CHECK comment flags this. Under the binding
> convention here (§3) that is the same `+1/2` power of the true metric. Do not
> propagate the `-1/2` literal into any test.

**Kinetic energy.** For velocities `u` (generalized) with `KE = ½ uᵀ M(φ) u`, or
momenta `p = M(φ)u` with `KE = ½ pᵀ M(φ)⁻¹ p`. The internal engine seeds
`u = √(RT)·M^{-1/2} g`, `g ~ N(0, I)` (`TestEquipartition.cpp` THEORY block).

---

## 2. Degrees of freedom (`n_dof`) — counting rules

`n_dof` is the number of momentum/velocity coordinates actually drawn and is the
`n` in every Gamma/chi-square/equipartition expression below.

- **Internal world, tree topology (acyclic):** `n_dof = Σ_joints mobility(joint)`.
  Joint mobilities: `Torsion` = 1; `Free` (6-DOF floating root) = 6; a **welded**
  root contributes 0. `RobotModel::nu` is this count (`TestEquipartition.cpp`
  uses `m.nu`).
- **Loop closure / holonomic constraints (RATTLE):** each independent distance
  (rod) constraint removes exactly one velocity DOF *in the M-metric*
  (`enforceVelocityConstraints`). So `n_dof = nu − n_C` (`TestEquipartition.cpp`
  `ConstrainedRing`: 8 − 1). This is the constraint-bookkeeping the Fixman
  loop term (`calcConstraintLogDet`, `World.cpp:860`) also depends on — the two
  MUST use the same `n_C`.
- **Cartesian / OpenMM world:** `n_dof = 3N − n_constraints − n_com`, where
  `n_com = 3` iff centre-of-mass motion removal is active and `n_constraints`
  counts SHAKE/SETTLE constraints. This is OpenMM's own DOF count and is a
  **distinct code path** from the internal one; a Cartesian-world KE test SHALL
  use the OpenMM DOF count, not `RobotModel::nu`.

> PRECONDITION. Any KE-distribution test SHALL obtain `n_dof` from the same engine
> that produced the KE samples. `n_dof` is exposed to Python as the read-only
> `World.n_dof` property (internal worlds: `nu − n_C`; Cartesian worlds: the OpenMM
> `3N − n_constraints − n_com` count). A Python-level KE-shape / equipartition test
> is therefore possible; it remains gated only on the KE-column semantics of
> `moves.csv` (Open Q1: resampled vs end-of-trajectory KE). The authoritative
> Gamma/equipartition tests live in C++ (Tier 0); the Python tiers otherwise use the
> DOF-free statistics (PE-slope, T-ratio).

---

## 3. The configurational marginal and the Fixman convention (binding)

Integrating the Gaussian momenta out of `ρ(φ,p) ∝ exp(−β[½ pᵀM(φ)⁻¹p + U])`
gives the momentum-integration factor `(det M(φ))^{+1/2}`:

    ρ(φ) ∝ (det M(φ))^{1/2} · exp(−βU(φ))          (binding; Spiridon 2017 eq:2;
                                                    TestFixmanBoltzmann.cpp THEORY)

**Production Fixman** (`World::calcFixman`, `World.cpp:835`):

    U_F(φ_f) = ½ RT · ( ln|M_{N_f}(φ_f)| − ln|M_3N| ),   |M_3N| = const

Adding `U_F` to the acceptance Hamiltonian multiplies the marginal by
`exp(−βU_F) = |M_{N_f}|^{-1/2}·const`, so:

    Fixman OFF:  ρ_OFF(φ_f) ∝ (det M_{N_f})^{1/2} · exp(−βU)      (metric-weighted)
    Fixman ON:   ρ_ON (φ_f) ∝                       exp(−βU)      (flat-metric)

Because `M_3N` (full Cartesian metric determinant) is torsion-independent, this is
the simplified pure-torsion Fixman of Kandel 2016 `eq:13`
(`kandel_2016_fixman_hybrid_icmd`): `ρ_f(α) ∝ exp(−U(α,q₀)/kT)`, the Boltzmann
distribution of the potential evaluated at the **frozen** bond/angle geometry `q₀`.

> INVARIANT (already covered by TestFixmanBoltzmann.cpp, restated for reuse). With
> `U=0` and a φ-dependent `det M`, ρ_OFF matches `(det M)^{1/2}` and rejects flat;
> ρ_ON matches flat and rejects `(det M)^{1/2}`. Tier 2 extends this to `U≠0`
> against an external reference.

---

## 4. The DOF-matching principle (the load-bearing correction)

Two samplers' **distributions may be compared for equality only if they span the
same degrees of freedom and the same accessible configuration space.** The Shirts
*slope* test (§5.1) is exempt — it validates one sampler against its own density
of states — but every histogram-equality / KS / moment-equality comparison is
subject to this rule.

Consequences that reshape the proposed tiering:

1. **A rigid-body torsional world does not have the same total-PE distribution as
   a fully flexible Cartesian model, even with Fixman ON.** Frozen bonds/angles
   contribute a *constant* to `U` (their geometry does not fluctuate), whereas a
   flexible Cartesian model carries ~`kT/2` of vibrational PE per bond/angle DOF.
   The internal-world PE distribution is therefore *narrower*; the totals cannot
   coincide. (`U` is always the full OpenMM potential, so this is exact, not an
   artefact.)
2. **Fixman corrects the metric (kinetic) bias only, not the rigid-vs-flexible
   configurational difference.** Kandel 2016 `eq:11`,`eq:18–19`: the Fixman-
   corrected constrained torsion marginal equals the fully-flexible marginal only
   in the *separable* limit (potential factorizes into torsion-only and
   bond/angle-only parts). Real force fields couple them through 1-4 and nonbonded
   terms (`eq:20` Coulomb test; `eq:22` FF form), leaving a residual **extrinsic**
   distortion that is physically real, not a bug. König 2014
   (`konig_2014_constraint_correction`) quantifies this residual as the free-energy
   cost of the constraint (`eq:21`,`eq:29`).

Therefore a cross-simulator distribution *equality* test requires DOF-matching. Two
recipes achieve it, but **neither is required by this suite** (see the NOTE below):

- **(M-full)** run the internal world with all bond/angle DOF flexible (full BAT,
  `det ratio ≈ 1`, Fixman ≈ const) and compare to *unconstrained* Cartesian MD; or
- **(M-rigid)** constrain the Cartesian reference (SHAKE) to the same rigid
  geometry the internal world freezes, and keep Fixman ON on the internal side.

For the common rigid-body torsional robot vs standard flexible OpenMM MD, neither is
in force, so the external comparison is **qualitative** (topology of minima,
symmetry, relative populations within a stated band), not an equality test.

> NOTE (optional; not a gate). The M-full and M-rigid recipes exist only to sharpen
> the *external* Tier-2 comparison from qualitative to exact. The suite does not
> depend on them, for two reasons: (i) an all-flexible BAT world is a test-only
> regime, not Robosample's production torsional sampler; (ii) Tier-2's T2.0 already
> provides an *exact internal* oracle (self-consistency vs `exp(−βU_full(φ;q₀))`),
> which catches Fixman/metric bugs with no external-simulator noise. M-full would
> only additionally check that the frozen-geometry approximation converges to full
> flexibility — an approximation-quality question, not a sampler-correctness one.
> Whether `build_flexibilities` can even make bonds/angles mobilizers is unconfirmed;
> treat M-full as an opt-in follow-up, not part of the initial suite.

---

## 5. Shirts ensemble-validation estimators

Shirts 2013, *Simple Quantitative Tests to Validate Sampling from Thermodynamic
Ensembles*, JCTC 9(2):909–926, DOI 10.1021/ct300688p (local
`shirts_2012_ensemble_validation`; preprint 2012 / published 2013). Reference
implementation: `physical_validation` (Merz, Volkhardt & Shirts), built on Merz &
Shirts 2018, PLoS ONE 13(9):e0202764, DOI 10.1371/journal.pone.0202764.

### 5.1 Two-temperature PE-slope test (configurational; DOF-agnostic)

For any Boltzmann ensemble the log-ratio of the PE histograms at two temperatures
is linear in `E` with a slope fixed by the temperatures and **independent of the
unknown density of states** (Shirts `eq:6`,`eq:20`; Spiridon 2017 `eq:6`):

    ln[ P(U|β₂) / P(U|β₁) ] = (α₀) − (β₂ − β₁)·U ,   slope α₁ = −(β₂ − β₁)

- **Null hypothesis H0:** both datasets are canonical at their stated `T`, i.e.
  `α₁ = −(β₂ − β₁)`.
- **Estimator A (weighted linear, binned):** fit `y_k = ln(n_{k,2}/n_{k,1})` vs bin
  centre `E_k` by weighted least squares over bins populated in both runs, weights
  `w_k = 1/var(ln r_k)`, `var(ln r_k) = 1/n_{k,1} + 1/n_{k,2} − 1/N₁ − 1/N₂`
  (Shirts `eq:30`). Slope covariance `(XᵀWX)⁻¹` (`eq:wls`).
- **Estimator B (maximum-likelihood, histogram-free; preferred):** logistic
  regression `ln L = Σ_{T₁} ln f(−α₀−α₁E_i) + Σ_{T₂} ln f(α₀+α₁E_j)`,
  `f(x)=1/(1+e^{-x})` (Shirts `eq:8`). Concave → unique maximum. Reparameterize
  with `β_ave = (β₁+β₂)/2` fixed so the two free parameters are `(β_ave ΔA, Δβ)`;
  report the fitted **effective ΔT** and its Fisher-information std error
  `var(α₁) = −1/(∂²lnL/∂α₁²)` (`eq:36`). This is the "effective-temperature MLE".
- **Sign discipline.** `eq:6` uses `P₂/P₁` (slope `−Δβ`); `eq:7` residual uses
  `P₁/P₂` (slope `+Δβ`). Fix one convention per implementation
  (`shirts .../equations.md` CHECK on eq:7).

### 5.2 KE-Gamma test (momentum-space; needs `n_dof`)

Freshly drawn momenta `p ~ N(0, kT·M(φ))` give the quadratic form
`pᵀM(φ)⁻¹p / kT ~ χ²_{n_dof}` **exactly and independent of φ** (any covariance
`Σ=kT M`: `pᵀΣ⁻¹p ~ χ²`). Hence

    KE ~ Gamma(shape = n_dof/2, scale = kT) ,  density ∝ E^{n/2−1} e^{−E/kT}

in **both** coordinate systems. Goodness-of-fit by Pearson χ² / G-test at
`alpha = 1e-4` (`StatTest.hpp`), with sub-sampled per-bin Gamma weights
(`TestEnsembleValidation.cpp::KEisMaxwellBoltzmann`, which must not use bin-centre
evaluation for this curved density). Merz & Shirts pair this with the equipartition
check to catch DOF-count and RNG errors.

> This tests the mass metric `M(φ)`, the `M^{-1/2}` momentum draw, the `kT`
> scaling and the `n_dof` count. It does NOT test accept/reject, the integrator,
> the potential, or detailed balance — those are Tiers 1/2.

### 5.3 Equipartition (first moment; §5.2 collapsed)

    ⟨2KE⟩ / n_dof = kT   (absolute; needs n_dof)
    ⟨KE(T₁)⟩ / ⟨KE(T₂)⟩ = T₁/T₂   (ratio; DOF-free, weaker)

`⟨KE⟩ = (n_dof/2)kT` holds **even for configuration-dependent `M(φ)`**: the
generalized equipartition theorem gives `⟨p_i ∂H/∂p_i⟩ = kT` per canonical
coordinate (Jain 2012 `jain_2012_icmd_equipartition`, Tolman derivation), and
`Σ_i p_i ∂H/∂p_i = 2KE`; equivalently the modal principle `⟨vvᵀ⟩ = kT·I` (Jain
2012 eq:26). The **per-coordinate** split fails in physical `(φ,p)` coordinates
(Jain's central result) — only the **total** KE obeys the simple law. State the
claim as total-KE equipartition, never per-torsion.

---

## 6. Statistical hygiene (applies to every distribution comparison)

- **Matched physics.** Compared runs SHALL share force field, `T`, nonbonded
  method + cutoff + PME tolerance, GBSA/implicit model, and constraint set.
  Constraints change **both** `n_dof` and the accessible configuration space
  (§4), so an unmatched constraint set silently changes the target ensemble.
- **Autocorrelation-aware error bars.** HMC/MD samples are correlated. Every error
  bar and every Shirts weight (§5.1) SHALL use the effective sample size
  `N_eff = N/g`, `g = 1 + 2τ_int` the statistical inefficiency, estimated per
  series (`python/robosample/batstat.py`, `autoblock.py`). Using raw `N` makes a
  correct sampler fail (error bar too tight). This is the single most common
  false-fail cause.
- **Choice of two-sample statistic.**
  - *Shirts slope* (§5.1): primary for "is this canonical / do two `T` agree".
    Robust to the unknown density of states; use for every per-sampler check.
  - *Moment comparison* (`MeanAccumulator`, means at `4·stderr` with `N_eff`):
    coarse pre-filter; cheap, catches gross errors before the distribution test.
  - *KS / χ² histogram equality*: only under DOF-matching (§4); use for
    DOF-matched sampler pairs (native OpenMM ↔ Robosample-OpenMM) and torsion
    marginals.
- **Tolerance derivation (not round numbers).**
  - Means: `|obs − exp| ≤ 4·stderr(N_eff)` (`StatTest.hpp`; two-sided false-fail
    ~6e-5).
  - Distribution GoF: χ²/G at `alpha = 1e-4`, with `N` (or `N_eff`) chosen so power
    > 0.99 against the *specific* wrong hypothesis (e.g. flat vs `(det M)^{1/2}`),
    documented at the call site (`StatTest.hpp` tolerance philosophy).
  - Shirts slope: z-score `z = (α̂₁ − α₁_true)/se(α̂₁)` with `se` from the
    covariance in §5.1 using `N_eff`; flag when `|z| > 3` **consistently across
    repeats** (Shirts checks.md line 134: "> 2–3 σ consistently"). The 25%
    magnitude bound in `TestEnsembleValidation.cpp::PEobeysEnsembleSlope` is a
    stand-in for runs that do not compute `se`; the z-score form supersedes it
    where `se` is available.
  - Temperature gap: choose `ΔT/T ≈ √(2k_B/C_V)` so the two PE means differ by
    ~2σ_E (Shirts `eq:gap-T`); too small → no slope signal, too large → no
    histogram overlap.

---

## 7. Corrections to the originating five framings (summary; details in §3–§5)

1. *Equipartition / `⟨KE⟩=(n_dof/2)kT` for `M(φ)`* — **correct**, cite Jain 2012;
   but it is **total-KE only**, not per-coordinate (§5.3). "Momenta Gaussian with
   covariance `M·kT`" is *conditionally* Gaussian in both systems; only in
   Cartesian is `M` constant so the *unconditional* momentum law is one Gaussian.
   The KE quadratic form is `χ²/Gamma` in both regardless (§5.2).
2. *KE `Gamma(n_dof/2, kT)` in both systems* — **correct** (§5.2). Scope: mass
   matrix + RNG + `T` + DOF; not dynamics. Confirmed.
3. *Energy-leak (OpenMM-only ↔ OpenMM+GCHMC same PE)* — **valid and DOF-matched**:
   in the composite the OpenMM world relaxes bonds/angles, so both span full DOF.
   This is the correct use of a PE-equality test (Tier 1 rung B-vs-A).
4. *Torsion target is the Cartesian-marginal PMF incl. `√det M`* — **partially
   correct, must be split** (§3,§4): the `√det M` factor is the *metric* piece
   Fixman handles; but the fully-flexible Cartesian marginal *also* contains
   bond/angle thermal coupling that a rigid-body torsional world lacks even with
   Fixman. Exact equality needs DOF-matching; otherwise the comparison is
   qualitative. The valid **self-consistency** reference is
   `exp(−βU_full(φ; q₀))` at frozen geometry (= what Fixman-ON GCHMC targets by
   construction, Kandel eq:13); the "bare **dihedral-term-only** scan" is
   generally **invalid** because it omits 1-4/nonbonded contributions to the
   torsion PMF.
5. *Native Cartesian PE dist ↔ Robosample pure-internal HMC PE dist* — **invalid
   as a total-PE equality** across DOF mismatch (§4.1): the rigid model has no
   bond/angle vibrational PE. Replace with: (a) per-sampler Shirts slope on each
   (both must give `−Δβ`), and (b) equality only on the shared torsion marginal
   under the Tier-2 caveats, or under DOF-matching.

The tiering is otherwise sound; the tier files below encode these corrections.

---

## 8. Open questions (status)

- **Q1 — `moves.csv` KE-column semantics (OPEN, gating one optional test).** Is the
  logged `KE` the freshly-resampled (start-of-move, `World.cpp:1703`) or the
  end-of-trajectory (`World.cpp:1772`) kinetic energy? The KE-Gamma / equipartition
  identity (§5.2–5.3) holds for the *resampled* draw. The C++ Tier-0 tests sidestep
  this by calling `seedMomenta` directly; only the optional Python-level KE smoke
  (Tier 0 T0.4) is blocked on confirming this.
- **Q2 — `n_dof` exposure to Python (RESOLVED).** Now available as the read-only
  `World.n_dof` property (§2 PRECONDITION). Unblocks Python-level equipartition.
- **Q3 — all-flexible (M-full) BAT world feasibility (OPEN, non-gating).** Whether
  `build_flexibilities` can make bonds/angles mobilizers determines whether the
  *optional* exact Tier-2 variant is buildable (§4 NOTE). The suite passes without
  it; this is an opt-in follow-up, not a prerequisite.
