# Spec: Torsional-dynamics efficiency and inter-basin crossing for the composed Gibbs/HMC chain

Status: **draft** — researcher artifact, ready for review. Author: derived from `references/papers/spiridon_2017_cdhmc_gibbs`, `references/papers/spiridon_2020_robosample`, the reference library `references/index.yaml`, the in-flight `docs/specs/rex-stationary-distribution-oracle.md` and `docs/specs/replica-exchange-nonequilibrium-work.md`, and the existing Robosample infrastructure (`python/robosample/autoblock.py`, `include/NMA.hpp`, `include/World.hpp`, `include/HMCSampler`/`RobotIntegrator`, `include/BatScaling`, the NCMC world).

Requirement keywords follow RFC 2119 (`styles/spec.md`): SHALL / SHOULD / MAY / NOTE.

Scope: this is a **research catalog and prioritization**, not an implementation ticket. It changes no source. It defines the candidate space, the correctness conditions each candidate SHALL satisfy to remain composable with the unbiased blocked-Gibbs/HMC chain, and the oracles that distinguish a correct implementation from a plausible-but-biased one. Individual candidates promoted from the shortlist SHALL each get their own implementation spec.

---

## Motivation

### Problem restatement (user's original phrasing)

> Produce an exhaustive, evaluated catalog of candidate methods for FASTER CONVERGENCE of the composed Markov chain in Robosample on BOTH axes — intra-basin decorrelation AND inter-basin rare-transition crossing. Currently the team optimizes only INTRA-basin mixing, via correlation-based block selection (Islam & Venugopal clustering of correlated torsions). They want to also accelerate INTER-basin rare-transition crossing (basin hopping), and to improve torsional-dynamics efficiency generally.

### The gap in the codebase's vocabulary

Robosample composes one Markov chain from a fixed scan of **worlds**. A world is a robot factorization exposing a subset of generalized (BAT) coordinates as free and welding the rest (`include/World.hpp`, `spiridon_2020_robosample` §2.1.1). Each sampled torsional world runs generalized-coordinate HMC (GC-HMC): draw generalized momenta `u` via `sqrt(M_φ^{-1})` at the internal-coordinate equipartition (`jain_2012_icmd_equipartition`), integrate constrained Velocity-Verlet on the unmodified potential, and Metropolize on the full Hamiltonian `H = PE + KE + U_F − ½ RT·ln(sin²…)` where `U_F` is the Fixman term and the last term is the torsion-measure Jacobian (`include/World.hpp`, `spiridon_2017_cdhmc_gibbs` §2.2). A Cartesian OpenMM world at low duty cycle supplies ergodicity over the welded DOF (`spiridon_2020_robosample` §2.1.1, condition 2). A ROUND is a scan of worlds; a TRAJECTORY is many rounds; only the composed trajectory converges to Boltzmann `π`.

The convergence rate of that composed chain has two independent bottlenecks, and the team currently addresses only the first:

1. **Intra-basin decorrelation** — the integrated autocorrelation time (IACT) of a within-basin observable. The current lever is correlation-based block selection: `python/robosample/autoblock.py` runs TICA → KMeans microstates → MSM → PCCA+ metastable states, then `context.chose_correlated_bonds` community-detects strong/weak/rogue torsion blocks from a dihedral correlation matrix (the "Islam & Venugopal clustering" referenced by the user). Highly-correlated torsions are co-blocked into one world so a single HMC move rotates them together.

2. **Inter-basin crossing** — the mean first-passage time (MFPT) or the number of independent A↔B crossings per unit compute. `spiridon_2020_robosample` Table 4 and `spiridon_2017_cdhmc_gibbs` §3.4 already measure this (MFPT to the isolated α_L basin). Reduced torsional/Ramachandran worlds cut MFPT across the φ=0 barrier by an order of magnitude versus fully-flexible MD. But nothing in the current scan is *designed* to cross barriers; crossing is an incidental benefit of long torsional moves.

**Reformulation.** The user's two asks map to two distinct estimator-variance problems on the composed round-level chain:

- *Torsional-dynamics efficiency* = maximize effective samples per unit compute of intra-basin observables = raise the per-world HMC timestep and move length without lowering acceptance, and equalize exploration across the retained soft DOF.
- *Inter-basin crossing* = minimize MFPT_AB (equivalently maximize independent A↔B crossings per unit compute) **without** changing the stationary populations `π(A)`, `π(B)` or `ΔF_AB`.

These are coupled by one shared correctness invariant: every world SHALL be exactly `π`-invariant, so that a systematic scan is `π`-invariant regardless of order (`chodera_2011_gibbs_replica_exchange` casts exactly this composition as Gibbs sampling on an augmented state). A method that raises crossing rate by biasing `π` is out-of-constraint (Critical severity, `CLAUDE.md`).

### Binding constraints (from the task; challenged only with published evidence)

1. **Unbiased.** No bias potential or reweighting in the base method; `π`-invariance / detailed balance is Critical. Biased methods are admitted only with exact validated reweighting, else classified out-of-constraint.
2. **Reduced-soft worlds are the crossing engine.** Sampled worlds expose only soft (torsional) DOF. Introducing a stiff (bond/angle) DOF into a sampled world caps the HMC timestep at the fastest retained mode and destroys the large-step advantage (`spiridon_2017_cdhmc_gibbs` §1: torsional MD stable to ~10 fs; `mazur_1991_icmd`: freezing fast DOF buys 9–13× timestep). Stiff/omitted DOF are handled by the OpenMM ergodicity world plus Fixman.
3. **Basin-dependent structure.** The all-vs-all torsion correlation matrix changes basin-to-basin. A block scheme fitted to one basin is mis-specified elsewhere, and within-basin clustering can SPLIT a concerted reaction coordinate across separate Gibbs blocks, manufacturing a kinetic barrier (see Behavior, family C).
4. **Composition correctness.** A systematic scan of exactly-`π`-invariant worlds is `π`-invariant regardless of order; order is an efficiency knob, not a correctness knob — *provided each world is exactly `π`-invariant*.
5. **Cost is not binding.** Torsional worlds and short OpenMM are cheap. The scarce resource is convergence speed per unit compute.
6. **Scale.** Up to 1M atoms, 100k rigid bodies, 10k robots. Global statistical estimation over all coordinates is infeasible; estimation SHALL be local/hierarchical.

---

## Behavior

The catalog below is the centerpiece. Candidates are grouped by mechanism family (A–H). Each candidate carries: **Mechanism**, **Axis** (intra / inter / both), **Correctness** (does it preserve unbiased `π`; if not, the exact fix), **Fit** to reduced-soft-torsion + worlds/Fixman/ABA, **Scale**, **Cost/pilot**, **Evidence tier**, **Composability**. A prioritized shortlist follows.

Evidence tiers: **[Published]** = theory in `references/index.yaml` or external DOI; **[Inference]** = derived here from published theory but not directly asserted by a cited source; **[Intuition]** = plausible, unverified.

### Standing correctness lemma (used by every candidate)

**L0 (dual-Hamiltonian freedom).** HMC unbiasedness requires only (i) momenta drawn from the Gaussian matching the kinetic metric used in the acceptance Hamiltonian, and (ii) a reversible, volume-preserving proposal, with acceptance `min(1, exp(−β ΔH_acc))` on the acceptance Hamiltonian `H_acc` (`duane_1987_hmc`; `spiridon_2017_cdhmc_gibbs` §1.3: "the integrator need not preserve Boltzmann — only the forward/reverse propagation ratio matters"). The *guidance* Hamiltonian used to propagate MAY differ from `H_acc`. Therefore any method that only changes the propagator (softened forces, fictitious inertia, SGLD guiding force, larger step) preserves `π` **iff** `H_acc` retains the true `U + U_F − ½RT ln sin²` and the momentum term is self-consistent with the drawn momenta. This is the single lever most crossing-accelerators pull, and its violation (Forrest & Suter omitting `U_F`, `forrest_1994_gchmc`; `spiridon_2017_cdhmc_gibbs` §1.2) is the canonical bias bug here. **[Published]**

---

### Family A — HMC kernel improvements (raise timestep / reduce random walk)

**A1. Look-ahead / rejection-suppressing HMC (LAHMC).**

- *Mechanism.* Replace the Metropolis accept/reject of the GC-HMC move with the LAHMC greedy transition operator: on a would-be rejection, extend the trajectory by additional leapfrog segments and accept a further state, choosing among candidate end-states by probabilities that satisfy a fixed-point (generalized detailed balance) equation. This suppresses momentum flips and the resulting random-walk back-tracking (`sohldickstein_2014_lahmc`).
- *Axis.* Both. Fewer flips → more ballistic excursions (helps cross) and lower IACT (helps decorrelate).
- *Correctness.* Preserves `π` exactly via the LAHMC fixed-point (generalized detailed balance), a published drop-in replacement for Metropolis HMC. No reweighting. Momentum-flip-on-reject bookkeeping SHALL be consistent with the constrained integrator's time-reversal.
- *Fit.* Drop-in at the acceptance layer of the torsional world; leaves Fixman/ABA untouched. Compatible with constraint (2): does not add DOF.
- *Scale.* Per-world cost multiplied by the look-ahead depth (small constant). Fine at 1M atoms because torsional worlds are cheap (constraint 5).
- *Cost/pilot.* No pilot; a step-size/look-ahead-depth tune.
- *Evidence.* **[Published]** `sohldickstein_2014_lahmc`.
- *Composability.* Combines with A2, B*, C*, E* (it only changes the inner kernel). Conflicts with nothing.

**A2. Partial momentum refreshment / generalized (underdamped) HMC.**

- *Mechanism.* Instead of full momentum resampling each move, refresh momenta partially: `u ← cosθ·u + sinθ·ξ`, `ξ ~ N(0, k_BT M_φ)`, small `θ`, with the Horowitz momentum-flip-on-reject. Directional persistence across successive short moves produces a longer ballistic path than repeated independent HMC moves of the same length.
- *Axis.* Both, weighted to inter-basin (ballistic transport suppresses the `sqrt(t)` random walk).
- *Correctness.* Preserves `π` with the flip-on-reject correction (standard generalized-HMC / Horowitz result). NOTE: partial refreshment interacts with the Gibbs world-swap — between worlds the coordinate set changes, so persistent momenta SHALL be re-projected or discarded at world boundaries; carrying momenta across a factorization change is not defined and SHALL NOT be done.
- *Fit.* Modifies momentum draw only; Fixman/ABA untouched.
- *Scale.* Free.
- *Cost/pilot.* Tune `θ`; no pilot.
- *Evidence.* **[Published]** general HMC theory; related to `sohldickstein_2014_lahmc`. **[Inference]** for the world-boundary re-projection rule.
- *Composability.* Combines with A1 (LAHMC subsumes some of this), B*, E*.

**A3. NUTS / Exhaustive-HMC (XHMC) adaptive trajectory length.**

- *Mechanism.* Terminate each trajectory dynamically (no-U-turn / virial-rate threshold) instead of a fixed step count `L` (`betancourt_2016_xhmc`). Removes the per-world `L` tuning that `spiridon_2020_robosample` §2.1.2 does by hand.
- *Axis.* Intra-basin primarily (picks the integration time that decorrelates within the current basin). NOT a barrier-crossing accelerator: the U-turn criterion terminates at the basin scale, not the barrier.
- *Correctness.* Preserves `π` (published; already partially wired via `use_nuts` in `python/robosample/run.py`).
- *Fit.* Already present as a sampler flag. Constraint-friendly.
- *Scale.* Fine; already used.
- *Cost/pilot.* None.
- *Evidence.* **[Published]** `betancourt_2016_xhmc`; already in-tree.
- *Composability.* Combines with B* (mass preconditioning changes the U-turn geometry favorably). Mild conflict with A2 (both manage trajectory length/persistence — pick one).

**A4. Delayed-rejection / look-ahead step-size ladder.**

- *Mechanism.* On rejection of a large-step proposal, attempt a smaller-step proposal from the same point before giving up, with the delayed-rejection acceptance correction.
- *Axis.* Intra-basin (recovers acceptance at aggressive steps).
- *Correctness.* Preserves `π` with the standard delayed-rejection ratio.
- *Fit.* Acceptance-layer only.
- *Scale.* Free.
- *Cost/pilot.* None.
- *Evidence.* **[Published]** (external, Mira 2001 delayed-rejection MCMC; not in `index.yaml`). **[Inference]** for the CDHMC mapping.
- *Composability.* Redundant with A1; choose one rejection-suppression scheme.

---

### Family B — Mass-metric / preconditioning (equalize frequencies, align moves with soft modes)

**B1. Constant mass-matrix preconditioner from the existing NMA Hessian.**

- *Mechanism.* Replace the equipartition momentum draw's metric with a *constant* fictitious metric `M̃` chosen to equalize the retained-mode frequencies: draw `u ~ N(0, k_BT M̃)`, propagate with `M̃`, accept on `H_acc` that uses `½ uᵀ M̃⁻¹ u` **plus the unchanged true Fixman term `U_F` built from the real `M_φ`**. `M̃` is built from `include/NMA.hpp` `computeRouteBNMA` (mass-weighted internal-coordinate Hessian eigenpairs already computed at a minimized `q0`): set each mode's fictitious inertia so `ω_k` are uniform, so a single timestep is no longer capped by the stiffest retained torsion.
- *Axis.* Intra-basin strongly (frequency equalization → larger uniform step, uniform exploration rate). Inter-basin **only if** the barrier lies along a retained soft mode (the preconditioner reshapes kinetic energy, not the potential barrier).
- *Correctness.* Preserves `π` by L0: a *constant* `M̃` has constant `det M̃`, which cancels in acceptance, so only `U_F` (true metric) governs the configurational marginal. This is exactly Forrest & Suter's constant fictitious inertia (`forrest_1994_gchmc`) made rigorous by keeping `U_F` (their omission was the bug, `spiridon_2017_cdhmc_gibbs` §1.2). **Correctness condition (INVARIANT):** the acceptance kinetic term SHALL use the same `M̃` used in the draw, and `U_F` SHALL remain the true-`M_φ` Fixman term; a position-dependent `M̃` would additionally require a `½ ln det M̃(q)` bookkeeping term (see B2) and SHALL NOT be introduced as if constant.
- *Fit.* Reuses `computeRouteBNMA`; ABA untouched; no new DOF (constraint 2 satisfied — same coordinate set, different inertia).
- *Scale.* NMA is `2ν` OpenMM force evals per world at build time (`include/NMA.hpp` cost note), ν = #retained torsions — cheap and local (constraint 6). Recompute per basin, not globally.
- *Cost/pilot.* One NMA build per world per basin.
- *Evidence.* **[Published]** `forrest_1994_gchmc`, `jain_2012_icmd_equipartition` (modal decoupling), `vitalis_2014_mixed_rigid_dihedral` (diagonal-metric decoupling). **[Inference]** for building `M̃` from the Route-B Hessian specifically.
- *Composability.* Combines with A1/A2/A3 and E*. Substitutes partially for finer blocking (see B-note).

**B-note (mass matrix vs finer blocking).** A constant preconditioner and finer Gibbs blocking both target the same pathology — different effective stiffnesses across retained DOF. The preconditioner equalizes them *inside one world* at O(ν) build cost and one acceptance term; finer blocking splits them across *more worlds* at the risk of the concerted-splitting failure (family C). **[Inference]:** where the slow structure is a stiffness spread rather than a topological reaction-coordinate split, B1 is the cheaper and safer lever and SHOULD be tried before subdividing blocks.

**B2. Riemannian-manifold HMC (position-dependent metric).**

- *Mechanism.* Use a full position-dependent metric (e.g. the true `M_φ`, or a softened Hessian) with the generalized (implicit) leapfrog, adding the `½ ln det M(q)` term to `H_acc` (`brubaker_2012_chmc` frames constrained/RMHMC in one family).
- *Axis.* Intra-basin.
- *Correctness.* Preserves `π` with the log-det-metric term and the implicit-integrator reversibility. NOTE: Robosample's true `M_φ` is *already* the physical metric, and its `U_F = ½ RT ln det M_φ` is *already* the log-det term — so "RMHMC with the physical metric" is essentially what the engine does. An *additional* artificial position-dependent metric would double the non-separable bookkeeping for little gain over B1.
- *Fit.* Poor: implicit integrator is expensive and fragile at scale; conflicts with constraint (5)'s preference for cheap large steps.
- *Scale.* Poor (fixed-point solves per step).
- *Cost/pilot.* High.
- *Evidence.* **[Published]** `brubaker_2012_chmc`. **[Inference]** for the redundancy with existing Fixman.
- *Composability.* Do not combine with B1 (competing metrics).

**B3. Diagonal-mass / modal decoupling (Fixman-free-by-construction propagator).**

- *Mechanism.* Propagate with the diagonal of the mass-metric tensor (per-DOF effective masses), integrating each generalized DOF independently (`vitalis_2014_mixed_rigid_dihedral`; modal coordinates `jain_2012_icmd_equipartition`). This decouples the equations of motion and removes `det G` from the *dynamics*.
- *Axis.* Intra-basin (large stable steps, no coupled inertia inversion).
- *Correctness.* As a *guidance* propagator inside CDHMC it is unbiased by L0 provided `H_acc` keeps the true `U_F`. WARNING: `vitalis_2014_mixed_rigid_dihedral` uses the diagonal metric to *avoid* the Fixman correction entirely ("configurational sampling free of MMT artifacts by construction, at the cost of correct dynamics"); adopting that claim *without* `U_F` in acceptance would bias `π` for coupled DOF. In Robosample the safe reading is diagonal-mass-as-guidance + true Fixman-in-acceptance.
- *Fit.* Good as guidance; ABA already yields the diagonal inertias.
- *Scale.* Good.
- *Cost/pilot.* Low.
- *Evidence.* **[Published]** `vitalis_2014_mixed_rigid_dihedral`, `jain_2012_icmd_equipartition`.
- *Composability.* Overlaps B1 (both are inertia-shaping); pick one guidance inertia.

---

### Family C — Gibbs blocking and scan design (the concerted-splitting failure and its fix)

**C1. Co-blocking the concerted set via a fixed union scan (the differential reaction block).**

- *Mechanism.* Diagnose the concerted set — the torsions whose coupling *flips sign or magnitude* between basin A and basin B — from a **differential** correlation/mutual-information analysis across basin endpoints, not a within-basin correlation. `python/robosample/autoblock.py` already computes transition-window correlation matrices (`transition_correlation_matrix`, `plot_per_transition_corr`) and exposes `context.compute_differential_correlation`. Put the concerted set into one world (one HMC move rotates them together). Then carry the **union** of {per-basin blocks, the differential reaction block, fine soft blocks} in a single *fixed* systematic scan. In a basin where a block does not match the local correlation structure, that block's move is a cheap near-no-op; in the matching basin it delivers the benefit.
- *Axis.* Inter-basin primarily (enables the concerted motion that a split-block scan forbids), with intra-basin benefit from the per-basin blocks.
- *Correctness.* Preserves `π` exactly: every world is `π`-invariant, so the fixed union scan is `π`-invariant regardless of order (constraint 4; `chodera_2011_gibbs_replica_exchange` composition). Crucially this avoids the **adaptive-selection reversibility hazard**: state-dependent choice of which block to run breaks reversibility unless done as a valid adaptive-MCMC scheme (diminishing adaptation; external: Roberts & Rosenthal 2007, Latuszynski et al. 2013). The fixed union scan needs no adaptation conditions.
- *Fit.* Excellent — reuses the autoblock pipeline; each block is still a reduced-soft torsional world (constraint 2).
- *Scale.* The differential analysis is local (per transition region), not global (constraint 6). Union size grows with #basins × #reaction-coordinates; bounded because near-no-op blocks are cheap (constraint 5).
- *Cost/pilot.* Requires at least the two basin endpoints (structures of A and B), and ideally one observed crossing to define the differential couplings — the bootstrap loop of family D.
- *Evidence.* **[Inference]** built from Gibbs-blocking theory (blocking correlated variables together shortens Gibbs mixing time; external: Liu 1994, Roberts & Sahu 1997 — not in `index.yaml`) and constraint 4. The differential-correlation machinery is already in `autoblock.py`.
- *Composability.* Combines with everything (it is the scan design). This is the backbone the other families plug into.

**C2. The concerted-coordinate-splitting failure mode (what C1 fixes) — stated as a hazard, not a candidate.**

- *Mechanism of failure.* If a concerted reaction coordinate (e.g. two torsions `θ1, θ2` that must move together to relieve a 1-4 clash) is split into separate Gibbs blocks by within-basin correlation clustering, then each block conditions on the other frozen at reactant values. The conditional barrier seen by block-`θ1` with `θ2` fixed is `≥` the true minimum-energy-path (MEP) barrier, because the MEP is off-axis to both single-coordinate cuts. Gibbs samplers of strongly correlated targets mix as `~1/(1−ρ²)` in the correlation `ρ`; near a concerted transition `ρ→1` and per-block conditional crossing is exponentially suppressed in the manufactured barrier height. The scan can thus be *slower* than fully-flexible MD at crossing, while still being unbiased in `π`.
- *Correctness.* This is not a `π` bug — populations remain correct in the ergodic limit — it is a *rate* bug (the composed chain is `π`-invariant but slowly mixing). It is Critical-adjacent because it defeats the software's primary goal (rare transitions).
- *Evidence.* **[Inference]** from Gibbs-sampling theory + `spiridon_2017_cdhmc_gibbs` §4 (1-4 vdW clashes gate low-T torsional transitions; "modified torsional terms or softened vdW can help").
- *Discriminator.* See Validation LEMMA V-C (alanine cis→trans / two-torsion clash gate).

**C3. Overlapping / partially-collapsed Gibbs blocks.**

- *Mechanism.* Let a concerted torsion appear in *more than one* block (overlap), or "collapse" nuisance coordinates by delegating them to the OpenMM ergodicity world (partial collapse). Overlapping blocks improve mixing of correlated targets without a single monolithic block.
- *Axis.* Both.
- *Correctness.* Overlapping systematic scans remain `π`-invariant (each world `π`-invariant). Partial collapse is exactly the existing "fully-flexible world for ergodicity" (`spiridon_2020_robosample` §2.1.1).
- *Fit.* Good; the OpenMM world is the collapse mechanism already.
- *Scale.* Fine.
- *Cost/pilot.* None beyond C1's diagnosis.
- *Evidence.* **[Inference]** (external: partially-collapsed Gibbs, van Dyk & Park 2008; not in `index.yaml`).
- *Composability.* A refinement of C1.

**C4. Adaptive scan order / adaptive block selection (admitted only under diminishing adaptation).**

- *Mechanism.* Choose which block to run based on the current state (e.g. run the reaction block more when near a transition).
- *Axis.* Inter-basin.
- *Correctness.* HAZARD: naive state-dependent selection breaks reversibility and can bias `π`. Admissible only as valid adaptive MCMC with **diminishing adaptation + containment** (external: Roberts & Rosenthal 2007; Latuszynski, Roberts & Rosenthal 2013; Chimisov et al. adaptive Gibbs 2018 — not in `index.yaml`). The fixed union scan C1 is the safe substitute and SHOULD be preferred unless a diminishing-adaptation proof is supplied.
- *Fit.* Risky; conflicts with the "order is not a correctness knob *provided each world is π-invariant*" framing (constraint 4) because adaptive selection makes the kernel state-dependent.
- *Evidence.* **[Published]** adaptive-MCMC correctness theory (external).
- *Composability.* Do not combine with C1 unless the adaptation demonstrably vanishes.

---

### Family D — Reaction-coordinate / slow-mode identification for the reaction block

**D1. tICA / VAMP / MSM / PCCA+ (already in `autoblock.py`).**

- *Mechanism.* Learn slow collective modes from trajectory data (`autoblock.py` TICA + MSM + PCCA+), then map the slowest mode's loadings onto torsions to define the reaction block.
- *Axis.* Inter-basin (defines the block that crosses).
- *Correctness.* Analysis only — does not touch `π`. Its *output* feeds C1's fixed scan, which is unbiased.
- *Fit.* Already implemented.
- *Scale.* Local per region (constraint 6).
- *Cost/pilot.* Needs trajectory data.
- *Evidence.* **[Published]** (deeptime/TICA/MSM; external). In-tree in `autoblock.py`.
- *Composability.* Feeds C1.
- **Key limitation (task-confirmed, and true):** tICA/VAMP fit to single-basin-trapped data cannot see the barrier mode — the slowest *observed* process within one basin is not the inter-basin coordinate. Therefore D1 alone cannot bootstrap crossing; it needs crossing data (D3) or endpoint contrast (D2).

**D2. Differential / contrastive analysis across known basin endpoints.**

- *Mechanism.* When both basin structures A and B are known (the usual case for a targeted transition, e.g. cis vs trans, active vs inactive GPCR), compute the *difference* of intra-basin correlation/MI matrices, or a contrastive projection (the directions that maximize A-vs-B variance ratio). This yields the concerted set *without* needing an observed crossing.
- *Axis.* Inter-basin.
- *Correctness.* Analysis only; feeds C1.
- *Fit.* `compute_differential_correlation` scaffold exists in `autoblock.py`.
- *Scale.* Local.
- *Cost/pilot.* Needs the two endpoints, not a crossing.
- *Evidence.* **[Inference]** (contrastive/differential PCA; external, e.g. Abid et al. 2018 contrastive PCA — not in `index.yaml`).
- *Composability.* Feeds C1; complements D1.

**D3. Bootstrap loop (seed endpoints → harvest crossings → refit RC → rerun).**

- *Mechanism.* Seed replicas at A and B; run the current scan (plus an unbiased crossing engine from family E); harvest the first crossings; pool the two-basin data; refit the reaction coordinate (D1/D2 on pooled data now *does* see the barrier mode); rebuild the reaction block; iterate.
- *Axis.* Inter-basin.
- *Correctness.* Each iteration's production run is unbiased (fixed union scan). The refit uses only past data → if block selection changes between whole runs (not within a run), each run is a separate unbiased chain; combine with care (do not concatenate as one Markov chain across a selection change).
- *Fit.* Matches the `for cycle in range(...)` loop already sketched in `autoblock.py`.
- *Scale.* Local per region.
- *Cost/pilot.* This *is* the pilot.
- *Evidence.* **[Inference]** (standard enhanced-sampling RC-refinement loop).
- *Composability.* Wraps C1 + E1/E2.

**D4. Committor / SPIB / diffusion-map reaction coordinates.**

- *Mechanism.* Learn a committor-based or information-bottleneck (SPIB) reaction coordinate to define the concerted set more sharply than linear tICA.
- *Axis.* Inter-basin.
- *Correctness.* Analysis only.
- *Fit.* Heavier ML dependency; deferrable.
- *Evidence.* **[Published]** external (SPIB: Wang & Tiwary 2021; committor: E & Vanden-Eijnden). Not in `index.yaml`.
- *Composability.* Optional upgrade to D1/D2.

---

### Family E — Unbiased generalized-ensemble crossing engines (composable with reduced-soft worlds)

**E1. Temperature replica exchange over the SAME reduced-soft Gibbs sweep.**

- *Mechanism.* Run `R` replicas of the identical world-scan on a geometric temperature ladder; propose neighbor swaps accepted with `min(1, exp(−(β_n−β_m)(E_i−E_j)))` (`sugita_okamoto_1999`; ladder spacing from `kofke_2002`). Hot replicas cross barriers using the *reduced torsional* representation — which stays numerically stable at high T where fully-flexible MD would blow up bonds/angles (`spiridon_2017_cdhmc_gibbs` §4 explicit future direction: "constrained dynamics could enable higher-temperature replicas... enhancing relative domain motion rather than unfolding"). Swaps carry crossings down to the cold target.
- *Axis.* Inter-basin (the primary unbiased crossing engine).
- *Correctness.* Preserves `π` exactly: the cold marginal equals target Boltzmann (`chodera_2011_gibbs_replica_exchange`; the in-flight `docs/specs/rex-stationary-distribution-oracle.md` INV-M1 gives the closed-form oracle). **Correctness condition for torsional worlds (INVARIANT, kT-convention hazard):** the Fixman term is `k_BT`-scaled (`U_F = ½ k_BT ln det M_φ`, `echenique_2006`/`go_1976`). At rung `k` the acceptance SHALL use the `T_k`-scaled Fixman and the `T_k`-scaled Jacobian so `β_k U_F` cancels the metric at every rung. Using a single-temperature `U_F` across the ladder biases every non-target rung. Separately, the swap acceptance compares configurational energies and SHALL be consistent about whether `U_F`/Jacobian are inside the swapped energy (the parent REX spec keeps Fixman out of the swap acceptance — this is correct only if the per-replica inner move already samples the `T_k`-Boltzmann-with-Fixman, so `U_F` cancels in the swap ratio; this cancellation SHALL be verified, not assumed).
- *Fit.* Excellent — Robosample has a REMC driver in flight (`docs/specs/rex-stationary-distribution-oracle.md`, `include/ReplicaExchange.hpp`); `spiridon_2020_robosample` §4 lists REX as future work; `larsen_2014_gneimo_casp_refinement` is the torsional-MD + REXMD precedent (32-replica GNEIMO REXMD).
- *Scale.* `R` replicas × cheap torsional scans; `R` grows as `sqrt(heat capacity)` ~ `sqrt(N_DOF)` for fixed acceptance (`kofke_2002`), which is the cost driver at 1M atoms — but the *reduced* worlds have far fewer active DOF than all-atom, shrinking the ladder. Explicit solvent inflates the effective heat capacity; consider H-REX (E2) instead.
- *Cost/pilot.* Ladder tuning (geometric, `kofke_2002`).
- *Evidence.* **[Published]** `sugita_okamoto_1999`, `kofke_2002`, `chodera_2011_gibbs_replica_exchange`, `larsen_2014_gneimo_casp_refinement`; in-flight REX spec.
- *Composability.* Wraps any inner scan (A*, B*, C*); the swap is orthogonal to the world design.

**E2. Hamiltonian REX with soft-core / BAT-scaling rungs (attacks the 1-4 clash directly).**

- *Mechanism.* Ladder in a Hamiltonian parameter rather than temperature: soften the 1-4 vdW / electrostatics (or scale torsional stiffness, `include/BatScaling`) at the hot rungs so the clash that gates the concerted transition is relieved, with the atomistic force field at the target rung. This is the internal-coordinate-force-field idea (`katritch_2003_icff`, `chen_2005_tamd` softened-vdW ICFF) turned into an *exact* H-REX ladder. `spiridon_2017_cdhmc_gibbs` §4: "use modified potentials as guidance while retaining the atomistic FF for acceptance"; the H-REX makes this rigorous instead of guidance-only.
- *Axis.* Inter-basin (directly lowers the barrier that gates crossing).
- *Correctness.* Preserves `π` exactly at the target rung (H-REX swap acceptance uses both replicas' energies under both Hamiltonians). Reweighting is *not* needed — the target rung is sampled directly; the soft rungs are auxiliary. The BAT-scaling driven-swap machinery already carries a Jacobian/work term (`RUN_TYPE::RENE`, `WTerm`, `docs/specs/replica-exchange-nonequilibrium-work.md`); its acceptance algebra is already tested (`tests/TestRexAcceptanceAlgebra.cpp`).
- *Fit.* Excellent — targets the exact physics `spiridon_2017` flags (1-4 clash) and reuses BAT-scaling.
- *Scale.* Fewer rungs than temperature REX because the softening is localized to the clashing interaction, not global thermal energy — better at 1M atoms / explicit solvent than E1.
- *Cost/pilot.* Choose the softening schedule and rung overlap.
- *Evidence.* **[Published]** `katritch_2003_icff`, `chen_2005_tamd` (softening), `minh_2019_algdock` (Hamiltonian REX precedent in the same lineage), `ballard_2009_rens`/`crooks_1998` (nonequilibrium swap), in-flight REX/BAT-scaling specs.
- *Composability.* Combines with E1 (2-D ladder T×λ) and C1 (soften along the reaction block).

**E3. RENS — replica exchange with nonequilibrium switches.**

- *Mechanism.* Propose swaps by a finite-time nonequilibrium switching protocol that drags each configuration toward the neighbor ensemble, accepted with the Crooks work criterion `min(1, e^{−w})` (`ballard_2009_rens`, `crooks_1998`). Increases phase-space overlap → fewer replicas. `τ=0` recovers ordinary REX; `τ→∞` gives `P_acc=1`.
- *Axis.* Inter-basin (raises swap acceptance so the ladder can be sparser / span more).
- *Correctness.* Exact (Crooks). This is essentially Robosample's *driven* REX (`RUN_TYPE::RENE`, `docs/specs/replica-exchange-nonequilibrium-work.md`).
- *Fit.* Already the design of the driven REX path (uncompiled per the REX spec's OQ-1).
- *Scale.* Switching cost per swap; amortized by fewer replicas.
- *Cost/pilot.* Protocol length `τ` tuning.
- *Evidence.* **[Published]** `ballard_2009_rens`, `crooks_1998`.
- *Composability.* A refinement of E1/E2.

**E4. NCMC — nonequilibrium candidate Monte Carlo intra-world.**

- *Mechanism.* Within a single chain (no replicas), drive a switching protocol that softens a clashing interaction or pushes along the reaction block, then accept with the NCMC work criterion (`nilmeier_2011_ncmc`). Robosample has an `add_ncmc_world` path.
- *Axis.* Inter-basin.
- *Correctness.* Exact (Nilmeier). HAZARD (documented in-tree, memory `ncmc-explicit-solvent-acceptance`): in explicit solvent the integrator shadow work scales `~dt²·N_DOF`, collapsing acceptance. In a *reduced torsional* world the propagated-DOF count is small, so NCMC acceptance is far healthier than in the all-atom Cartesian world — this is the natural home for NCMC in Robosample.
- *Fit.* Good in reduced worlds; poor in the all-atom Cartesian world.
- *Scale.* Protocol-length cost; bounded by reduced-DOF count.
- *Cost/pilot.* Protocol tuning.
- *Evidence.* **[Published]** `nilmeier_2011_ncmc`, `crooks_1998`.
- *Composability.* Combines with C1 (drive along the reaction block); alternative to E1/E2 when replicas are undesirable.

**E5. Simulated tempering / expanded ensemble.**

- *Mechanism.* A single chain wanders in temperature (or λ) with per-rung weights; unbiased with correct weights (`chodera_2011_gibbs_replica_exchange` casts as Gibbs on `(x,k)`).
- *Axis.* Inter-basin.
- *Correctness.* Exact *only* with correct free-energy weights, which require a pilot estimate. Wrong weights do not bias the *conditional* target-rung samples but destroy rung mixing.
- *Fit.* Single-chain (no replica memory) — attractive at 10k-robot scale where holding `R` full replicas is heavy.
- *Scale.* Better memory footprint than E1 at scale.
- *Cost/pilot.* Weight estimation (Wang-Landau / MBAR).
- *Evidence.* **[Published]** `chodera_2011_gibbs_replica_exchange`.
- *Composability.* Alternative to E1; combines with C1.

**E6. RXSGLD — self-guided replica exchange.**

- *Mechanism.* Ladder in a *self-guiding* temperature rather than temperature; exchange probability derived from the SGLD partition function so the base stage samples canonical (`wu_2012_rxsgld`). Higher exchange acceptance than temperature REX for large systems.
- *Axis.* Inter-basin.
- *Correctness.* Base stage canonical by construction; SGLD stages auxiliary. See F1 correctness caveat.
- *Fit.* Large-system oriented (the 1M-atom regime).
- *Evidence.* **[Published]** `wu_2012_rxsgld`.
- *Composability.* Alternative ladder for E1.

---

### Family F — Guidance-only potential softening (unbiased by L0)

**F1. SGLD guiding force as the CDHMC propagator.**

- *Mechanism.* Propagate the torsional world with a self-guided Langevin equation that boosts low-frequency (conformational) motion, then Metropolize on the true `H_acc` (L0). `spiridon_2017_cdhmc_gibbs` §4 explicitly suggests "possible use of the SGLD integrator in CDHMC moves"; `spiridon_2017` §3.5 already attributes CDHMC's efficiency to an SGLD-like mechanism (suppress high-frequency motion, focus KE on low-frequency DOF).
- *Axis.* Both (focuses transport on slow modes → intra decorrelation + inter crossing along soft modes).
- *Correctness.* Unbiased by L0 (guidance ≠ acceptance): the SGLD bias lives only in the propagator; acceptance uses true `U + U_F`. Alternatively SGLD-GLE (`wu_2016_sgld_gle`) restores exact detailed balance *directly* (no reweighting, exact NVT) and could be a stand-alone unbiased integrator — but inside CDHMC the guidance route is simplest and provably unbiased regardless of SGLD's own ensemble.
- *Fit.* Good — a propagator swap; ABA/Fixman untouched.
- *Scale.* Per-atom-tunable guiding factor is size-extensive (`wu_2016_sgld_gle`).
- *Cost/pilot.* Guiding-factor tuning.
- *Evidence.* **[Published]** `wu_2011_sgld`, `wu_2016_sgld_gle`, `spiridon_2017_cdhmc_gibbs` §3.5/§4.
- *Composability.* Combines with E*and B* (guiding force + preconditioner + REX).

**F2. Internal-coordinate force field (ICFF) as guidance.**

- *Mechanism.* Use a softened-vdW/CMAP internal-coordinate force field (`katritch_2003_icff`, `chen_2005_tamd`) as the *guidance* Hamiltonian (10–20 fs stable steps), atomistic FF in acceptance.
- *Axis.* Both.
- *Correctness.* Unbiased by L0 if atomistic `H_acc` is retained. If the ICFF is used *as the target* (no atomistic acceptance) it changes `π` → out-of-constraint unless reweighted.
- *Fit.* Good as guidance; the ICFF must be built once per system.
- *Evidence.* **[Published]** `katritch_2003_icff`, `chen_2005_tamd`.
- *Composability.* Overlaps E2 (softening); E2 makes it a rigorous ladder, F2 makes it guidance-only.

---

### Family G — Non-reversible / lifted / event-chain (research-grade, poor architecture fit)

**G1. Lifted / non-reversible MCMC (skew detailed balance).**

- *Mechanism.* Augment the state with a direction variable; use skew-detailed-balance to suppress diffusive back-tracking, giving `~N` vs `~N²` mixing on some targets.
- *Axis.* Both, but the theoretical speedup is target-specific.
- *Correctness.* Unbiased under skew/modified detailed balance, but the momentum-flip and lifting bookkeeping across a *world change* is undefined (same boundary problem as A2, worse).
- *Fit.* Poor — invasive; interacts badly with the blocked-Gibbs world composition.
- *Scale.* Unclear at 1M atoms.
- *Evidence.* **[Published]** external (Turitsyn et al. 2011; not in `index.yaml`).
- *Composability.* Conflicts with the Gibbs-scan structure.

**G2. Event-chain Monte Carlo (ECMC).**

- *Mechanism.* Rejection-free lifted moves for particle systems.
- *Axis.* Inter/intra for liquids.
- *Correctness.* Exact, but designed for pairwise/hard-core Cartesian systems, not articulated-body torsional FF with bonded terms.
- *Fit.* Poor — no clean map to BAT/ABA.
- *Evidence.* **[Published]** external (Bernard, Krauth, Wilson 2009; not in `index.yaml`).
- *Composability.* Not composable here.

---

### Family H — Biased methods (out-of-constraint unless exactly reweighted) and MTS

**H1. Metadynamics / umbrella sampling / ABF.**

- *Mechanism.* Add a history-dependent or windowed bias along a chosen CV.
- *Correctness.* Biases `π`. Admissible only with exact reweighting (WHAM/MBAR for umbrella; metadynamics reweighting), which (a) needs a pre-chosen CV = the reaction coordinate we do not yet know, and (b) inflates estimator variance. Classified **out-of-constraint** for the base method (constraint 1); the unbiased substitute is E2 (soft-core H-REX along the reaction block).
- *Evidence.* **[Published]** external.
- *Composability.* Only via full reweighting; not recommended as base.

**H2. Accelerated MD / Gaussian-accelerated MD (aMD/GaMD).**

- *Mechanism.* Boost the potential below a threshold; reweight by cumulant expansion.
- *Correctness.* Reweighting is *approximate* (cumulant truncation); variance grows with boost. Out-of-constraint as unbiased base.
- *Evidence.* **[Published]** external.
- *Composability.* Not recommended.

**H3. Multiple-timestep integration (RESPA / r-RESPA), mollified impulse.**

- *Mechanism.* Split forces by frequency; integrate fast forces on inner steps, slow forces on outer steps, to enlarge the outer step (`vangunsteren_1977` predictor-corrector + SHAKE lineage for the constraint side).
- *Axis.* Intra-basin (timestep gain).
- *Correctness.* Symplectic/reversible → unbiased as a guidance integrator (L0).
- *Fit.* WEAK here, and this is the key negative result: the reduced torsional world has *already removed* the fast bond/angle DOF (constraint 2), so the intra-world force spectrum is compressed and there is little frequency separation left to exploit. The remaining fast content (1-4 / short-range nonbonded on the retained torsions) caps the RESPA outer step at the resonance limit `≈ 2×` the fastest retained period (~5 fs), which the reduced world's large uniform step already approaches. Mollified impulse (MOLLY) pushes the resonance limit but adds complexity. The memory note `two-robot-contact-campaign` records MTS was ruled out for the clash case. **[Inference]:** MTS is dominated by B1 (preconditioning equalizes frequencies, achieving the same "no single fast mode caps the step" outcome without the RESPA resonance ceiling).
- *Evidence.* **[Published]** RESPA (external Tuckerman 1992; not in `index.yaml`); resonance limit (external Schlick/Izaguirre). `vangunsteren_1977_md_constraints` for the constraint interface.
- *Composability.* Redundant with B1; deprioritized.

---

## Prioritized shortlist (highest leverage on BOTH axes)

Ordered by leverage-per-risk. Each entry states why it helps both axes, what it composes with, and its validation oracles (tagged PRECONDITION / INVARIANT / LEMMA).

### S1. Soft-core Hamiltonian REX along the reaction block (E2), wrapping the fixed union scan (C1)

Rationale. This is the only combination that attacks *both* the potential barrier (soft rungs relieve the 1-4 clash that gates the concerted transition, `spiridon_2017_cdhmc_gibbs` §4) and the concerted-motion requirement (C1 co-blocks the torsions the barrier couples), while remaining exactly unbiased at the target rung and reusing the in-flight REX + BAT-scaling drivers. It is the rigorous version of `spiridon_2017`'s own stated future direction. Inter-basin from the ladder; intra-basin from the reduced-world large steps at every rung.

Validation.

- **PRECONDITION (per-world π-invariance).** Every world in the scan, at every rung, SHALL pass the existing single-temperature Fixman/ensemble oracle: the Shirts energy-histogram log-ratio slope test (`shirts_2012_ensemble_validation`; `spiridon_2017_cdhmc_gibbs` §3.2) equals `−(β_2−β_1)` within one standard error, and the C4 idealized-chain torsion marginal is uniform (`spiridon_2017` §3.1). Runtime guard, not a rate test.
- **INVARIANT (target-rung marginal).** The target rung's cold marginal SHALL equal the reference (long fully-flexible run or the closed-form oracle of `docs/specs/rex-stationary-distribution-oracle.md` INV-M1). The soft rungs SHALL NOT shift `π(A)/π(B)` at the target rung.
- **INVARIANT (kT-scaled Fixman across rungs).** At rung `k` the acceptance SHALL use the `T_k`-scaled `U_F = ½ k_BT_k ln det M_φ` and `T_k`-scaled Jacobian; a discriminating test SHALL corrupt this (single-`T` Fixman across the ladder) and confirm the target marginal breaks.
- **LEMMA (crossing gain).** Independent A↔B crossings per unit compute with S1 SHALL exceed the flat torsional scan at equal compute, while `ΔF_AB` converges to the *same* value (unbiasedness). Discriminating structure: the run must exercise the concerted clash together (both reaction-block torsions + the softened 1-4) — softening without co-blocking, or co-blocking without softening, SHALL each cross more slowly than the pair.
- **LEMMA (swap-rate overlap).** Neighbor swap acceptance matches the closed-form energy-overlap (`kofke_2002` eq:10), per `docs/specs/rex-stationary-distribution-oracle.md` O6.

### S2. Fixed union scan with differential reaction block (C1 + D2/D3)

Rationale. Directly removes the concerted-splitting rate bug (C2) that within-basin clustering can introduce, at zero correctness risk (systematic scan, constraint 4) and reusing the `autoblock.py` differential-correlation machinery. Inter-basin (enables concerted motion) with retained intra-basin benefit. This is the safe alternative to adaptive block selection (C4).

Validation.

- **PRECONDITION.** Each block is a reduced-soft torsional world passing the S1 PRECONDITION.
- **INVARIANT (order-invariance oracle).** Run the union scan in two different fixed orders; stationary basin populations and `ΔF_AB` SHALL agree within CI. This is the executable correctness oracle the task requests — it catches Fixman/constraint/Jacobian/momentum-resampling defects that break per-world `π`-invariance (a defect makes populations order-dependent). **[Published composition]** `chodera_2011_gibbs_replica_exchange`.
- **LEMMA (concerted-splitting discriminator — alanine cis→trans / two-torsion clash gate).** On a system whose transition is gated by a 1-4 clash requiring two torsions to move together (alanine cis→trans, or a two-torsion toy), a scan that SPLITS the two torsions into separate blocks SHALL show conditional barrier `≥` MEP barrier and near-zero crossings, while the co-blocked union scan crosses; **both SHALL converge to the same populations**. Failing to reproduce equal populations indicates a `π` bug, not a rate bug — the two must be separated (Critical vs rate).
- **LEMMA (differential vs within-basin).** The differential (D2) reaction block SHALL contain the concerted set that a within-basin correlation block misses; verify by seeding at A and measuring crossing rate under each block definition.

### S3. Constant mass-matrix preconditioner from the existing NMA Hessian (B1)

Rationale. Cheapest intra-basin lever with existing infrastructure (`include/NMA.hpp` already computes the mode spectrum); equalizes retained-torsion frequencies so the uniform HMC step is not capped by the stiffest torsion, and aligns longer moves with soft modes (inter-basin help when the barrier is soft-mode-aligned). Substitutes for finer blocking (B-note) without the concerted-splitting risk.

Validation.

- **PRECONDITION (metric self-consistency).** Momenta drawn with `M̃` and acceptance kinetic term SHALL use the same `M̃`; `U_F` SHALL remain the true-`M_φ` Fixman. Guard.
- **INVARIANT (unbiasedness under preconditioning).** The preconditioned world SHALL reproduce the un-preconditioned world's marginal: C4 uniform torsion (`spiridon_2017` §3.1) and the alanine φ/ψ PMF (`spiridon_2017` §3.3) within CI. A preconditioner that changes the marginal is a `π` bug.
- **INVARIANT (equipartition).** `⟨½ uᵀ M̃⁻¹ u⟩ = (ν/2) k_BT` for the metric actually used (`jain_2012_icmd_equipartition`). This ties the draw to the acceptance metric and catches a mismatched preconditioner.
- **LEMMA (frequency equalization).** Post-preconditioning, the retained-mode `ω_k` spread SHALL narrow (from `computeRouteBNMA` eigenvalues) and the acceptance-preserving timestep SHALL rise versus the un-preconditioned world at equal acceptance (~0.651, `spiridon_2020_robosample` §2.1.2).
- **OPEN (NMA momentum bias audit).** The existing `DistortOption::NMA` momentum bias (`uScaleFactors`, `nmaBias`) SHALL be audited: biasing the momentum draw toward the soft mode is unbiased *only* if it corresponds to a consistent mass metric in both draw and acceptance. If the current path scales drawn momenta along `û` without a matching acceptance-metric change, it biases `π`. This is listed in Open questions.

### S4. Rejection-suppressing inner kernel: LAHMC or partial-momentum-refresh (A1/A2)

Rationale. Reduces the random-walk component of the reduced-world move, converting the same trajectory length into more ballistic transport — helps both axes at near-zero cost and drops into the existing acceptance layer, orthogonal to S1–S3.

Validation.

- **INVARIANT (generalized detailed balance).** The LAHMC transition SHALL satisfy its fixed-point equation (`sohldickstein_2014_lahmc`); a unit test SHALL verify the composed transition preserves the C4 uniform torsion marginal.
- **INVARIANT (world-boundary momentum handling).** Persistent momenta (A2) SHALL be discarded/re-projected at every world change; a test that carries momenta across a factorization change SHALL fail the order-invariance oracle (S2 INVARIANT), proving the boundary rule is enforced.
- **LEMMA (IACT reduction).** Intra-basin IACT of a within-basin torsion SHALL drop versus plain Metropolis-HMC at equal compute.

### S5 (conditional). SGLD guiding-force propagator (F1)

Rationale. `spiridon_2017` already attributes CDHMC's efficiency to an SGLD-like mechanism and names SGLD-in-CDHMC as future work; as a guidance-only propagator it is unbiased by L0 and focuses transport on the slow modes that cross barriers. Promote after S1–S4 if soft-mode transport remains the bottleneck.

Validation.

- **PRECONDITION.** Acceptance uses true `U + U_F` (guidance ≠ target).
- **INVARIANT.** Marginal unchanged vs non-guided world (C4 uniform torsion, φ/ψ PMF).
- **LEMMA.** Soft-mode (tICA IC-1) autocorrelation drops without changing populations.

---

## Invariants (correctness conditions after any candidate lands)

- **INV-1 (per-world π-invariance).** Every sampled world, at every rung/temperature, leaves its intended Boltzmann-with-Fixman measure invariant. Verified by the Shirts slope test and the C4 uniform-torsion oracle. Any candidate that changes only the guidance (softening, preconditioner inertia, SGLD force, larger step) SHALL keep the true `U + U_F − ½RT ln sin²` in acceptance (L0). Critical severity if violated.
- **INV-2 (composition / order-invariance).** The composed fixed scan is `π`-invariant regardless of order; two distinct fixed orders yield equal stationary populations and `ΔF_AB` within CI. This is the load-bearing oracle for detecting per-world defects.
- **INV-3 (Fixman kT-scaling).** Where a temperature or Hamiltonian ladder is introduced, the Fixman and Jacobian terms are scaled to the rung so the metric cancels at every rung. (Guards the kT-convention hazard flagged in `CLAUDE.md`.)
- **INV-4 (world-boundary momentum).** Momenta are never carried across a robot factorization change; they are resampled or re-projected at each world boundary.
- **INV-5 (unbiased crossing).** Every accepted crossing accelerator leaves `π(A)`, `π(B)`, `ΔF_AB` unchanged; rate improves, populations do not move. Distinguishes a rate improvement (allowed) from a `π` bias (Critical).
- **INV-6 (rate ≠ correctness separation).** The concerted-splitting failure (C2) is a *rate* defect, not a `π` defect; a slow scan and a biased scan SHALL be diagnosed differently — order-invariance (INV-2) tests `π`, MFPT/crossing-count tests rate.

---

## Interface — Touch list and conventions at risk

Components that change (per promoted candidate; none change in this catalog itself):

- **Sampler / inner kernel** (`include/HMCSampler`, `RobotIntegrator`): A1 (LAHMC accept), A2 (partial refresh), A3 (NUTS, present), A4 (delayed rejection), F1 (SGLD guidance), B3 (diagonal-mass guidance).
- **Momentum draw / mass operator** (`RobotEngine::multiplyBySqrtMInv`, `include/NMA.hpp`): B1 (constant preconditioner), B2 (RMHMC, discouraged), the existing `DistortOption::NMA` path (S3 audit).
- **Acceptance Hamiltonian** (`include/World.hpp` `H = PE+KE+U_F−½RT ln sin²`): must retain true `U_F`/Jacobian under every guidance change (INV-1); must kT-scale under ladders (INV-3).
- **Scan / world scheduler** (`python/robosample/autoblock.py`, `context.chose_correlated_bonds`, `build_gibbs_blocks_from_trajectory`, `compute_differential_correlation`): C1/C3 (union & overlapping blocks), D1/D2/D3 (reaction-block RC), C4 (adaptive, discouraged).
- **REX driver** (`include/ReplicaExchange.hpp`, `src/Context.cpp`, `include/BatScaling`): E1 (temperature), E2 (soft-core H-REX), E3 (RENS/driven), E5/E6 (ST/RXSGLD).
- **NCMC world** (`add_ncmc_world`): E4.

Conventions at risk (must be checked whenever these components move):

- The Fixman term and its `kT`-scaling; the `−½RT ln sin²` torsion-measure Jacobian (`include/World.hpp`). Convention ambiguity between `U_F = ½ k_BT ln det M_φ` (`echenique_2006`) and `Vc = ½ ln det M` in kT units (`jain_1997`) — the acceptance SHALL fix one convention and the ladder SHALL scale it consistently (INV-3).
- The equipartition momentum metric (`jain_2012`) vs any preconditioner metric — must match in draw and acceptance (INV-1).
- The label-swap vs coordinate-swap REX conventions and by-thermodynamic-state recording (`docs/specs/rex-stationary-distribution-oracle.md` INV-M3/M6).
- Frame F/M mobilizer placement (`spiridon_2020_robosample` §2.1.3: random joint-frame placement degrades transition rates, Suppl. Table S3) — reaction-block worlds SHALL place mobilizer frames along bonds.

---

## Validation strategy (oracles the coder SHALL implement per promoted candidate)

Global oracles (apply to any promotion):

- **O-INV-2 (order-invariance).** [INVARIANT] Two fixed scan orders → equal basin populations / `ΔF_AB` within CI on alanine dipeptide φ/ψ (`spiridon_2017` §3.3 reference values). Discriminator: inject a per-world Fixman defect → populations become order-dependent.
- **O-C4 (idealized-chain uniform torsion).** [INVARIANT] Any guidance/preconditioner/kernel change reproduces the C4 uniform torsion marginal (`spiridon_2017` §3.1); non-uniformity ⇒ missing/incorrect `U_F` (the Forrest-Suter bug).
- **O-Shirts (energy-histogram slope).** [INVARIANT] Two-temperature log-ratio slope `= −(β_2−β_1)` within one SE (`shirts_2012_ensemble_validation`, `spiridon_2017` §3.2, expected `0.13363` for butane).
- **O-MFPT (rate, not correctness).** [LEMMA] MFPT_AB / independent-crossing-count per unit compute improves vs the flat torsional scan (`spiridon_2020` Table 4 methodology; `autoblock.py` MSM). Reported *alongside* O-INV-2 so a rate gain is never mistaken for correctness.
- **O-cis/trans (concerted-splitting discriminator).** [LEMMA] Alanine cis→trans (1-4-clash-gated, needs both torsions): split-block scan fails to cross with conditional barrier `≥` MEP; co-blocked/soft-core scan crosses; both converge to the same populations. This is the single oracle that exercises the concerted set *together* and fails on either torsion or softening alone.

Per-family additions are stated inline in each shortlist entry (S1–S5).

Estimator methodology (both axes):

- Intra-basin: IACT / effective sample size of a within-basin observable (end-to-end distance or a per-basin torsion), per unit compute (`spiridon_2020` §2.1.2 autocorrelation-time tuning).
- Inter-basin: MFPT_AB and independent A↔B crossings per unit compute (MSM/PySAL transition matrix, `autoblock.py`; `spiridon_2017` §3.4).
- Composed-chain convergence: round-level population / `ΔF_AB` vs a trusted reference (long fully-flexible run or REX cold marginal), Gelman-Rubin across independent seeds, and O-INV-2.

---

## Consequences and trade-offs

- Choosing the **fixed union scan** (C1) forecloses adaptive block selection (C4) as the base method; adaptive selection is admitted only with a diminishing-adaptation proof. This trades a modest efficiency ceiling for correctness safety.
- Choosing a **constant preconditioner** (B1) forecloses position-dependent RMHMC (B2) for the same DOF; B2's extra log-det bookkeeping is not worth it given the physical `M_φ`/Fixman already is the position-dependent metric.
- **Soft-core H-REX** (E2) is preferred over **temperature REX** (E1) in explicit solvent because the ladder length scales with the softened-interaction free energy, not the global heat capacity (`kofke_2002`) — but E2 requires choosing which interaction to soften (the reaction block), coupling it to D2/D3.
- **MTS/RESPA** (H3) is deprioritized: the reduced world already removed the fast modes, so the resonance ceiling leaves little to gain over B1.

---

## Open questions

- **OQ-1 (NMA momentum-bias correctness).** Does the existing `DistortOption::NMA` path (`uScaleFactors`/`nmaBias`, `include/NMA.hpp`, `include/World.hpp`) bias the momentum draw with a matching acceptance-metric change, or does it scale drawn momenta along the soft mode without a consistent kinetic-energy metric? The former is unbiased (a preconditioner, S3); the latter biases `π` (Critical). Adjacent invariant at stake: INV-1 for every torsional world already using NMA distortion. Needed to proceed: confirmation of whether the acceptance Hamiltonian's kinetic term uses the same anisotropic metric as the NMA-biased draw. Blocking for S3 promotion; non-blocking for S1/S2.
- **OQ-2 (REX swap Fixman cancellation).** The in-flight REX spec keeps Fixman out of the swap acceptance. For torsional worlds this is correct only if each replica's inner move already samples the `T_k`-Boltzmann-with-Fixman so `U_F` cancels in the swap ratio. Adjacent invariant: INV-3. Needed to proceed: a derivation (or the `rex-stationary-distribution-oracle` extended to a torsional, metric-varying rung) confirming `U_F` cancels in the label-swap ratio at differing temperatures. The current oracle uses pure-translational worlds precisely to avoid this (its C2/OQ), so it does not yet answer this. Blocking for S1 (soft-core H-REX over torsional worlds).
- **OQ-3 (differential correlation without crossings).** D2 assumes both basin endpoints are known structures. For discovery problems (unknown product basin) the reaction block cannot be defined by contrast and must come from the D3 bootstrap, which needs at least one crossing — a chicken-and-egg that only an unbiased engine (E1/E2/E4) breaks. Needed to proceed: confirmation of the target use case (targeted transition with known endpoints vs blind discovery). This determines whether S2 can start from D2 or must start from S1 to harvest the first crossings.
- **OQ-4 (kT convention).** The acceptance's Fixman convention (`½ k_BT ln det M_φ` vs `½ ln det M` in kT units) must be pinned before any ladder is built (INV-3). Needed: the exact convention in `include/World.hpp`/`FixmanCorrection`. Non-blocking for reading, blocking for E1/E2 implementation.

---

## References (index.yaml keys and external)

- `spiridon_2017_cdhmc_gibbs`, `spiridon_2020_robosample` — CDHMC-as-Gibbs, worlds, MFPT, 1-4-clash gating, SGLD-in-CDHMC and REX as stated future directions.
- `duane_1987_hmc`, `forrest_1994_gchmc`, `sohldickstein_2014_lahmc`, `betancourt_2016_xhmc`, `brubaker_2012_chmc` — HMC kernel and dual-Hamiltonian freedom (L0), fictitious inertia, LAHMC, XHMC/NUTS, constrained/RMHMC.
- `jain_2012_icmd_equipartition`, `jain_1997_fixman_compensating_potential`, `jain_2013_fixman_branched`, `kandel_2016_fixman_hybrid_icmd`, `echenique_2006_stiff_rigid_constraints`, `go_1976_rigid_flexible_partition`, `vitalis_2014_mixed_rigid_dihedral`, `mazur_1991_icmd`, `chen_2005_tamd`, `katritch_2003_icff`, `vaidehi_1996_neimo_thermostat`, `larsen_2014_gneimo_casp_refinement` — internal-coordinate dynamics, modal/diagonal decoupling, Fixman, ICFF softening, torsional REXMD.
- `sugita_okamoto_1999`, `kofke_2002`, `chodera_2011_gibbs_replica_exchange`, `ballard_2009_rens`, `crooks_1998`, `nilmeier_2011_ncmc`, `wu_2012_rxsgld`, `wu_2011_sgld`, `wu_2016_sgld_gle`, `minh_2019_algdock` — generalized-ensemble crossing engines and guidance softening.
- `shirts_2012_ensemble_validation` — ensemble-correctness slope oracle.
- `docs/specs/rex-stationary-distribution-oracle.md`, `docs/specs/replica-exchange-nonequilibrium-work.md` — in-flight REX drivers and oracles this catalog composes with.
- External (not in `index.yaml`, cited as inference sources): Liu 1994 / Roberts & Sahu 1997 (Gibbs blocking of correlated variables); Roberts & Rosenthal 2007, Latuszynski et al. 2013 (adaptive-MCMC correctness); van Dyk & Park 2008 (partially-collapsed Gibbs); Tuckerman 1992 / Izaguirre resonance (RESPA); Wang & Tiwary 2021 (SPIB); Abid et al. 2018 (contrastive PCA); Turitsyn 2011 (lifting), Bernard-Krauth-Wilson 2009 (ECMC).

---

## Addendum: molecules in contact (docking, membrane lipids, explicit solvent)

This section extends the catalog to the case where two or more molecules are in contact. It reuses the standing lemma **L0** (dual-Hamiltonian freedom), the mass-metric family **B**, the blocking family **C**, the reaction-coordinate family **D**, the generalized-ensemble family **E**, the guidance-softening family **F**, and the MTS analysis **H3**. Candidates are family **K**; the contact shortlist is **SK**. Evidence tiers as before: **[Published]** / **[Inference]** / **[Intuition]**.

### Problem restatement (contact regime, in the codebase's vocabulary)

Three concrete regimes, all sharing one diagnosis:

- **(a) Docking** — receptor + ligand. The ligand is a reduced-soft robot (internal torsions + a 6-DOF root); the receptor is held rigid or partly flexible. The rate-limiting motion is the ligand rotating/translating against the pocket wall.
- **(b) GPCR in a membrane** — protein + densely packed lipids. Each lipid tail is a torsional chain; tails are in permanent side-by-side contact; a tail torsion cannot sweep without pushing a neighbor tail.
- **(c) Explicit solvent** — water + ions with PME + PBC (already supported, `python/robosample/run.py` Cartesian OpenMM world). Waters carry no useful torsions; each is a Free-root rigid body pushed against its solvation shell.

### K0 — Core diagnosis: contact stiffness sits *inside* the soft sampled coordinate

In the intramolecular case, constraint (2) works because freezing bond/angle DOF removes the fast modes and leaves a soft torsional spectrum that tolerates a large HMC step (`mazur_1991_icmd`, `spiridon_2017_cdhmc_gibbs` §1). Contact defeats this. The **non-bonded** terms — the steric `r^-12` wall and electrostatics — inject high curvature *along the very torsional coordinate being sampled*: a small increment `Δξ` of a ligand rotation or a lipid tail torsion drives an atom into a neighbor's repulsive wall, so the effective potential along that soft coordinate has a locally enormous second derivative. The maximum efficient HMC timestep is set by the stiffest curvature the integrator sees over the trajectory (acceptance falls as the leapfrog energy error grows, `duane_1987_hmc`; the ~0.651 acceptance target, `spiridon_2020_robosample` §2.1.2). Under contact that stiffest curvature is a **non-bonded** effective frequency living inside the retained soft DOF, not an intramolecular fast mode. Excluding bonds/angles from the world does not touch it.

This is exactly the two-robot-contact campaign's finding, expressed there as the win metric **ℓ > Δξ**: the per-step displacement `Δξ` along the contact-sensitive coordinate SHALL stay below the contact length scale `ℓ` (the width of the wall region) or the move slams into the wall and is rejected (`tests/TestTwoRobotContact.cpp`, `docs/specs/two-robot-contact/` INV-REV/INV-KE per the MTS integrator doc and memory `two-robot-contact-campaign`). Every contact candidate below is a strategy to restore `ℓ > Δξ` without lowering the step of the *other*, non-contacting worlds. **[Inference]** from `spiridon_2017` §4 (1-4 vdW clashes gate low-T torsional transitions) generalized from intramolecular 1-4 to intermolecular contact.

### K-fidelity — The high-timestep / atomistic-fidelity settlement (extends L0)

The user asked whether large HMC timesteps push contact systems into a "hydrodynamic regime" and drop atomistic character. Settlement, defended:

1. **HMC preserves the exact atomistic Boltzmann `π` at any timestep.** Leapfrog (constrained or not) is volume-preserving and time-reversible for *any* `dt`; the Metropolis step corrects the `O(dt^p)` integrator energy error exactly (`duane_1987_hmc`; `brubaker_2012_chmc` proves detailed balance for the constrained RATTLE leapfrog; `spiridon_2017_cdhmc_gibbs` §1.3: "the integrator need not preserve Boltzmann — only the forward/reverse ratio matters"). Large `dt` costs **acceptance**, not distributional fidelity. Past the stability limit the chain stalls (acceptance → 0, coordinates stop moving) but the samples it *does* accept are still drawn from the exact `π`. There is no biased-but-running regime from `dt` alone. **[Published]**
2. **A hydrodynamic/continuum regime is a deliberate model change, not an emergent one.** It arises only if you *change the target*: implicit solvent replaces explicit waters by a solvent PMF (a different `π`, the AlGDock/implicit-ligand lineage `minh_2012_implicit_ligand_theory`, `minh_2019_algdock`), or Brownian/Stokesian dynamics replaces the Hamiltonian by overdamped Langevin. Rigid-chain Brownian dynamics is exactly `pear_1979_brownian_rigid_chain`: it still carries the metric determinant `√g` and a Fixman potential `U = kT ln√g` to recover the flexible distribution — i.e. even the continuum model must reproduce the same configurational `π` if done correctly, and it is adopted by *choice*, never induced by a timestep. **[Published]** `pear_1979`, `minh_2012_implicit_ligand_theory`.
3. **The Mori-Zwanzig picture affects kinetics, not the sampled equilibrium.** Integrating out fast/solvent DOF gives the slow torsions a friction + memory kernel (a generalized Langevin equation). That is a statement about *dynamics*; the stationary distribution of a correctly constructed GLE is still the marginal Boltzmann of the slow coordinates (`wu_2016_sgld_gle` restores exact NVT under a GLE). For a *sampler* — which needs only the equilibrium `π`, and whose torsional-world "time" is already fictitious (`spiridon_2020_robosample` §4 states rigid-body simulation cannot reproduce time-dependent quantities) — MZ friction is irrelevant to correctness. It only reminds us that contact-world MFPT is a *sampling-efficiency* metric, not a physical rate. **[Published]** `wu_2016_sgld_gle`, `spiridon_2020` §4.

Consequence for the acceptance layer: **INV-1** (per-world `π`-invariance) is unaffected by `dt`. Every contact candidate that only touches the propagator or the timestep inherits L0 and cannot bias `π`; only candidates that change the *target* (K5 GCMC ensemble choice, K7 implicit solvent) can, and those are flagged explicitly.

### Reconciliation with the two-robot-contact / MTS-ruled-out finding (REQUIRED before K1)

The two-robot-contact campaign ruled out multiple-timestep (MTS/RESPA) integration and chose a **mixed torsional/Cartesian "contact world"** instead; the Fixman correction needed no change (INV-FIX), and the win metric was `ℓ > Δξ` (memory `two-robot-contact-campaign`; `docs/specs/refactor/DOC-MTSIntegrator.md` records "MTS was ruled out for the contact-world campaign, but it exists as an option"; `tests/TestTwoRobotContact.cpp`). The reason MTS was ruled out, reconstructed and stated so K1 can be judged:

RESPA sub-cycles a *fast force* that acts on a *separable fast coordinate* — you integrate that coordinate at small `dt` while advancing the slow coordinate at large `dt`. In a reduced torsional world the contact stiffness is curvature along the **same** generalized torsion that carries the soft motion; there is no separate fast DOF to sub-cycle. Splitting the force by frequency does not decouple, and the RESPA nonlinear-resonance limit still caps the outer step at `≈ 2×` the fastest period represented in the outer force (`≈` the contact-wall oscillation projected onto the torsion). So integrator-level MTS **inherits** the K0 ceiling.

Critically, the chosen mixed contact world **already provides** the timescale separation MTS was meant to provide, but at the *Gibbs-block level*: the stiff contact direction becomes a **separate coordinate in a separate world** (the Cartesian branch of the mixed world), which can run its own small `dt` while the non-contacting torsional worlds keep their large `dt`. Different worlds already carry different timesteps (`spiridon_2020_robosample` §2.2.2 gives per-world `L×ε`). So world-level MTS is exact, avoids the resonance ceiling, and makes integrator-level MTS redundant. **This is the decisive point:** interaction-split RESPA (K1) does *not* escape the rule-out reason on its own; the escape is to first make the stiff contact direction its own coordinate (K4) or slow it with an anisotropic mass (K2) — after which the separation is already a world/block split (K6), not an integrator MTS. **[Inference]** grounded in the campaign finding and RESPA resonance theory (external Tuckerman 1992; Izaguirre resonance).

### Family K — contact candidates

**K1. RESPA / MTS split by interaction (large outer step for soft/long-range/intramolecular, inner steps for close-contact steric).**

- *Mechanism.* Force-group split the contact `r^-12`/short-range term to inner sub-steps, everything else to the outer step, to enlarge the outer step under contact.
- *Axis.* Intra-basin (timestep under contact).
- *Correctness.* Symplectic/reversible → unbiased as a guidance integrator (L0); Fixman unchanged (INV-FIX).
- *Fit.* **Poor — inherits the rule-out (see reconciliation above).** In the reduced world the contact stiffness is not a separable fast DOF, so the outer step stays capped at the contact-resonance limit; the win over the reduced world's existing large step is marginal. Escapes the ceiling only when combined with K4 (separate contact coordinate), at which point the separation is a world/block split (K6), not integrator MTS.
- *Scale.* Neutral.
- *Cost/pilot.* Force-group tuning.
- *Evidence.* **[Inference]** (campaign rule-out + RESPA resonance, external). Redundant with K6, dominated by H3's conclusion.
- *Composability.* Only meaningful atop K4/K6; do not deploy standalone.

**K2. Anisotropic mass reflecting contact curvature (contact-aware metric HMC).**

- *Mechanism.* Give the momentum draw and kinetic term a metric that is **heavy along the clash normals and the ligand-approach direction, light along the soft torsions**, so a single HMC trajectory takes small displacement along the stiff contact direction and large displacement along the soft DOF — restoring `ℓ > Δξ` on the contact coordinate only. Two variants: (i) a *constant* per-move metric `M̃` re-estimated between moves (frozen during the trajectory) — the contact analog of **B1**; (ii) a *position-dependent* metric that tracks the contact geometry as molecules move — the contact analog of **B2**.
- *Axis.* Intra-basin strongly (raises the contact-limited step); inter-basin when the crossing runs along a soft slide/rotate against a wall (docking).
- *Correctness.* Variant (i) is unbiased by L0 with the same argument as B1: a metric held constant over the trajectory has constant `det M̃` that cancels in acceptance, and the true `U_F` (from the physical `M_φ`, unchanged by contact — INV-FIX) stays in `H_acc`. **INVARIANT:** draw and acceptance SHALL use the same `M̃`; `U_F` SHALL remain the true-metric term. Variant (ii) is a genuinely position-dependent metric → non-separable Hamiltonian → requires the `½ ln det M̃(q)` bookkeeping term in acceptance *and* an implicit (generalized-leapfrog) integrator, whose cost is the reason B2 was discouraged; a contact metric that changes as molecules move is expensive to recompute and to keep symmetric-positive-definite. **[Published]** (RMHMC cost, external Girolami & Calderhead 2011; `brubaker_2012_chmc`).
- *Fit.* Variant (i) reuses `RobotEngine::multiplyBySqrtMInv` / `include/NMA.hpp` with the metric built from the *contact Hessian* (the non-bonded second derivatives along the ligand root + interface torsions), not the intramolecular Hessian — a modest extension of Route-B NMA. Good.
- *Scale.* Variant (i): one contact-Hessian build per move region, local (constraint 6). Variant (ii): poor at 1M atoms (per-step fixed-point solves).
- *Cost/pilot.* Contact-Hessian estimation (finite-difference of the non-bonded generalized force, as `computeRouteBNMA` already does for the bonded force).
- *Evidence.* **[Published]** `brubaker_2012_chmc`, `jain_2012_icmd_equipartition` (metric decoupling); **[Inference]** for the contact-Hessian construction. Relates to **B1/B2**.
- *Composability.* Variant (i) combines with K3 (soft-core removes the singular spike; anisotropic mass slows the residual approach), K4, K6, and E2. Do not combine the two variants (competing metrics, cf. B1/B2).

**K3. Soft-core proposal + exact-potential acceptance (the highest-leverage contact primitive).**

- *Mechanism.* Propagate the contact world on a **de-singularized** non-bonded potential (soft-core `r^-12`/electrostatics with a shift parameter, so the wall has finite height and finite curvature over the whole approach), then Metropolize against the **true** atomistic energy. The trajectory never sees the `r^-12` singularity, so it takes large steps without being kicked by the wall; the true energy still gates acceptance, so a move that ends inside a real clash is still rejected. This is the direct application of L0's dual-Hamiltonian freedom to the non-bonded contact term, and the rigorous form of `spiridon_2017_cdhmc_gibbs` §4's "use modified potentials as guidance while retaining the atomistic FF for acceptance." It effectively widens `ℓ` (the guidance wall is soft) while keeping the acceptance wall sharp.
- *Axis.* Both — larger contact steps (intra) and the ability to slide a ligand/tail *through* a transient soft clash to the far side (inter).
- *Correctness.* Unbiased by L0 **provided** (a) the soft-core enters only the guidance forces, not `H_acc`; (b) the leapfrog on the soft-core forces stays reversible and volume-preserving (soft-core changes force magnitudes, not the integrator's symplectic structure, so this holds); (c) the momentum draw/kinetic term are self-consistent; (d) `U_F` (true metric, unchanged) stays in `H_acc` — INV-FIX. **No reweighting is needed** because the soft-core is guidance-only, not a modified target. The forward/reverse proposal-density ratio is unity for the reversible soft-core leapfrog (the standard HMC condition, `duane_1987_hmc`, `spiridon_2017` §1.3). **INVARIANT:** the volume-preservation/reversibility bookkeeping of the soft-core leapfrog SHALL be verified exactly as the intramolecular CDHMC integrator is.
- *Fit.* Excellent — a guidance-force swap in the contact world; ABA and Fixman untouched. Reuses the alchemy/soft-core force machinery already present for BAT-scaling/alchemical worlds (`include/AlchemyForceFactory`, `docs/specs/refactor/DOC-AlchemyForceFactory.md`).
- *Scale.* Free (one extra force-group evaluation for guidance; the true energy is already computed for acceptance).
- *Cost/pilot.* Soft-core shift parameter tuning (target ~0.651 acceptance).
- *Evidence.* **[Published]** L0 (`duane_1987_hmc`, `spiridon_2017` §1.3/§4), soft-core alchemy (external Beutler et al. 1994), `katritch_2003_icff`/`chen_2005_tamd` (softened non-bonded in internal coordinates). Relates to **E2** (soft-core as a REX rung) and **F2** (soft-core as guidance).
- *Composability.* Combines with K2 (metric + soft-core), K6 (soft-core inside the contact world), and — most importantly — with **E2**: the *same* soft-core parameter that is guidance-only here can instead be a REX ladder rung, giving a continuous choice between guidance-only (K3, single chain) and exact-ladder (E2, replicas).

**K4. Contact-respecting coordinates (pocket-aligned docking DOF; intermolecular relative joint with directional stiffness).**

- *Mechanism.* Choose the sampled coordinates so the stiff and soft contact directions are *separate* generalized coordinates. For docking: sample the ligand 6-DOF in a **pocket-aligned frame** so the stiff "approach/insertion" translation is one coordinate (small step) and the soft "slide/rotate in the plane of the pocket" are others (large step). For lipids: exploit chain-torsion locality (a tail torsion's contact is dominated by a few near neighbors). Generally, promote the **intermolecular relative rigid-body DOF to a first-class sampled joint** with a per-DOF (directional) effective stiffness, so K2's anisotropy is expressed structurally rather than through a dense metric.
- *Axis.* Both (docking pose search is inter-basin; pocket-local decorrelation is intra).
- *Correctness.* Unbiased — a coordinate choice only; the Fixman term and the `−½RT ln sin²` Jacobian already handle whatever metric the new joint induces (INV-FIX, and the World acceptance already carries the joint Jacobian, `include/World.hpp`). Mobilizer frame placement matters: `spiridon_2020_robosample` §2.1.3 (Suppl. Table S3) shows off-bond joint frames degrade transition rates — the pocket-aligned frame SHALL be placed on the physical approach axis, not arbitrarily.
- *Fit.* Partial infrastructure exists: the **docking world** already does a rigid-body "KICK" of the ligand (a symmetric Cartesian proposal; `include/World.hpp` header, `docs/specs/refactor/DOC-DockingMove.md`); K4 upgrades that from an isotropic kick to a pocket-aligned anisotropic move. Robosample's Cylinder/Ball/Slider joints (`spiridon_2020` §2.1.3) already provide directional intermolecular joints.
- *Scale.* Good — local per interface.
- *Cost/pilot.* Identify the pocket/approach axis (one-time per complex; from the bound-pose geometry).
- *Evidence.* **[Published]** `spiridon_2020_robosample` §2.1.3 (joint types + frame placement), `flores_2011_rnabuilder` (mobilizer frame alignment). **[Inference]** for the pocket-aligned docking frame.
- *Composability.* The structural expression of K2; combines with K3 (soft-core along the approach coordinate) and K6.

**K5. Do not torsionally sample water/ions — GCMC for number/placement, rigid-body MC for configuration.**

- *Mechanism.* Water and ions have no useful torsions; sampling them through the torsional machinery is the source of the explicit-solvent inefficiency. Instead: (i) sample each solvent molecule's **configuration** by rigid-body HMC/MC (already done — every water is a Free-root rigid body, `python/robosample/run.py` comment); (ii) sample the **number and placement** of structural/binding-site waters and ions by **grand-canonical Monte Carlo (GCMC)** insertion/deletion in the pocket/interface region, decoupling solvent occupancy from the torsional worlds entirely.
- *Axis.* Both — GCMC fills/empties buried sites that gate ligand/loop motion (inter), and equilibrates the solvation shell without pushing waters against walls (intra).
- *Correctness.* **This candidate can change the ensemble and SHALL be flagged.** GCMC samples the grand-canonical (μVT) distribution in the sampled region. It is unbiased **iff** the excess chemical potential `μ` is calibrated to the bulk so that insertion/deletion balance reproduces bulk density (the Adams-`B` calibration). Regional GCMC coupled to the rest of the system as a reservoir leaves the whole system's configurational distribution (including a fluctuating buried-water count) at its correct equilibrium; an *uncalibrated* `μ` biases the local water count → biases `π` (Critical). If the global target is strictly fixed-`N` NVT, GCMC SHALL be used only as a *regional* semi-grand move with a reservoir, not as a global `N`-changing move. GCMC composes as its own Gibbs block ("solvent world") and, like every world, SHALL be `π`-invariant for the chosen ensemble (INV-1). **[Published]** external GCMC water theory (Adams 1975; Ross, Bodnarchuk & Essex 2015; GCNCMC Ben-Shalom et al. 2019); relates to the AlGDock implicit-water lineage (`minh_2012_implicit_ligand_theory`).
- *Fit.* Decouples solvent from ABA/Fixman entirely — the solvent world is not a torsional robot. Directly addresses the team's open structural/binding-site-water problem (memory `reaction-force-monitoring-campaign`, FFAR1).
- *Scale.* GCMC insertion acceptance falls in dense/buried regions; pair with NCMC-GCMC (E4-style switched insertion) at scale. Local per region (constraint 6).
- *Cost/pilot.* `μ` calibration against bulk (one-time per solvent model/temperature).
- *Evidence.* **[Published]** external GCMC; `minh_2012`, `minh_2019` (implicit-water alternative).
- *Composability.* A parallel Gibbs block alongside the torsional worlds; combines with K3/K6 (GCMC frees a buried site, soft-core lets the ligand slide into it). Does *not* compose with K7 implicit solvent (mutually exclusive solvent models).

**K6. Dedicated contact world/block (the multi-molecule analog of the reaction block C1).**

- *Mechanism.* Define one world whose free DOF are the **intermolecular relative rigid-body DOF plus the interface torsions of *both* partners**, run at its **own small timestep**, scanned alongside the cheap large-step intramolecular torsional worlds. This is the mixed torsional/Cartesian contact world the two-robot campaign already chose, generalized: the contact world co-blocks exactly the coordinates whose coupling is the contact (the multi-molecule version of co-blocking the concerted set in C1), so a single move can relieve a clash by moving the interface of both partners together — which separate per-partner worlds (each seeing the other frozen) cannot.
- *Axis.* Inter-basin (concerted interface motion / pose change), with intra benefit at the interface.
- *Correctness.* Unbiased — a world in the systematic scan, `π`-invariant, Fixman unchanged (INV-FIX, the campaign's confirmed result). The order-invariance oracle (**O-INV-2**) applies unchanged. This candidate carries the same **concerted-splitting hazard C2** across molecules: splitting the two partners' interface DOF into separate worlds manufactures a conditional contact barrier ≥ the true one. K6 is the fix.
- *Fit.* Excellent — it *is* the existing mixed contact world (`tests/TestTwoRobotContact.cpp`, `docs/specs/two-robot-contact/`); this catalog places it in the general framework and connects it to C1/D2/D3.
- *Scale.* The contact world's small `dt` is paid only on the interface DOF (few), while the bulk of the system stays in large-step worlds — the world-level MTS that makes integrator MTS (K1) redundant. Local per interface.
- *Cost/pilot.* Identify the interface DOF (from a contact map at the endpoints — the D2 differential analysis applied to *inter*molecular contacts).
- *Evidence.* **[Published]/[Inference]** two-robot campaign + `chodera_2011_gibbs_replica_exchange` composition. Relates to **C1** (reaction block), **D2/D3** (interface identification), **H3** (world-level vs integrator MTS).
- *Composability.* The container for K2/K3/K4/K5: anisotropic mass (K2) and soft-core guidance (K3) go *inside* the contact world; contact-respecting coordinates (K4) *are* its coordinate choice; the GCMC solvent world (K5) is a sibling block. This is the backbone of the contact shortlist.

**K7. Solvent-model choice as a deliberate lever (implicit vs explicit), not an emergent effect.**

- *Mechanism.* Implicit solvent (GB/SA, as in `spiridon_2020` FFAR1 implicit-solvent runs and `larsen_2014_gneimo_casp_refinement` GB/SA REXMD) removes the water-contact timestep problem entirely by removing explicit waters — the (c) regime collapses to a solute-only reduced-soft problem. This is legitimate but changes `π` (to a solvent-PMF target) and SHALL be justified as a modelling decision, per K-fidelity point 2, not adopted as an efficiency trick. If explicit solvent is scientifically required (structural waters, specific ion effects), the unbiased path is the explicit-solvent stack: **K6 contact world + K3 soft-core proposals + K2 preconditioner + K5 GCMC for buried water**.
- *Axis.* Both (removes the dominant contact-stiffness source in (c)).
- *Correctness.* Implicit solvent is a *different, self-consistent* `π`; unbiased *for that target*, not for the explicit one. Not a `π` bug — a target choice. Must be reported as such.
- *Fit.* Implicit path already supported (`spiridon_2020` FFAR1 implicit, `minh_2019_algdock`).
- *Scale.* Implicit is far cheaper at 1M atoms (no water); the trade is loss of explicit-water structure.
- *Cost/pilot.* None beyond model selection.
- *Evidence.* **[Published]** `spiridon_2020_robosample`, `larsen_2014_gneimo_casp_refinement`, `minh_2019_algdock`.
- *Composability.* Implicit (K7) and explicit-with-GCMC (K5) are mutually exclusive solvent choices; both compose with K6/K3/K2 on the solute side.

### Contact shortlist (SK) with oracles

Ordered by leverage-per-risk. All inherit the global oracles **O-INV-2** (order-invariance), **O-C4** (uniform-torsion `π` guard), **O-Shirts** (ensemble slope), and **O-MFPT** (rate, reported alongside `π`).

**SK1. Soft-core-guidance contact world (K3 inside K6), optionally promoted to a soft-core REX rung (E2).**
Rationale. The single highest-leverage contact primitive: it lifts the contact timestep ceiling (K0) by softening the guidance wall while the true wall still gates acceptance, is unbiased by L0 with no reweighting, reuses the alchemy/BAT-scaling machinery, and is exactly `spiridon_2017` §4's stated direction made rigorous. It attacks all three regimes (ligand-into-pocket, tail-into-tail, water-into-shell). Continuous with E2: the same softening is either guidance-only (single chain) or a ladder rung (replicas).

- **PRECONDITION (per-world π-invariance).** The contact world at full (unsoftened) acceptance SHALL pass O-C4 and O-Shirts; `U_F` retained (INV-FIX). Runtime guard.
- **INVARIANT (soft-core is guidance-only).** With soft-core in the propagator and true energy in acceptance, the contact world's marginal SHALL equal the un-softened world's marginal (O-C4 on a contact toy; INV-5). A soft-core leak into `H_acc` SHALL be caught by a marginal shift.
- **INVARIANT (reversibility/volume bookkeeping).** The soft-core leapfrog SHALL satisfy the same reversibility/volume test the intramolecular CDHMC integrator passes (`brubaker_2012_chmc` detailed-balance condition); a deliberately non-reversible soft-core step SHALL fail it.
- **LEMMA (contact timestep gain).** At fixed ~0.651 acceptance, the soft-core contact world SHALL sustain a larger `dt` (larger `Δξ`, i.e. `ℓ` effectively widened) than the hard-core contact world, with equal populations.

**SK2. Dedicated contact world with contact-respecting coordinates and anisotropic mass (K6 + K4 + K2 variant (i)).**
Rationale. The structural backbone: co-blocks the interface DOF of both partners (removing the cross-molecule concerted-splitting hazard C2), gives the stiff approach its own small-step coordinate, and shrinks `Δξ` along the clash normal via a constant contact-Hessian metric — all unbiased, all reusing existing joints/NMA. This is the generalized, framework-placed version of the two-robot campaign's chosen design.

- **PRECONDITION.** Each contact-world DOF is `π`-invariant; the constant metric `M̃` is used identically in draw and acceptance; `U_F` is the true-metric term (INV-1, INV-FIX).
- **INVARIANT (order-invariance, O-INV-2).** Two fixed scan orders with the contact world in different positions SHALL give equal `ΔF_bind` / interface populations within CI. A Fixman or metric defect in the contact world makes them order-dependent.
- **INVARIANT (frame placement).** Mobilizer frames on the physical approach/bond axes (`spiridon_2020` §2.1.3); a randomized-frame control SHALL show degraded transition rate (Suppl. Table S3) — confirming the coordinate choice, not `π`, is what improves.
- **LEMMA (docking/lipid clash discriminator — the contact analog of O-cis/trans).** On a system whose transition is gated by an intermolecular clash requiring the interface DOF of *both* partners to move together (a ligand pose flip that needs a pocket side-chain torsion + the ligand approach to move jointly; or a lipid tail rotamer flip that needs the neighbor tail to yield), a scan that puts the two partners' interface DOF in **separate** worlds SHALL show conditional contact barrier ≥ the true barrier and near-zero crossings, while the co-blocked contact world crosses — and **both SHALL converge to the same bound/unbound (or rotamer) populations**. Failing equal populations is a `π` bug (Critical), not a rate bug (INV-6). This oracle SHALL fail on either half alone (co-blocking without a small enough `dt`, or small `dt` without co-blocking), proving it exercises the concerted contact together.

**SK3. GCMC solvent block for buried/binding-site water and ions (K5), calibrated to bulk.**
Rationale. Removes the dominant explicit-solvent contact-stiffness source (water pushed against its shell) by taking water out of the torsional machinery, and fills/empties buried sites that gate solute motion — the team's open structural-water problem. Unbiased only with calibrated `μ`, so it carries a real correctness gate.

- **PRECONDITION (μ calibration).** In bulk, GCMC insertion/deletion SHALL reproduce the reference bulk water density within CI at the chosen `μ`; an uncalibrated `μ` is a runtime-rejected configuration. This is the ensemble-correctness guard.
- **INVARIANT (ensemble consistency).** With calibrated `μ`, a control region's water-count distribution and the solute configurational marginal SHALL match a long reference explicit-solvent run (INV-5). A miscalibrated `μ` SHALL shift the buried-water occupancy — the discriminating failure.
- **LEMMA (occupancy-gated crossing).** A solute transition gated by a buried water SHALL cross faster with the GCMC block active than with fixed water count, while `ΔF` between solute states converges to the same value.

**SK4 (dominated / deprioritized). Integrator-level MTS (K1) and moving-metric RMHMC (K2 variant (ii)).**
Rationale for deprioritization, stated so it is not silently dropped: K1 inherits the two-robot rule-out (no separable fast DOF in the reduced world; resonance ceiling) and is dominated by the world-level timescale split already provided by K6; K2(ii) needs an implicit integrator and per-step metric recomputation that conflict with the cheap-large-step premise (constraint 5) and are dominated by K2(i)+K3. Promote only if a contact regime is found where the stiff direction is genuinely inseparable *and* constant-metric preconditioning provably fails — an outcome none of the cited evidence predicts.

### Contact-section additions to Interface, Invariants, Open questions

- **Interface (additions to the Touch list).** Contact candidates touch: the **contact/docking world** (`include/World.hpp` docking KICK, `docs/specs/refactor/DOC-DockingMove.md`) — K4/K6; the **alchemy/soft-core force factory** (`include/AlchemyForceFactory`, `docs/specs/refactor/DOC-AlchemyForceFactory.md`) — K3/E2; the **mass operator/NMA** with a *contact* Hessian — K2; a **new GCMC solvent block** (no ABA/Fixman) — K5; the **OpenMM PME/PBC path** and solvent-model selection — K7. Conventions at risk: mobilizer frame placement on the approach axis (`spiridon_2020` §2.1.3); `U_F` is unchanged by contact (INV-FIX) and SHALL NOT be recomputed from non-bonded curvature; soft-core SHALL stay out of `H_acc`.
- **Invariants (contact additions).** **INV-FIX-C:** contact never modifies the Fixman/Jacobian terms (contact stiffness is potential curvature, not metric change) — reaffirms the two-robot INV-FIX. **INV-ENS-C:** any candidate that changes the ensemble (K5 GCMC, K7 implicit) SHALL be labelled a target change and validated against a reference of that target, never against the explicit-NVT reference.
- **Open questions (contact additions).** **OQ-5 (contact-Hessian metric cost).** Is a constant per-move contact-Hessian metric (K2(i)) cheap and stable enough at membrane/solvent density (10^5 lipids/waters), or does the interface curvature change too fast between moves? Adjacent invariant: INV-1 for the contact world. Needed: a pilot on FFAR1+nanodisc measuring metric staleness vs move length. **OQ-6 (GCMC μ transferability).** Does a bulk-calibrated `μ` remain correct inside a crowded pocket/interface where the local environment differs from bulk, or is a position-dependent excess chemical potential required? Adjacent invariant: INV-ENS-C. Needed: a buried-site occupancy check against a long explicit reference. Blocking for SK3 in tight pockets. **OQ-7 (regime of MTS validity).** Is there any contact regime with a genuinely separable fast contact DOF where K1 beats K6's world-level split? If not, K1 SHALL be formally retired. Needed: the two-robot campaign's quantitative rule-out data (referenced in memory but the `docs/specs/two-robot-contact/` files are not on this branch).

### Contact-section references (additions)

- `pear_1979_brownian_rigid_chain` — Brownian/continuum rigid-chain dynamics still carries the metric/Fixman factor; supports "hydrodynamic regime is a deliberate model change."
- `minh_2012_implicit_ligand_theory`, `minh_2019_algdock` — implicit-solvent / integrated-out-water lineage (K7) and Hamiltonian REX precedent.
- `spiridon_2020_robosample` §2.1.3 (joint types, frame placement), §2.2.2 (per-world timestep), §4 (rigid-body dynamics ≠ real kinetics) — K4, K6, K-fidelity point 3.
- `brubaker_2012_chmc`, `duane_1987_hmc`, `spiridon_2017_cdhmc_gibbs` §1.3/§4 — L0 dual-Hamiltonian basis for K3 (soft-core guidance) and the any-`dt` fidelity argument.
- `katritch_2003_icff`, `chen_2005_tamd` — softened internal-coordinate non-bonded (K3/E2).
- `larsen_2014_gneimo_casp_refinement` — implicit-solvent torsional REXMD (K7).
- External (not in `index.yaml`): Adams 1975 / Ross, Bodnarchuk & Essex 2015 / Ben-Shalom et al. 2019 (GCMC and GCNCMC water, K5); Beutler et al. 1994 (soft-core, K3); Girolami & Calderhead 2011 (RMHMC moving-metric cost, K2(ii)); Tuckerman 1992 and Izaguirre-Reich-Skeel resonance (RESPA ceiling, K1); Zwanzig 1961 / Mori 1965 (Mori-Zwanzig, K-fidelity point 3).
- `docs/specs/two-robot-contact/` (not on the `disasm` branch; recovered via `tests/TestTwoRobotContact.cpp`, `docs/specs/refactor/DOC-MTSIntegrator.md`, `DOC-DockingMove.md`, and memory `two-robot-contact-campaign`) — the MTS rule-out, mixed contact world, INV-FIX, and `ℓ > Δξ` metric this section builds on.
