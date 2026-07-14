# Robosample research program: what to solve, in what order, what to publish

Status: draft synthesis from a Stage-0 orientation pass (2026-07-14). Literature is discovery-tier
(scout-level), not yet verified by `researcher`/`verifier`. Treat every citation below as a lead to
confirm, not an established fact.

## Organizing principle

The program has been anchored to the hardest unsolved problem (efficient concurrent explicit-solvent
robotics) and to one method (static Gibbs blocks) that is stuck. Neither should gate publication.
Separate two axes:

- **Next paper** — what is differentiated, ready, and low algorithmic risk.
- **Hardest problem** — the multi-quarter research bet.

The two are not the same thread, and conflating them is why the program feels stalled.

---

## Thread 1 — Gibbs block selection (the founding goal)

**Diagnosis.** Your blocking builds an **equal-time** correlation/covariance matrix from a within-basin
trajectory and partitions coordinates by it. By construction that objective encodes only co-fluctuation
*inside* the basin that generated the data. It has no mechanism to see which coordinates move together
*during* a transition, because transitions are absent from that trajectory. "The matrix is non-stationary
across basins" is therefore not a defect to repair — it is the expected behavior of the wrong objective.

**Provenance gap.** "Islam & Venugopal" could not be located in scite or on the web under any query
variant. Closest published method of that class: *Correlation-Based Feature Selection to Identify
Functional Dynamics in Proteins*, JCTC 2022, `10.1021/acs.jctc.2c00337` (instantaneous-covariance
coordinate grouping). ACTION: supply the real citation you have been working from, or accept that
surrogate.

**Reframes (each with a published anchor, none currently in `references/`):**

1. **Time-lagged objective (direct fix).** Replace equal-time covariance with time-lagged C(τ):
   tICA / VAMP / variational approach to conformational dynamics. Slow eigenvectors of C(τ) provably
   align with the slowest relaxations — the barrier-crossing coordinates. Anchors: VAMPnets
   `10.1038/s41467-017-02388-1`; Noé & Nüske VAC `10.1007/s00332-019-09567-y`. Extends your existing
   NMA slow-mode plumbing (`DOC-NMA`, `DOC-VelocityDistortion`) and self-guided Langevin
   (`wu_2011/2012/2016`) rather than replacing it.
2. **Adaptive Gibbs blocking (honest answer to non-stationarity).** Let the block partition depend on
   the current state. Roberts & Rosenthal adaptive-Gibbs theory (`10.1214/11-aap806`; companion
   arXiv:1801.09299) gives conditions (diminishing adaptation, containment) under which this preserves
   ergodicity. Turns the biology problem into a solved MCMC-theory problem with convergence proofs to
   check against. NOTE: verify whether the current fixed-per-sweep factorization interface
   (`RobotModel`/`Context`) can accommodate a state-dependent partition without a large interface change.
3. **Transition-relevance objective.** Choose blocks to maximize the Markov-model spectral gap or
   committor discrimination, not equilibrium correlation. Anchors: SGOOP `10.1073/pnas.1600917113`,
   multi-CV SGOOP `10.1063/1.5064856`.

**Cross-link.** `chodera_2011_gibbs_replica_exchange` — replica exchange *is* a Gibbs sampler, and your
REMC driver is mature. REX across basins is your existing multi-basin mechanism; connect it to the
block-non-stationarity story rather than treating it as separate.

**Publishable outputs:**

- **Paper 2 (characterization).** "Optimal within-basin coordinate blocking does not transfer across the
  conformational landscape." Map how the optimal grouping varies across folds; the CATH NMR subset is the
  right instrument because it spans folds. This is the honest negative result and motivates Paper 3.
- **Paper 3 (method).** Time-lagged and/or adaptive block selection for rare transitions. De-risked by
  Paper 2's map.

---

## Thread 2 — GPCR per-body spatial forces, holo vs apo (the differentiated observable)

**Why it is the strongest near-term paper.** Your per-body wrench (net force + net torque about a body
frame) is the *generalized reaction force conjugate to your chosen rigid-body factorization*. It exists
only because the ABA recursion already runs for HMC, so it is free, and it reports **torque about a body
frame** — a dimension the closest prior art does not naturally expose. It sits between Force Distribution
Analysis (finer, force-only, no torque) and elastic-network/PRS methods (coarser, model-derived). Much
of the plumbing exists: `DOC-ForceReducer` (per-atom→per-body reduction), `DOC-ReactionReporter`
(`captureReactionSnapshot`, `TestReactionForces`), and the Vaidehi GNEIMO/NEIMO GPCR lineage
(`vaidehi_2015_icmd_gneimo`, `vaidehi_1996_neimo`) are already in the repo.

**Prior art to confront:** Force Distribution Analysis (Stacklies & Gräter, `10.1186/1471-2105-12-101`,
`10.1371/journal.pcbi.1000306`) and Kalescky–Tao dynamic-allostery-from-interatomic-forces
(`10.1371/journal.pcbi.1004358`, `10.1371/journal.pcbi.1000574`). These already claim ligand-binding
allostery signatures from forces. Precedent, not a blocker — but the paper must benchmark against them.

**Three requirements before it is reviewer-defensible (all tractable, none a research risk):**

1. **Fix the torque reference frame.** Torque about a body is origin-dependent; adopt one physically
   motivated convention (helix centroid, Cα axis, or joint frame) before any holo/apo comparison. A
   reviewer asks this first. → `researcher` question.
2. **Benchmark against FDA/PRS on one shared system**, and show the rigid-body torque captures a
   rotational/lever-arm effect at helix scale that pairwise atom forces do not expose. That is the
   incremental-contribution argument.
3. **Show convergence of net torque.** It nearly cancels for near-rigid domains, so it is tail-sensitive;
   demonstrate statistical convergence over finite trajectories.

**Publishable output: Paper 1.** GPCR holo-vs-apo per-helix wrench signatures over your database. Ships
first and independently of every other thread.

---

## Thread 3 — Explicit solvent + robotics (the hard research bet)

**Diagnosis, refined.** There are two *separable* acceptance killers, and the current diagnosis may
conflate them:

1. **Integrator shadow work** ~ dt²·nDOF (the NCMC-specific mechanism you already localized).
2. **Steric / excluded-volume clash** — a perfect symplectic integrator with zero shadow work would
   still reject because a concerted torsional move collides with frozen solvent (cavity mismatch).

Distinguishing these is the first research question, not an implementation detail.

**Field consensus (important).** No production system was found that runs native BAT/articulated-body HMC
concurrently with unmodified explicit solvent at acceptable acceptance. Every production precedent either
uses implicit solvent for the search phase, or Cartesian/rigid-receptor MD. The implicit-search /
explicit-rescore posture is likely the field's actual consensus, not a fallback.

**Solution directions (tradeoffs):**

1. **Implicit-solvent search, explicit-solvent rescore** (MM/GBSA, MM/PBSA; PlaceWaters, WScore docking).
   Defensible as *scoring* methodology; weaker as *sampling* (GB/PB can rank poses wrong vs explicit free
   energy). Lowest risk, publishable near-term.
2. **NCMC with configurational freezing of nearby solvent** (Sindhikara et al. `10.1021/ct500340b`) +
   GCNCMC water handling. Closest fit to your existing NCMC skeleton; directly budgets the shadow-work
   term. Adds protocol-tuning and work-reweighting machinery.
3. **Mixed-resolution hybrid MC/MD with a solvent-only relaxation sub-step** between torsional proposal
   and acceptance (Ribeiro et al. `10.1002/jcc.22925`). Most mechanistically direct fix for the clash
   mechanism; requires a new Gibbs-block type and detailed-balance validation.

**Publishable output: Paper 4.** Implicit-search / explicit-rescore internal-coordinate docking & free
energy surfaces — the tractable slice. Full concurrent explicit-solvent robotics stays a longer-horizon
research track behind it.

---

## Thread 4 — Loop modeling (HCV E2) and applications

A consumer of a working sampler, not a source of method novelty. Becomes a strong application paper after
Thread 1's method matures, and a good stress-test dataset in the meantime. Do not lead with it.

---

## Decisions taken (2026-07-14)

- **Lead: both in parallel.** Paper 1 (GPCR spatial forces, Thread 2) and the Gibbs-block papers
  (Thread 1) advance simultaneously.
- **Docking posture: chase concurrency now.** Attack native articulated-body HMC in explicit solvent
  directly (Thread 3) rather than shipping the implicit-search / explicit-rescore posture. This is the
  high-risk choice with no near-term docking paper; it promotes Thread 3 to an active track and shelves
  the Paper 4 (implicit-rescore) form for now. The rescore posture stays on file as the fallback if the
  concurrency attack stalls.

## Recommended sequence

| # | Paper | Thread | Why here | Main gate |
|---|-------|--------|----------|-----------|
| 1 | GPCR holo-vs-apo per-body spatial wrench | 2 | Differentiated, plumbing exists, no algorithmic risk | Torque-frame convention; FDA/PRS benchmark; convergence |
| 2 | Non-stationarity of optimal blocking across folds (CATH NMR) | 1 | Honest result you already have data for; motivates Paper 3 | Confirm/replace the Islam–Venugopal citation |
| 3 | Time-lagged / adaptive Gibbs block selection for rare transitions | 1 | The founding-goal method paper, de-risked by Paper 2 | tICA/VAMP + adaptive-Gibbs interface feasibility |
| 4 | Implicit-search / explicit-rescore internal-coordinate docking & FES | 3 | The tractable docking slice; field-consensus posture | Separate shadow-work vs clash; pick a solvent scheme |
| — | HCV E2 loop modeling; full concurrent explicit-solvent robotics | 4, 3 | Applications and the long research bet | Downstream of Papers 3–4 |

Narrative arc: Papers 2→3 are the core methods contribution (the founding goal). Paper 1 is the
differentiated observable that ships first and independently. Paper 4 turns the sampler toward
applications.

## Immediate next actions

1. Supply the real "Islam & Venugopal" citation (or accept the JCTC-2022 surrogate). Blocks Paper 2's
   framing.
2. Decide the lead paper (see open decisions).
3. On the chosen lead, commission a `researcher` spec: for Paper 1, the torque-reference-frame convention
   and FDA/PRS benchmark design; for Paper 3, the time-lagged-vs-adaptive objective and the
   `RobotModel`/`Context` interface question.
4. Ingest the new literature branch — none of the tICA/VAMP/SGOOP/adaptive-Gibbs/FDA/NCMC-freezing
   sources above are in `references/` yet. Run `/ingest-papers` once PDFs are gathered.

## Open decisions (need the user)

- **Lead paper:** GPCR forces (ready, differentiated) vs the Gibbs method (the founding goal). Default
  recommendation: GPCR forces first, Gibbs papers in parallel behind it.
- **Docking posture:** accept implicit-search / explicit-rescore as the near-term publishable methodology,
  or hold docking until concurrent explicit-solvent robotics is solved.
