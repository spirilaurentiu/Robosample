# NCMC explicit-solvent acceptance — Diagnosis and scaling

## 1. Problem restatement

User phrasing: "NCMC acceptance is 0 for 2ala in TIP3P explicit water with `add_ncmc_world(..., ncmc_steps=20, hold_fraction=0.1, use_fixman=True, relax_solvent=True)`; why, and how do we keep acceptance usable as the moving region and the bath grow?"

Codebase restatement. `World::ncmcMove()` (`src/World.cpp:1892`) is one HMC proposal whose proposal map is the entire alchemical `λ:1→0→1` switching trajectory, accepting on the endpoint Hamiltonian difference

```
ΔH = H_end − H_start
H   = V_λ + K(u;q) + K_s(v_s) + U_F(q) − ½·RT·ln sin²γ₂ − nmaCorr
```

(`reinitialize()` `:1592`, `currentTotalEnergy()` `:1761`, acceptance `:2004–2045`). The switch is a composition of `ncmc_steps` unadjusted fixed-λ velocity-Verlet steps (`RobotEngine::verletStep`/`stepTo`, `include/RobotIntegrator.hpp`) each advancing (i) solute internal DOF via the O(n) articulated-body integrator + implicit-trapezoid velocity corrector, and (ii) — under `relax_solvent` — every welded water atom via flat-Cartesian velocity-Verlet on the same OpenMM force (`verletStep` lines 121–272; `setCartesianSolvent` `:1543`; `drawSolventVelocities` `:1564`).

The reformulated question is an **efficiency** question about the endpoint-ΔH acceptance, not a correctness question:

> `min(1,e^{−βΔH})` sees the integrator's accumulated shadow work `Q` of every propagated DOF, because propagation is unadjusted and there is no intermediate correction. In explicit solvent the propagated-DOF count is `N_int + 3·N_water`, and `Q` grows with it until `e^{−βΔH}→0`.

Two readings, treated separately because they need different fixes:
- **(R-effic)** keep the exact endpoint-ΔH proposal, drive `Q` down at the source (integrator quality, `dt`, region size). Bounded — cannot beat the bath-size scaling (§4).
- **(R-construct)** change the acceptance *construction* so bath shadow work never enters acceptance (Metropolized-dynamics NCMC). The route that scales to FFAR1 + explicit solvent.

## 2. Binding vocabulary (defined once, reused across files)

- **Region A** — atom block whose *intermolecular* nonbonded coupling is λ-scaled (`configureNcmc`, `enableAlchemy`, `createAlchemyDecouplingForces`). Currently the whole solute.
- **Bath B** — everything not in A (here all TIP3P water).
- **Protocol work `W`** — `Σ` over fixed-`q` λ-perturbation substeps of `V_{λnew}(q) − V_{λold}(q)` (Nilmeier et al. 2011 "generalized work"; the `work` accumulator in `ncmcMove`; openmmtools `_add_alchemical_perturbation_step`). Zero-coordinate-motion ⇒ unit coordinate Jacobian (`TestNCMCWork` Part 3 proves `ln|M(q)|` and the pitch term are bitwise invariant across a perturbation).
- **Shadow work / heat `Q`** — fixed-λ propagation's Hamiltonian contribution, `Σ` over propagate substeps of `H_λ(after) − H_λ(before)`. For exact continuous `π_λ`-dynamics `Q` is exchanged with a heat bath and does not bias sampling; for a discrete unadjusted integrator `Q≠0` is discretization error. `TestNCMCWork.WorkPlusHeatEqualsDeltaH` establishes `ΔH = W + Q` (for `H=V+K`).
- **`U_F(q)`** — Fixman compensating potential, `calcFixman()` `:835` (Spiridon & Minh 2017, `references/papers/spiridon_2017_cdhmc_gibbs`, Eq. 3). Pure function of `q`.
- **`gap = W − ΔH`** — the logged diagnostic (`ncmcMove:2042`). `gap = −(Q + ΔU_F + ΔJ)`, with `ΔU_F = U_F(q_end)−U_F(q_start)`, `ΔJ = −½RT·[ln sin²γ₂(q_end)−ln sin²γ₂(q_start)]`. The gap folds shadow work together with the Fixman/Jacobian state-function drift.

## 3. Diagnosis

### 3.1 The `run.py` half-diagnosis is correct as far as it goes
`run.py:119–128`: `dH(5fs)=+70, dH(2fs)=+15, dH(1fs)=+2.5 kJ/mol → acc 0%,0%,~55%`, and `ncmc_steps 20→40 @1fs → acc 0`. Both consistent with "acceptance set by integrator shadow work" and both reproduced by §4. `relax_solvent=True` (`context.py:1010`) is the right *physical* idea (the frozen cage cannot accommodate the moved solute) but cannot fix acceptance under endpoint-ΔH, because each water-relaxation step is itself a shadow-work source landing in the same `ΔH`. Central conflict: **more reversible switching + more solvent relaxation lowers reorganization work but raises shadow work, and endpoint-ΔH charges for both.**

### 3.2 `ncmc_steps 20→40 kills acceptance` is a *drift* signature
**CLAIM D1.** The propagation has a secular, sign-definite energy drift; that drift, not reorganization, is the proximate cause of 0 acceptance in the reported config. For an *exactly* symplectic integrator the shadow Hamiltonian is conserved, so `⟨ΔH⟩` at fixed `dt` is independent of trajectory length — adding steps would not raise `⟨ΔH⟩` and (smoother schedule) tends to lower reorganization. Monotone rise of rejection with `ncmc_steps` at fixed `dt` ⇒ non-symplectic propagator with per-step heat `Q_step>0` accumulating as `Q≈n·⟨Q_step⟩`. Candidate sources, attributed by the reproducer:
- **(S-b)** internal velocity corrector not reaching its fixed point — `verletStep` takes the step anyway and warns (`RobotIntegrator.hpp:396`). Non-converged velocity ⇒ non-symplectic, injects heat. **Fixable** (tighten corrector / lower `dt`): see `20-...`.
- **(S-a)** near-symplectic shadow of a *convergent* corrector. The internal integrator is reversible + volume-preserving (`checkReversibility`) but only approximately symplectic (Newmark-type position update using `qddot`, not clean leapfrog); reversible-but-not-symplectic can drift. Reducible only by `dt`, better integrator, or Metropolization.
- **(S-c)** the flat-Cartesian solvent Verlet is symplectic (its own energy bounded) but exchanges energy with the solute every step; endpoint `ΔH` sees the total and the solute-side drift is fed by `3·N_water` DOF.

NOTE. D1 does not contradict the `~dt²` reading: a single logged `dH` is one realization whose *magnitude* is the fluctuation `~√N·dt²` (§4), while the *mean* that sets acceptance is `⟨ΔH⟩`. The `ncmc_steps` dependence separates bounded (`⟨ΔH⟩` flat in `n`) from drifting (`⟨ΔH⟩∝n`).

## 4. Scaling law

Let `N = N_int + 3·N_water`, `h = dt`, `n = ncmc_steps`.

### 4.1 Bounded-symplectic floor
Order-2 symplectic: `ΔH = −h²·[H₂(end) − H₂(start)]`, `H₂` a sum of per-DOF `O(1)` terms. Weakly correlated endpoints ⇒
- per-realization magnitude `|ΔH| ~ √N·h²` (matches the `run.py` "~dt²"),
- ensemble mean `⟨ΔH⟩ = (β/2)·Var(ΔH) ~ N·h⁴` (from `⟨e^{−βΔH}⟩=1` ⇒ `⟨ΔH⟩≈½β·Var`).

Acceptance `~ erfc(√(⟨ΔH⟩/2))`. Fixed acceptance as bath grows needs `N·h⁴=O(1)`, i.e. `h∝N^{−1/4}` (Beskos–Pillai–Roberts–Sanz-Serna–Stuart, *Bernoulli* 19(5A):1501, 2013; Neal, *Handbook of MCMC* ch.5, 2011). For FFAR1 + water (`N≈10⁵–10⁶`) this alone makes endpoint-ΔH untenable.

### 4.2 Drifting regime (the reported config)
`⟨ΔH⟩ ≈ n·⟨Q_step⟩ ~ n·N·h²·(non-symplectic prefactor)` — linear in `n` and `N`; strictly worse than §4.1; matches `20→40→0`.

### 4.3 Consequence
**CLAIM D2.** In both regimes endpoint-ΔH acceptance is extensive in bath size. No `dt`/`ncmc_steps`/`mass_scale` choice removes the `N`-dependence; they move the prefactor. `relax_solvent` under endpoint-ΔH is self-limiting; the construction cannot scale to a receptor in explicit solvent. The scalable fix removes bath shadow work from the *construction*, not by tuning `dt`.

## 5. Assumptions / invariants
- **(A1)** endpoint-ΔH is unbiased: deterministic, volume-preserving, F-reversible map, palindromic schedule pinned to `λ=1` at both ends (`NcmcProtocol.hpp`; `TestNCMCWork` Parts 1&4). This spec does not modify that; any construction change SHALL preserve exactness under its own stated invariant.
- **(A2)** `U_F(q)` and the pitch term are pure functions of `q`, invariant under a fixed-`q` perturbation (`TestNCMCWork` Part 3).
- **(A3)** welded solvent contributes `det=1` to `|M(q)|`, so `U_F` is bath-count-independent (`TestNCMCWork` Part 5). The Metropolized construction keeps solvent welded (internal DOF frozen), adding only Cartesian relaxation, so the Fixman marginal is untouched.

## 6. Open provenance gaps (not blocking)
- NCMC paper absent from `references/index.yaml`. Cited by DOI: Nilmeier, Crooks, Minh, Chodera, *PNAS* 108(45):E1009, 2011, `10.1073/pnas.1106094108`. RECOMMEND adding a `meta.yaml`.
- `docs/specs/ncmc_solvent_relax.md` (referenced by `run.py:142`, `context.py:1002`, `RobotIntegrator.hpp:125`) does not exist. Reconciliation here is against the *implementation* of `relax_solvent`: a right-idea/wrong-acceptance-frame partial fix — necessary under the Metropolized construction, self-defeating under endpoint-ΔH.
