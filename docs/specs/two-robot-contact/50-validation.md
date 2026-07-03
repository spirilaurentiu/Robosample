# Two-robot contact world — Validation plan and oracles

Consolidated checkable oracles. Terms from `00-...` §2. Each tagged PRECONDITION
(runtime guard) / INVARIANT (must hold; test can fail) / LEMMA (derived value +
tolerance). The discriminating question for every oracle: what distinguishes a
correct implementation from a plausible-but-biased one.

## 1. Preconditions (runtime guards, not tests)
- **P1 (atom-disjoint partition + welded E, `10-...`).** R's atoms and E's atoms
  disjoint; E flagged `cartSolventAtoms` and assigned only to **0-DOF `Weld`-rooted
  all-`Rigid` bodies**; no R generalized coordinate moves an E atom.
  `setCartesianSolvent` SHALL `throw` if a flagged atom's body has `bodyNU > 0` or a
  flagged body is partially masked (enforced, not assumed).
- **P2 (momentum flip under partial refresh, `10-...`).** Under partial
  refreshment, reject-flip negates BOTH `u` and `v_E`; corrector non-convergence
  is a reject/reduce-`dt` condition, never silently taken.
- **`ξ` selection.** The slow coordinate(s) and their `u_ξ` projection are declared
  per system (input, not inferred).

## 2. Correctness invariants (the mixed integrator, `10-...`)
- **INV-WELD (E is 0-DOF; guards the whole M1/M2/M3 chain, PRIMARY).** Build a
  two-robot world, flag E as `cartSolvent`, and assert every flagged atom's body has
  `bodyNU == 0` and flagged bodies are fully covered by the mask.
  `setCartesianSolvent` SHALL `throw` on a `Free`-rooted or partially-masked E (P1).
  Wrong-arm: a deliberately `Free`-rooted E SHALL be rejected at construction — if
  instead allowed, KE is double-counted and the sampled density is not (T).
  Discriminates the single silent-corruption path the spec's own INV-KE does NOT
  catch (INV-KE passes trivially for correctly-welded E).
- **INV-FIX (Fixman mixed-manifold, PRIMARY correctness for question 3).** For a
  fixed E, `ΔU_F` over any `(φ,x_E)`→`(φ',x_E')` displacement is bitwise identical
  whether or not E's atoms are in the Cartesian reference `lnDetMCartesian_`.
  Discriminates: a build that mistakenly adds a *configuration-dependent* E term
  (a spurious cross-Fixman) fails; the correct build (constant `−½RT ln|M_E|`
  offset only) passes. Must exercise a move of BOTH `φ` and `x_E`.
- **INV-DRAW (momentum block-diagonality, `10-...` M3).** Over many draws, the
  sampled covariance of `(u, v_E)` matches `RT·diag(M_φ⁻¹, M_E⁻¹)` within
  statistics, with cross-block covariance ~0. Discriminates: a build that couples
  the blocks (or draws `v_E` with the wrong `σ`) fails.
- **INV-REV (joint reversibility, `10-...` M4).** `checkReversibility` on the
  JOINT R+E map at the two-robot working `dt`: round-trip residual `≤ 1e-6`. Must
  be run with E present (not R alone). Discriminates a shared-force coupling
  asymmetry.
- **INV-KE (Hamiltonian assembly).** `ke + keSolvent = ½uᵀM_φu + ½Σ_E m|v|²`
  matches `½pᵀM_mixed⁻¹p` for the drawn momenta (consistency of `H` with the
  draw). Discriminates a metric mismatch between draw and acceptance.

## 3. The reference-free ensemble gate (does it sample (T)?)
- **INV-ENS (Shirts energy-histogram test, adapted).** Run the mixed move on a
  small two-robot GB/OBC system at two temperatures; the log-ratio of PE
  histograms SHALL be linear in `U` with slope `−(β₂−β₁)` (Spiridon eq:5-6,
  `references/papers/spiridon_2017_cdhmc_gibbs`). Discriminates any Boltzmann-bias
  from a wrong Fixman/measure without an external free-energy reference.
- **INV-MARGINAL (cross-representation equality — the decisive (T) check).** The
  target robot's torsional marginal (φ/ψ or `ξ` histogram) from the MIXED move
  SHALL agree, within statistics, with the SAME robot's marginal from a
  fully-Cartesian HMC of the whole two-robot system (both target `e^{−βU}`, T).
  Discriminates: a build with the spurious cross-Fixman (INV-FIX would also catch
  it) or a wrong environment measure produces a *displaced* interface marginal.
  Include a deliberately-wrong arm (E drawn with `σ=√RT` ignoring `1/m`, or Fixman
  omitted) that SHALL diverge.

## 4. Transport oracles (does it win? `20-...`)
- **LEMMA T-WIN (question 1, make-or-break).** At the collapsed contact `dt` and the
  WORKING trajectory length, measure `ℓ_traj = v_th·min(L·dt, τ_p, t_Uturn)`,
  `Δξ_barrier`, and Cartesian's OWN diffusion constant `D_cart`. The method is
  predicted to beat Cartesian MD/HMC per force-eval when `ℓ_traj > Δξ_barrier`
  (T3, revised). Report `ℓ_traj/Δξ_barrier` (`>1` is the go/no-go) AND the measured
  speedup against `D_cart` — not against the mixed move's own `D`. NOTE: predictive,
  not a pass/fail on the code — it decides whether the method is worth shipping.
- **LEMMA T-MSAD.** `⟨Δξ²⟩ ~ t^α` with `α→2` (ballistic) for the mixed move on a
  clean single NVE trajectory, degrading toward `α→1` as `L` exceeds the Lyapunov
  time (T2). A mixed move that is `α≈1` even at short `L` and small `dt` is
  transport-dead (or defect-dominated — cross-check `30-...`).
- **INV-EFFIC (the actual payoff, must exercise a real barrier).** An
  alanine-dipeptide-style transition gated by a 1–4 clash requiring angle+torsion
  TOGETHER (C7eq↔C7ax / cis↔trans): the mixed move's transition rate SHALL exceed
  Cartesian MD/HMC at matched wall-clock (force-eval count), while the recovered
  basin free-energy difference matches within statistics. Discriminates a move
  that "accepts" by distorting the internal landscape from one that genuinely
  crosses. (Mirrors `docs/specs/ncmc-explicit-solvent/40-...` INV3.)

## 5. Friction-attribution oracles (`30-...`)
- **LEMMA DR1 (NVE energy drift).** Refresh OFF: `|⟨H(t)−H(0)⟩|` bounded and
  `O(dt²)`, no secular term; drift rate `→0` as `dt→0` at fixed physical time.
  A surviving secular, root-localized term matching the Simbody differential
  oracle is γ-pump (bug), not γ-phys.
- **LEMMA DR2 (KE-flow localization).** Per-body KE source integrating `∝L` is a
  defect; elastic exchange integrates to ~0 by recurrence. Attribution rule in
  `30-...` §3.

## 6. Pass condition
The mixed contact world is accepted as CORRECT when it passes INV-WELD, INV-FIX,
INV-DRAW, INV-REV, INV-KE (integrator), INV-ENS and INV-MARGINAL with the wrong-arm
divergence (samples (T)), and DR1/DR2 attribute any residual friction (physical
vs numerical). It is accepted as WORTH SHIPPING when T-WIN reports `ℓ/Δξ_barrier
> 1` and INV-EFFIC shows a wall-clock speedup over Cartesian MD/HMC on the
clash-gated transition. Correctness (samples (T)) and viability (beats Cartesian)
are separate gates; a build can be correct yet not worth shipping if `ℓ ≤
Δξ_barrier`.
