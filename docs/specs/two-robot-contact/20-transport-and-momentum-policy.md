# Two-robot contact world — Transport observable and momentum policy

Answers questions 1 and 2. Terms from `00-...` §2. This file decides whether the
method *can* win (§2-3) and how to run it if it can (§4).

## 1. The make-or-break framing

At a stiff contact both representations share the same collapsed `dt` and pay ~1
force eval per `dt` (`00-...` §1). So the mixed move wins ONLY if torsional
transport decorrelates `ξ` in less PHYSICAL time than Cartesian diffusion — the
number of force evals is (physical time)/`dt` for both. The decision reduces to a
single measurable inequality (CLAIM T3 below).

## 2. The discriminating observable (question 1)

Define, along a SINGLE NVE trajectory (momentum persisted, no refresh) and
averaged over trajectories:

- **Generalized-momentum autocorrelation of the slow coordinate**
  `C(t) = ⟨u_ξ(t₀)·u_ξ(t₀+t)⟩ / ⟨u_ξ²⟩`, with `u_ξ` a component/projection of
  `state_.u()` (engine state, logged per step). Ballistic ⇒ `C` stays near +1 /
  oscillates coherently over `L·dt`; diffusive ⇒ `C` decays exponentially with
  short `τ_p`.
- **Momentum correlation time** `τ_p = ∫₀^{L·dt} C(t) dt`.
- **Mean-square angular displacement** `⟨Δξ(t)²⟩ = ⟨(ξ(t)−ξ(0))²⟩`, `ξ` the
  dihedral from atom positions (engine state). Fit exponent `α`:
  `⟨Δξ²⟩ ~ t^α`; `α≈2` ballistic (Green–Kubo inertial regime), `α≈1` diffusive.
- **Thermal generalized velocity** `v_th = √⟨u_ξ²⟩ = √(RT·(M_φ⁻¹)_{ξξ})`
  (equipartition, Jain et al. 2012 `references/papers/jain_2012_icmd_equipartition`
  eq:26 `⟨v v*⟩ = kT I`).
- **Ballistic mean free path** `ℓ = v_th · τ_p` (angular distance traversed
  coherently before the momentum decorrelates).

Transport-theory basis: velocity autocorrelation / MSD exponent
(standard Green–Kubo); HMC-as-ballistic-transport and the U-turn/optimal-length
picture (Betancourt 2016 XHMC/NUTS, `references/papers/betancourt_2016_xhmc`,
eq:virial, eq:nuts); internal-coordinate ballistic torsional transport
(Forrest & Suter 1994 GC-HMC, `references/papers/forrest_1994_gchmc`; Chen et al.
2005 TAMD, `references/papers/chen_2005_tamd`, 10–20 fs steps from removing stiff
modes).

## 3. Condition for ballistic survival across a clash (question 1)

**CLAIM T1 (elastic vs thermalizing contact).** Coherent `u_ξ` survives a clash
(elastic exchange) rather than thermalizing when (i) the contact is impulsive —
contact duration short vs the internal relaxation time of the DOF that receive the
energy — and (ii) the effective number of strongly-coupled receiving DOF is small
(low-dimensional bath ⇒ energy recurs rather than disperses). In GB/OBC two-robot
contact the bath is exactly this low-dimensional conservative set (`00-...` §4),
so genuine thermalization is WEAK and `τ_p` should be long relative to a single
contact. A *strong* measured friction (short `τ_p`, `α→1`) is therefore prima
facie numerical (`30-...`), not physical.

**CLAIM T2 (chaotic-mixing caveat).** The one physical route to fast
thermalization even with few DOF is chaotic (Lyapunov) mixing at a hard clash;
its rate grows with clash stiffness and trajectory length. This is the physical
ceiling on `L`: beyond the Lyapunov time the coherent momentum is scrambled and
extra steps are wasted. Measured as the `L` at which `α` crosses from 2 toward 1.

**CLAIM T3 (win condition — MAKE OR BREAK).** The mixed move beats Cartesian
MD/HMC per force-eval when a single coherent flight spans the slow-dihedral
barrier — i.e. when the *per-refresh* ballistic mean free path exceeds the barrier
width:

```
ℓ_traj = v_th · min(L·dt, τ_p, t_Uturn)  >  Δξ_barrier.                 (WIN)
```

*Why `ℓ_traj`, not the intrinsic `ℓ = v_th·τ_p`.* Under full refresh (the default;
§4) coherent transport is capped at the trajectory length `L·dt`, not the intrinsic
`τ_p`: the refresh randomizes direction between trajectories, so cross-trajectory
motion is a random walk with step `ℓ_traj = v_th·min(L·dt, τ_p, t_Uturn)`. In the
collapsed-`dt` regime `L·dt ≪ τ_p`, so `ℓ_traj ≪ ℓ` and the intrinsic `ℓ` is
over-optimistic; the go/no-go SHALL use `ℓ_traj` at the working `L·dt`.

*Threshold, not a closed-form speedup.* `ℓ_traj > Δξ_barrier` is a *condition* (one
flight clears the barrier), not a speedup formula. A quantified speedup requires
comparing the mixed move against Cartesian's OWN diffusion constant `D_cart`
(measured separately), because in the winning regime the mixed move is ballistic —
its own `D = v_th²τ_p` does not describe it. NOTE: the earlier "faster by
`ℓ/Δξ_barrier`" ratio is WITHDRAWN — it mis-inverted `t_d/t_b = Δξ_barrier/ℓ` and
reused the mixed move's `D` in a regime where the move is not diffusive. □ Viability
is decided by the *measured* inequality above plus a `D_cart` cross-check, not a
formula.

## 4. Trajectory length and momentum policy (question 2)

**Full resampling (current).** `(u, v_E)` redrawn every trajectory
(`reinitialize`). Per-trajectory transport is ballistic up to the U-turn, but the
refresh randomizes direction between trajectories ⇒ the chain of endpoints is a
random walk with step `ℓ_traj = v_th·min(L·dt, τ_p, t_Uturn)`. At the collapsed
`dt`, `L·dt` is short, so `ℓ_traj` is short and the chain looks diffusive — this
is exactly the reported "shortened ballistic / lengthened diffusive" symptom, and
it is partly a POLICY artifact, not physics.

**Optimal length / U-turn (mixed system).** The trajectory SHOULD stop at the
U-turn of the SLOW coordinate `ξ`, not of the full system (whose NUTS criterion is
dominated by the fast contact mode and would truncate far too early). Slow-`ξ`
U-turn: stop when `(ξ(t)−ξ(0))·u_ξ(t) < 0` (Hoffman & Gelman NUTS, Betancourt
2016 eq:nuts, applied to the `ξ`-projection). Equivalently choose
`L·dt ≈ min(t_Uturn(ξ), τ_p)`: no benefit past `τ_p` (momentum already
decorrelated, T2) or past the `ξ` half-period (U-turn re-correlates).

**Partial momentum refreshment (Horowitz/GHMC) — the transport-preserving
option.** Replace the full redraw with a Horowitz rotation that mixes only a small
noise fraction (Horowitz, *Phys. Lett. B* 268(2):247, 1991; Kennedy & Pendleton
2001):

```
(u, v_E) ← cos φ_H · (u, v_E) + sin φ_H · ξ_MB,   ξ_MB ~ N(0, RT·M_mixed⁻¹),
```

applied CONSISTENTLY to both blocks with the correct block metrics (M3). Small
`φ_H` keeps momentum coherent ACROSS trajectory boundaries, extending `ℓ` over
many trajectories and suppressing the random walk — the direct antidote to the
policy-induced friction above.

**CLAIM P1 (the GHMC tension).** Partial refreshment REQUIRES momentum reversal on
rejection (P2 in `10-...`), and each rejection *reverses* the coherent momentum,
undoing ballistic progress. So GHMC helps ONLY when inner acceptance is high
(short trajectory per refresh, or small `dt`) — the momentum-flip cost grows with
rejection rate (partial-refreshment chains "explore more slowly due to momentum
reversals on rejection"). Reduced-flip variants (Sohl-Dickstein / "HMC with
reduced momentum flips", arXiv:1205.1939) mitigate this.

**RECOMMENDATION (policy is data-driven, not assumed).**
- Whether a single MH over ~1 ps (`L=500`, `dt=2fs`) is right depends on the
  measured `τ_p` and slow-`ξ` U-turn time. If `τ_p ≳ 1 ps` and the `ξ` U-turn is
  `≳ 1 ps`, the single long trajectory is fine. If `τ_p ≪ 1 ps` (fast
  thermalization) the extra steps are wasted and shorter trajectories with partial
  refreshment dominate. So `L` SHALL be set from §2's measurement, not fixed a
  priori.
- Partial refreshment is RECOMMENDED specifically to combat the collapsed-`dt`
  policy friction, CONDITIONAL on high inner acceptance (P1) and on the mixed-block
  flip (P2). Under high clash-rejection it degrades below full-refresh HMC.

**Why NOT multiple-time-stepping (RESPA) at contact — context for question 2.**
Splitting the stiff contact force onto a small inner `dt` (RESPA) is rejected
because impulse-RESPA has a linear resonance instability near half the fastest
period and nonlinear instabilities at ~1/3–1/4 of it (Barth & Schlick 1998; Ma,
Izaguirre & Skeel, *SIAM J. Sci. Comput.* 24:1951, 2003, "Verlet-I/r-RESPA/Impulse
is Limited by Nonlinear Instabilities"). The stiff contact mode is precisely the
fast mode, so RESPA cannot lift `dt` past that resonance — it buys nothing over
the collapsed Verlet `dt`, while breaking the single-map reversibility M4 needs.
The mixed ballistic-transport route (this file) and local soft-core NCMC (`40-...`)
are the two admissible responses.

## 5. Touch list
- New diagnostic logging of `u_ξ` and `ξ` per NVE step (diagnostic, not physics)
  to compute `C(t)`, `τ_p`, `α`, `ℓ`. Selecting `ξ` (which `u` component / which
  dihedral) is a per-system input.
- `World::reinitialize` (`:1608`): OPTIONAL partial-refreshment path (Horowitz
  rotation of `(u, v_E)` with block metrics) behind a flag; retains full refresh
  as default. Convention at risk: applying `φ_H` to `u` (metric `M_φ`) and `v_E`
  (metric `M_E`) with the correct per-block noise; reject-flip of both blocks
  (`10-...` P2).
- Slow-`ξ` U-turn stop as an OPTIONAL trajectory-termination criterion in the
  `L`-loop (`World.cpp:1457`); default remains fixed `L`.
