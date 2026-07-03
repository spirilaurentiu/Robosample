# NCMC explicit-solvent acceptance — Acceptance construction

## 1. Scope
Two exact NCMC acceptance constructions; the exact acceptance ratio and unbiasedness sketch for each; the matched-pair theorem; the recommended target for explicit solvent. Terms from `00-diagnosis-and-scaling.md` §2.

## 2. The two constructions and the matched-pair theorem

**CLAIM C1 (matched pair).** Exactness requires acceptance rule and inner-propagation type to match:
- **Construction I — endpoint-ΔH (current).** Inner propagation = unadjusted deterministic Verlet. The whole switch is one deterministic, volume-preserving, F-reversible map `T` on `(q,u,x_s,v_s)`; accept `min(1,e^{−βΔH})`. Exact by the HMC involution argument (A1). Shadow work `Q` is *in* `ΔH`.
- **Construction II — Metropolized-dynamics NCMC (target).** Inner propagation `Φ_λ` = a kernel that is **exactly `π_λ`-invariant** (per-step GHMC accept/reject against `H_λ`, with momentum reversal on reject). Accept the outer move `min(1,e^{−βW})` on the protocol work alone. Exact by the NCMC path-space MH argument (Nilmeier et al. 2011, Eqs. 12–13). Shadow work is absorbed into inner rejections and is **not** in the acceptance.

**Crossing them biases the sampler:** accept-on-`W` with plain Verlet omits the shadow work that genuinely perturbs the distribution (biased); accept-on-`ΔH` with stochastic non-deterministic inner steps breaks the single-map involution (biased). SHALL NOT mix.

## 3. Construction II — exact ratio and unbiasedness sketch

Extended state `x=(q,u,x_s,v_s)`; target `π(x) ∝ exp(−β·H_{λ=1}(x))` with `H` spelled EXACTLY as `currentTotalEnergy()` (`World.cpp:1788`) assembles it:

```
H = V + K(u;q) + K_s(v_s) + U_F(q) − ½·RT·ln sin²γ₂(q) − nmaCorr
U_F(q) = +½·RT·ln( |M_{Nf}(q)| / |M_{3N}| )   (calcFixman, World.cpp:864)
```

NOTE (F2): do NOT read this as `|M(q)|^{1/2} · exp(−β[V+K+K_s])` with `|M|^{1/2}` a separate explicit multiplicative factor. The `|M(q)|^{1/2}` weighting is the *realized* configuration-space marginal produced BY `U_F` together with the Maxwell–Boltzmann velocity resample; it is already inside `exp(−βH)` via `U_F`, and the density also carries the pitch term `−½RT ln sin²γ₂` and `−nmaCorr` that the shorthand drops. A coder who trusts the shorthand over this expansion implements exactly the F1 bug (inner accept omitting `U_F`/pitch). The inner GHMC accept SHALL evaluate this full `H_λ`.

Protocol = alternating kernels

```
K₀ , Φ_{λ0} , K₁ , Φ_{λ1} , … , K_{n-1} , Φ_{λ_{n-1}}
```

with `K_i` the deterministic fixed-`x` perturbation `λ_{i-1}→λ_i` (accumulates `w_i = V_{λi}(x) − V_{λi-1}(x)`, moves nothing) and `Φ_{λi}` a `π_{λi}`-invariant propagation. Schedule palindromic with `λ₀=λ_n=1` (`NcmcProtocol.hpp`, unchanged).

**Acceptance:** `a = min(1, exp(−β·W))`, `W = Σ_i w_i`.

**Unbiasedness sketch.** In NCMC path space the composite move is Metropolis–Hastings; the acceptance is the ratio of reverse-path to forward-path weights. Each `Φ_{λi}` satisfies detailed balance w.r.t. `π_{λi}`, so its forward/reverse contribution is exactly `π_{λi}(x')/π_{λi}(x)`; chaining these across the protocol telescopes, and every term except the perturbation energy bookkeeping cancels, leaving `exp(−β Σ_i w_i)`. Crucially, `Φ`'s accept/reject makes it `π_λ`-invariant *regardless of the integrator's discretization error*: a rejected inner step contributes zero net displacement and zero net work, so the integrator's shadow work never enters `W`. The palindromic, cyclic (`λ:1→1`) protocol has `ΔF=0`, so by Jarzynski `⟨e^{−βW}⟩ = e^{−βΔF} = 1` — the mean acceptance weight is 1, and acceptance is set by `Var(W)`, i.e. the *dissipation* `⟨W_diss⟩ = ⟨W⟩ − ΔF = ⟨W⟩ ≥ 0`. □

**INVARIANT (matched-pair, must-check).** A `Φ` that is only unadjusted Verlet (not `π_λ`-invariant) with accept-on-`W` SHALL fail a known-ΔF check; the same protocol with a genuine `π_λ`-invariant `Φ` SHALL pass. This is the discriminating oracle (`40-reproducer-and-oracles.md`).

## 4. Scaling of the Construction-II residual

**CLAIM C2.** Construction II removes the `N_bath`-shadow scaling from acceptance. The residual acceptance cost is `Var(W)` where `W` is the reorganization work of switching Region A's *intermolecular* coupling off and back on. For a fixed-size Region A in a bath larger than the correlation length, the solvation-response work is set by A's solvation-shell reorganization and is **not** extensive in box size. The bath DOF still cost *inner-move rejections* (mixing time), but no longer crush the outer acceptance. This is the viability argument for FFAR1 + explicit solvent.

**Residual that does NOT vanish.** `⟨W_diss⟩ > 0` and grows with (i) coupling strength, (ii) the shell that must reorganize, (iii) how far the λ=0 basin-hop displaces A relative to the shell's relaxation during `Φ`. Mitigated by more switching stages and by giving `Φ` enough inner steps to relax the shell (now free of shadow-work penalty) and by minimum-dissipation protocols (Sivak, Chodera, Crooks, *PRL* 108:190602, 2012; *J. Chem. Theory Comput.* 10:2803, 2014). This is the physically irreducible cost NCMC exists to pay down reversibly; it is categorically different from shadow work (pure numerical waste).

## 5. Reconciliation with `relax_solvent`
`relax_solvent` supplies the shell relaxation `Φ` needs, but under Construction I every relaxation step is charged as shadow work (`00-...` §3.1). Under Construction II the *same* Cartesian-Verlet relaxation becomes pure benefit (lowers `W_diss`, costs nothing in acceptance) **once each inner step is Metropolized**. So `relax_solvent` is retained and its inner loop is wrapped in the GHMC accept/reject of `20-inner-integrator.md`.

NOTE (F3, scope limit on "pure benefit"): "costs nothing in acceptance" holds ONLY where the inner GHMC *proposal* map stays volume-preserving and F-reversible. A non-converged velocity corrector (`RobotIntegrator.hpp:396`, "take the step anyway") is not guaranteed reversible, so as an inner proposal it breaks `π_λ`-invariance and silently biases the target — this is a *correctness* failure under Construction II, not the self-correcting shadow-work of Construction I. See `20-inner-integrator.md` §2–§3: under Construction II a non-converged inner step is a reject / fail-loud / reduce-dt condition.

## 6. Touch list (Construction II)
- `World::ncmcMove` (`:1892`): acceptance switches from `metropolis(Hstart,Hend)` to `min(1,e^{−βW})`; `W` is the already-computed `work`. Momentum treatment: outer full resample stays (Gibbs velocity update); inner GHMC owns the reject-flip.
- New inner kernel wrapping `RobotEngine::stepTo` per substep: evaluate `H_λ` before/after the joint (solute+solvent) step, accept/reject, on reject restore `(q,u,x_s,v_s)` and negate `(u,v_s)`. Requires per-substep `H_λ` including `U_F(q)` and `K_s(v_s)`.
- `include/NCMCProtocol.hpp`: unchanged (palindrome still needed so the perturbation ledger `W` is the self-reverse work; the reversibility that mattered for Construction I's single map is now supplied per-`Φ` by GHMC, but the palindrome keeps `ΔF=0` and the `λ=1` gate).
- Conventions at risk: `H_λ` must use the *same* Fixman/pitch/solvent-KE terms in the inner accept as the outer `H` (`currentTotalEnergy`), or the inner kernel is not `π_λ`-invariant. Frame F/M and angular-over-linear conventions unchanged.

NOTE. Construction I remains correct and SHOULD be retained (a `use_metropolized_inner` flag), because for a small rigid ligand with a clean corrector it is cheaper (no inner rejections). Construction II is REQUIRED for explicit-solvent / large-bath viability.
