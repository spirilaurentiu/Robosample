# Per-rigid-body spatial wrench as a GPCR allostery observable (holo vs apo)

Status: **COMPLETE** (except discovery/empirical open questions OQ2–OQ5, none of which blocks a grounded
normative claim). Every normative claim is grounded in authority (`references/`, codebase) or in a
CAS-closed derivation. The transform lemmas L1–L4 are now **proof-grade** (Wolfram Language kernel). The
seven prior-art papers (FDA, PRS, strain, GPCR-MD review) are verified to exist in peer-reviewed venues via
scite but are **not in `references/`**; they stay discovery-tier and are listed for ingestion.

Revision (2026-07-14, post spec-review re-derivation): resolved two Blocking findings at the root by
committing to the reaction-side reading. B1 — the `⟨F^B⟩→0` equivalence collapse is withdrawn; the CoM
choice rests on L1 (factorization-independence) alone, and the mean reference-dependence residual is a
reported quantity. B2 — CAS proves `â·N(Bo)` equals the ABA generalized force only for a tip body, so I5 is
corrected to the own-atom reaction-side reading and C-TM6 is re-grounded as a restrained-ensemble
hinge-torque sign test, separated from the headline `N_c^B`. S1–S4 folded into V6/V8/I6 and the C-TM6
pinning. L1–L4 CAS-closed proof-grade.

The key words MUST, MUST NOT, SHALL, SHALL NOT, SHOULD, SHOULD NOT, MAY are used per RFC 2119.

Evidence tier legend on each claim: **[A]** authority (`references/` verified equation or codebase
symbol), **[proof]** CAS-closed identity, **[D]** discovery-tier (existence verified, not yet authority).

---

## Motivation

### Problem restatement

**User's original phrasing (verbatim):**

> Robosample Paper 1: the per-rigid-body spatial wrench (net force + net torque about a body frame) as a
> GPCR allostery observable, compared holo vs apo. This observable is a byproduct of the articulated-body
> (ABA) recursion the sampler already runs. […] (1) TORQUE REFERENCE-FRAME CONVENTION […] fix ONE
> physically motivated, holo/apo-comparable convention […] (2) BENCHMARK AGAINST FORCE DISTRIBUTION
> ANALYSIS (FDA) AND PRS […] show the rigid-body TORQUE captures a rotational / lever-arm effect at
> helix/domain scale that pairwise atomic forces (FDA) do not directly expose […] (3) STATISTICAL
> CONVERGENCE […] Net torque nearly cancels for near-rigid domains, so it is tail-sensitive.

**Restatement in codebase vocabulary.** Robosample already reduces OpenMM per-atom applied forces to one
spatial wrench per rigid body via `reduceAtomForcesToBodies` (`src/bridge/ForceReducer.cpp`), producing
`RobotState::bodyForceG[b]`, a `SpatialVec` whose `[0]` (angular) is the net moment about the body origin
and `[1]` (linear) is the net force, in the Ground frame (INV-1; `include/bridge/ForceReducer.hpp:52-57`)
**[A]**. This is exactly Jain's effective cluster spatial force `f̂_c(k)`, eq:4.6, "torque part on top, net
force on bottom" **[A]** (`references/papers/jain_1993_on_recursive_dynamics/equations.md`, eq:4.6).
`World::captureReactionSnapshot` (`src/world/ReactionReporter.cpp:100`) already emits this per selected
body at DCD cadence, optionally plus the static (u=0) mobilizer reaction, both referenced at the body
origin Bo, Ground frame **[A]**.

The scientific problem: define, from this already-computed wrench, a **per-body observable comparable
between two ligand states of the same receptor (holo = agonist-bound, apo = ligand-free) and across a
database of GPCRs**, such that a change in the observable reports an allosteric, rotational,
helix/domain-scale mechanical signal — specifically the moment tending to drive a known GPCR activation
motion (TM6 outward movement; NPxxY / DRY / connector microswitch rearrangements).

The implementation problem that follows: the raw `bodyForceG` is **not** yet such an observable, for three
reasons that are the three reviewer objections:

1. **Reference-point dependence.** The moment `Σ_i (r_i − p) × f_i` depends on the reference point `p`. The
   current convention uses the body origin `p = X_GB[b].p()` — the inboard mobilizer frame — which is a
   **modeling choice** (where the joint was placed by the robot factorization), not a physical invariant.
   Two factorizations of the same molecule, or two receptors in a database, place the mobilizer frame at
   different atoms, so the reported moment is not directly comparable.
2. **No FDA/PRS differentiation stated.** Prior interatomic-force allostery methods exist and are cited to
   make binding-signature claims. The spec must state precisely what the rigid-body torque adds.
3. **Statistical fragility.** For a near-rigid domain in near-equilibrium the mean net force is small, so
   the mean torque is a small residual of a large-variance quantity; convergence must be certified, not
   assumed.

### Quantified constraints and assumptions

- **Target systems.** GPCRs at the scales Robosample targets (CLAUDE.md): FFAR1 (a GPCR) in implicit
  solvent and in an explicit membrane nanodisc; up to 100k rigid bodies / 10k robots / 1M atoms. A GPCR TM
  bundle is naturally factorized as ≈7 helical rigid bodies (torsional/rigid hinges), consistent with the
  GNEIMO rigid-cluster + torsional-hinge model **[A]** (`references/papers/vaidehi_2015_icmd_gneimo/paper.md`
  §2; `larsen_2014`).
- **Units (INV-3).** Engine MD system: length nm, mass dalton, time ps, energy kJ/mol
  (`include/Units.hpp:26-27`) **[A]**. Therefore net force is in **kJ·mol⁻¹·nm⁻¹**, and the moment `r×f` is
  in **kJ·mol⁻¹** (nm × kJ·mol⁻¹·nm⁻¹). These units SHALL be stated on every reported column.
- **Frame algebra.** `SpatialVec[0]=angular (moment)`, `[1]=linear (force)` (`include/robot_math.hpp:656-664`)
  **[A]**; the spatial-force reference-point shift is the codebase `PhiMatrix` operator
  `φ(offset)·f = {N + offset×F, F}` (`include/robot_math.hpp:862-875`) **[A]**, which is Jain's
  `φ(O_x,O_y)` **[A]** (`jain_1993/notation.md`, eq:A.2). Per-body CoM in Ground (nm) is available as
  `RobotState::comG()` **[A]** (`include/RobotState.hpp:311-317`).
- **Cost.** The observable is a post-processing contraction of quantities the ABA recursion already
  produces; it SHALL add no forward-dynamics cost beyond the already-optional snapshot.

---

## Behavior

### The observable (defined once)

Let body `b` be a rigid cluster with real (mass ≠ 0) atoms `i`, Ground positions `r_i`, and OpenMM applied
forces `f_i = −∂U/∂r_i`. Define, in Ground:

- Net force **F = Σ_i f_i** (`bodyForceG[b].linear`).
- Net moment about reference point `p`: **N(p) = Σ_i (r_i − p) × f_i**. The engine computes `N(Bo)` with
  `Bo = X_GB[b].p()` (`bodyForceG[b].angular`).

Let `R = X_GB[b].R` be the Ground←Body rotation and `c = comG[b]` the body CoM in Ground.

**Fixed convention.** The per-body allostery observable is the pair

> **`(F^B, N_c^B)`** = the net force and the net moment about the body center of mass, both rotated into
> the body's own frame:
> **`F^B = Rᵀ F`**, **`N_c^B = Rᵀ [ N(Bo) + (Bo − c) × F ]`**,
> ensemble-averaged **in the body frame**.

Reported alongside:

- the two screw invariants **`|F|`** and **`π = N(p)·F̂`** (reference-point-independent by L2, algebra
  alone, carrying no zero-mean assumption);
- the hinge-axis scalar **`T = â · N(Bo)`** = the own-atom applied moment about the hinge axis (equal to the
  ABA generalized force only for a tip body; see L4);
- the mean reference-dependence residual **`−(Rᵀ(c−Bo)) × ⟨F^B⟩`** (L3), so the reference-point dependence
  of the mean torque is explicit, not assumed away.

The allostery signal is the holo−apo difference **`ΔN_c^B`** per body (and `ΔF^B`), evaluated only after
the convergence criterion (V8) is met.

### CAS-closed lemmas (proof-grade)

All four reduce `lhs − rhs → 0` symbolically (Wolfram Language kernel; oracles in Validation). Vectors are
generic; sums over atoms are linear so per-atom / two-atom checks suffice. `R ∈ SO(3)` is a generic
Euler-angle rotation.

- **L1 (Varignon transform).** `N(p_B) = N(p_A) + (p_A − p_B) × F`, and `F` is reference-point-independent.
  CAS: `Cross[r−pB,f] − (Cross[r−pA,f] + Cross[pA−pB,f]) = {0,0,0}`. Independent witness: identical to the
  codebase `PhiMatrix` operator `φ(offset)·{N,F} = {N+offset×F, F}` (`robot_math.hpp:874`) **[A]**.
  **Grade: proof.**
- **L2 (screw invariants).** `N(p)·F` is reference-independent; hence `|F|`, the pitch `h=(N·F)/|F|²`, and
  `N_∥=hF` are reference-clean. CAS (two atoms): `Ntot[pB].Ftot − Ntot[pA].Ftot = 0`. **Grade: proof.**
  These are reference-clean by algebra alone — they carry **no** zero-mean assumption.
- **L3 (body-frame moment transform).** Averaging MUST be done in the body frame. Using
  `R·(a×b)=(R·a)×(R·b)` (CAS: `{0,0,0}` for the Euler-angle `R`), with `g=c−Bo`, `g^B=Rᵀg`, `F^B=RᵀF`:
  `N^B(Bo) = N^B(c) + g^B × F^B`, i.e. **`⟨N^B(c)⟩ − ⟨N^B(Bo)⟩ = −g^B × ⟨F^B⟩`** (CAS:
  `Rᵀ·Cross[g,F] − Cross[Rᵀg,RᵀF] = {0,0,0}`). **Grade: proof.**
- **L4 (own-atom hinge moment vs generalized force).** With the child spatial force `f(k−1)={N_c,F_c}`,
  offset `l`, own-cluster applied wrench `{N_own,F_own}`, and `f(k)=φ(l)·f(k−1)+{N_own,F_own}` (eq:4.7,
  inertial/gyro set aside): `T(k)=â·f(k)_ang` and `â·N(Bo)=â·N_own` differ by
  **`â·(N_c + l×F_c) = â·(φ(l)·f(k−1))_ang`** — the hinge-axis projection of the outboard-transmitted
  moment. They are **equal only for a tip body** (`f(k−1)=0`; CAS: `tip = 0`). **Grade: proof.**

### L3 — no equivalence-class collapse (resolves B1)

The premise `⟨F^B⟩ → 0` is **false** for the target bodies. It holds only for a free 6-DOF body, whose
ensemble mean force vanishes on its three translational DOFs. A TM helix is a torsional-hinge (or
rigid-hinge) body: its net atomic force `F` is the constraint reaction holding it against the bundle, with
**generically nonzero mean** `⟨F^B⟩ ≈ −⟨hinge reaction force⟩`. NOTE: Ground-frame momentum stationarity
`⟨d/dt(spatial momentum)⟩=0` is exact, but the body-frame average carries a rotational correction
(`⟨Rᵀ ṗ⟩ ≠ d/dt⟨Rᵀ p⟩`); the relation is used only to support the conservative claim that `⟨F^B⟩` is
generically nonzero, with its empirical magnitude deferred to OQ4, not as an exact equality.
"Near-rigid" bounds the *fluctuation* of the internal geometry, not the *mean* of the external reaction.

Consequently the CoM and mobilizer-frame mean torques **differ**, by the CAS-closed residual

> **`⟨N_c^B⟩ − ⟨N_Bo^B⟩ = −(Rᵀ(c−Bo)) × ⟨F^B⟩ ≠ 0`** (L3).

The equivalence-class-collapse argument is withdrawn. **The CoM reference point is chosen on L1 alone:** the
CoM is fixed by the body's atom set and masses, hence factorization-independent and identical holo/apo (same
atoms), which the mobilizer frame is not (its L1 residual `(Bo−c)×F` is a modeling artifact). This
justification needs no statement about `⟨F⟩`. The residual `−(Rᵀ(c−Bo))×⟨F^B⟩` SHALL be reported per body
alongside `⟨N_c^B⟩`.

### L4 / I5 — reaction-side reading (resolves B2)

`â·N(Bo)`, with `N(Bo)=bodyForceG[k].angular`, is the **applied moment about the hinge axis from body k's
own atoms** (eq:4.6, `Σ_{i=1..r(k)}`) — a reaction-side quantity. By L4 it equals the ABA generalized force
`T(k)=H(k)f(k)` (eq:4.7, which accumulates the whole outboard subtree plus inertial/gyroscopic terms)
**only for a tip body**. The earlier claim "`â·N(Bo)` is exactly the generalized force" is corrected to hold
**only in the tip case**, and is otherwise a partial (own-atom) contribution differing by the CAS-closed
transmitted term.

**Reconciliation with equipartition (the shared root).** For a free periodic torsion `θ(k)` at canonical
equilibrium, `⟨T(k)⟩ = ⟨−∂U/∂θ(k)⟩ = 0` (`∮ ∂U/∂θ · e^{−βU} dθ = 0` by periodicity). Therefore:

- The **free hinge-axis component** `â·N(Bo)` is zero-mean at free equilibrium — exactly for a tip body, up
  to the transmitted term otherwise. It CANNOT carry an allostery signal as a free-equilibrium mean; it is
  used only as a **stationarity sanity check** (I6/V8.3).
- The **constrained-direction reactions** — the net force `F` and the two off-hinge-axis moment components
  — are structurally **nonzero-mean** (`= −⟨hinge reaction⟩`), because the hinge resists those directions.
  These carry the reaction-side signal.

The spec commits to the reaction-side reading throughout: the allostery signal lives in the nonzero-mean
constrained-direction reactions; the reference-clean invariants `|F|`, `N·F̂` (L1/L2) are algebraically
reference-independent regardless of mean and carry a generically nonzero reaction-side mean.

### C-FDA — incremental contribution (unchanged in substance)

`N_c^B` is a signed, lever-arm-weighted vector contraction of the vector force field (eq:4.6). It is not a
function of the FDA scalar magnitudes `{|F_ij|}` or of the symmetric strain tensor: two force fields with
identical per-pair magnitudes but a sign-flipped tangential component give identical `σ_i` and identical
strain yet opposite `N_c^B`. The incremental contribution is the residual of `ΔN_c^B` after regressing out
the within-body `{Δσ_i}` (operationalized in V6). Rests on authority-tier eq:4.6 plus elementary algebra;
unaffected by B1/B2.

### C-TM6 — falsifiable, equipartition-consistent (resolves B2, S2)

TM6 outward movement is the class-A activation coordinate; NPxxY/DRY/connector are coupled microswitches
**[D]** (`10.1016/j.sbi.2019.03.016`). Because the free hinge-axis torque is zero-mean at free equilibrium,
C-TM6 SHALL be posed as one of two reaction-side, nonzero-mean tests, **explicitly distinct from the
headline observable `N_c^B`**:

> **C-TM6 (primary, restrained-ensemble hinge-torque sign test).** In a **reference-restrained ensemble**
> with `θ_TM6` held near the inactive (reference) geometry — so the torsion-periodicity identity no longer
> forces zero mean —
> the agonist-induced change in the hinge-axis applied moment,
> `Δ⟨â·N(Bo)⟩_TM6 = ⟨·⟩_holo − ⟨·⟩_apo`, SHALL have the sign that drives TM6 **outward**. Restraint breaks
> the periodic zero-mean; the sign is the predicted drive direction.

> **C-TM6 (complementary, free-equilibrium constrained-reaction test).** At free equilibrium in each state,
> the change in the constrained-direction reaction on TM6 — `Δ⟨F^B⟩` and the off-hinge-axis components of
> `Δ⟨N_c^B⟩` — is a nonzero-mean signal reporting the changed mechanical loading between states. It is
> descriptive, not a signed drive predictor.

**Falsifiability pinning (S2), all SHALL be fixed before the test:**
1. **Axis orientation.** `â` is the mobilizer rotation axis, oriented by the fixed rule: `+`rotation about
   `â` increases the distance from the TM6 cytoplasmic tip to the bundle axis in the reference (inactive)
   structure. Removes the `â` vs `−â` ambiguity.
2. **Torque→motion map.** With that orientation, `+â·N(Bo) ⇒` outward-driving moment. The predicted
   activation signal is `Δ⟨â·N(Bo)⟩ > 0`.
3. **Frame.** `â` is constant in the mobilizer/body frame; the scalar is rotation-invariant, but the
   ensemble mean SHALL be formed as `â_body · ⟨N^B(Bo)⟩` (body-frame average, L3) to avoid rotational
   washout.
4. **Scope.** The falsifiable anchor tests the **hinge scalar `T` at the mobilizer frame**, NOT the
   headline CoM/body-frame observable `N_c^B`. The two SHALL NOT be conflated: `N_c^B` is the primary
   reported convention; `T` (restrained) is the pre-registered activation predictor.

*Failure:* a converged, significant `Δ⟨â·N(Bo)⟩` of the wrong sign, or null when TM6 is known to move,
falsifies the observable as an activation reporter for that system (report as a quantified bound).

### Holo/apo and cross-database comparability — the precondition it rests on

The observable is comparable between two states iff they share the **same atom→body assignment and the same
robot factorization**. For holo vs apo of one receptor this is automatic. For a **database**, comparability
requires a shared factorization protocol: SHALL fix one canonical body per TM helix with the mobilizer at a
conserved position (e.g. the GPCRdb generic-numbered residue), so that `â` and the body atom set map across
receptors (PRECONDITION-DB, OQ3). Without it, only the reference-invariant pair (`|F|`, `π`) from L2 is
comparable across receptors, and `N_c^B` is comparable only within a receptor.

---

## Invariants

- **I1 (force invariance).** `F` is reference-independent, equals `bodyForceG[b].linear`. **[proof, L1/A]**
- **I2 (transform law).** `N(p_B)=N(p_A)+(p_A−p_B)×F`, via `PhiMatrix`; Bo→CoM uses
  `offset = X_GB[b].p() − comG[b]`. **[proof, L1/A]**
- **I3 (screw invariance).** `|F|` and `N·F̂` identical about Bo or `c` to machine precision. **[proof, L2]**
- **I4 (body-frame averaging).** Means accumulated on `Rᵀ`-rotated components; `R` orthogonal preserves
  `|F|`, `|N|`, `N·F`. Ground-frame averaging is prohibited (biases the mean toward zero, L3).
  **[proof, L3/A]**
- **I5 (hinge-axis identity — CORRECTED).** `â·N(Bo)` is the applied hinge-axis moment from body k's **own
  atoms**; it equals the ABA generalized force `T(k)=H(k)f(k)` **only for a tip body**, otherwise differing
  by the outboard-transmitted moment `â·(φ(l)f(k−1))_ang` (L4). For a tip body it is zero-mean at free
  equilibrium (torsion-periodicity identity `∮ ∂U/∂θ · e^{−βU} dθ = 0`, not the equipartition theorem).
  **[proof, L4/A]**
- **I6 (regime / stationarity check — CORRECTED, resolves S3).** I6 does **not** assert
  reference-independence. It is a stationarity sanity check applied **only to the free-DOF projection**:
  `|â·⟨F^B⟩_hinge-projected|` (and, for a tip body, `|⟨â·N(Bo)⟩|`) SHALL be within its SE of the
  torsion-periodicity value 0. Passing it certifies the sampled DOF is stationary; it does **not** buy
  reference-independence of the full-3-vector mean (which L3 shows is generically reference-dependent,
  residual `−g^B×⟨F^B⟩`). The full-3-vector `⟨F^B⟩` is expected nonzero for constrained TM bodies and is
  reported, not gated out.
- **I7 (read-only w.r.t. sampling).** Computing the observable SHALL NOT perturb the accepted `q`/`u` or the
  sampling distribution — it consumes `bodyForceG` (velocity-independent) and, if the optional static
  reaction is included, restores `u` and all `u`-derived caches, exactly as `captureReactionSnapshot`
  already does (`ReactionReporter.cpp:139-164`). **[A]**
- **I8 (mass/virtual-site consistency, INV-1/INV-2).** The reduction skips massless slots; `F` equals the
  system net applied force to machine precision. **[A]**

---

## Interface

- **Touch — post-processing only (no dynamics change).** The observable is a contraction of already-emitted
  quantities `bodyForceG[b]` (SpatialVec), `X_GB[b]` (Transform), `comG[b]` (Vec3), all present in
  `RobotState`.
  - `world/ReactionReporter` (`src/world/ReactionReporter.cpp`, `captureReactionSnapshot`): after buffering
    the Bo-referenced wrench, additionally compute and buffer `(F^B, N_c^B, |F|, π, â·N(Bo))` and the
    residual `−(Rᵀ(c−Bo))×F` per interesting body, using `PhiMatrix(X_GB[b].p() − comG[b]) * bodyForceG[b]`
    then rotating by `X_GB[b].R.transpose()`. Existing `robot_math` operations; no new math kernel.
  - `OutputWriter` (CSV): new columns per body, with explicit units (force kJ·mol⁻¹·nm⁻¹; moment kJ·mol⁻¹)
    and a header naming the frame (body frame) and reference point (CoM). The existing per-frame reaction
    CSV is the carrier.
- **Touch — none in the ABA solver, mass operator, integrator, or acceptance.** Diagnostic only. It SHALL
  NOT enter the Hamiltonian, the Metropolis–Hastings acceptance, or the Fixman term.
- **`guidance / acceptance` split.** The new term belongs entirely on the **reporting** side, on **neither**
  the guidance nor the acceptance side of any move. No implementer routes a force/torque diagnostic into
  `-dU/dq` or into `ΔE`.
- **Analysis (offline, Python).** Body-frame accumulation, block-averaging / autocorrelation, holo−apo
  differencing, the FDA-residual regression, and the restrained-ensemble TM6 test are offline analysis over
  the CSV/DCD; they touch no engine code.
- **Preconditions surfaced at the boundary.** `setReactionReporter`/`enableReactionReporter` already reject
  a Cartesian world (`ReactionReporter.cpp:22-32`). The database-comparability precondition (shared
  factorization) is an input-construction contract, enforced by the world builder.

---

## Validation strategy

CAS oracles are closed (proof-grade). V6/V8 revised per S4/S1; I6/V8.3 per S3.

- **V1 — INVARIANT, proof.** L1: `Cross[r−pB,f] − (Cross[r−pA,f]+Cross[pA−pB,f]) = {0,0,0}` (closed).
  Numeric witness: `PhiMatrix(Bo−c)*bodyForceG[b]` reproduces direct `Σ_i(r_i−c)×f_i`, `< 1e-10`.
  **model_independent: true** (numeric leg).
- **V2 — INVARIANT, proof.** L2: `Ntot[pB].Ftot − Ntot[pA].Ftot = 0` (closed). *Fails on:* treating any
  component other than `N_∥` as reference-clean.
- **V3 — INVARIANT, proof + falsification.** L3: `Rᵀ·Cross[g,F] − Cross[Rᵀg,RᵀF] = {0,0,0}` (closed).
  Falsification leg (**model_independent: true**): inject a constant body-frame couple on a body sampled
  over a wide rotational spread; the Ground-frame mean decays toward zero while the body-frame mean is
  unbiased — catches Ground-frame averaging.
- **V4 — INVARIANT, proof.** L4: `â·f(k)_ang − â·N_own = â·(N_c+l×F_c)` and `tip = 0` (closed). *Fails on:*
  asserting `â·N(Bo)` is the generalized force for a non-tip body; asserting a nonzero free-equilibrium mean
  for a tip-body hinge component. Cross-checked against `TestReactionForces` free-mobilizer law. **[A]**
- **V5 — PRECONDITION.** Net-force closure `|Σ_i f_i − F_system,b| < 1e-9`, virtual sites skipped
  (INV-1/INV-2), via `TestForceReducer` host-vs-CUDA parity. **[A]**
- **V6 — LEMMA, falsification (revised per S4). model_independent: true.** Per TM body compute FDA
  punctual-stress change `Δσ_i` and the reaction-side `ΔN_c^B` (constrained-direction components). Regress
  `ΔN_c^B` on `{Δσ_i}`; **null model:** permutation null on the `{Δσ_i}` labels (≥10³ permutations) **and**
  a bootstrap CI on the residual variance across trajectory blocks. Accept "incremental" only if the
  observed residual exceeds the 95th percentile of the permutation null **and** the residual CI excludes 0
  **and** the residual's projection on the TM6 outward direction has a fixed sign at a stated threshold
  (e.g. `|cosθ| > 0.5`). *Fails on:* residual explained by estimator noise. Deterministic witness: V7
  counterexample (identical `{|F_ij|}`, opposite `N_c^B`), CAS-checkable.
- **V7 — LEMMA, falsification (deterministic).** Two synthetic force fields on one rigid body with identical
  `{|F_ij|}` and identical symmetric strain tensor but opposite tangential sign; confirm equal `σ_i`, equal
  strain, opposite `N_c^B`. **model_independent: true**.
- **V8 — INVARIANT, falsification (revised per S1). model_independent: true.** Tail-sensitive near-cancelling
  mean; the dominant failure is a **false blocking plateau** when trajectory length is not `≫ τ_int` (the
  methodology that sank the Fixman slow tier). Acceptance, all SHALL hold:
  1. **Blocking plateau** — block SE stable within 10% over a decade of block sizes;
  2. **Length-vs-`τ_int` stability guard** — `τ_int` by both Flyvbjerg–Petersen blocking and the
     initial-positive-sequence estimator, agreeing within 30%, **and** stable (within 20%) under truncating
     the trajectory to its first and second halves;
  3. **Tail-dominance / jackknife** — jackknife the largest-magnitude single-frame contributors; the
     estimate SHALL NOT move by more than its SE when the top 1% of contributors are dropped;
  4. **Effective sample size** — `N_eff = N/(2τ_int) ≥ 200` per state per component (at `N_eff=50` the
     relative SE on σ is ~14%, marginal for a 2σ claim; `N_eff≥200` gives ~7%);
  5. **Stationarity** — I6 free-DOF check passes.
  A holo−apo difference is a signal only if `|Δ| > 2√(SE_holo²+SE_apo²)`. A null is reported as the
  quantified upper bound `|Δ| < 2·SE_combined`, not as "no effect."
- **V9 — PRECONDITION.** Block-order invariance of the state-function observable at fixed accepted `q`,
  recomputed on the same frames in two block orders. **[A]**

---

## Consequences and trade-offs

- **Forecloses** reporting the raw `bodyForceG` moment as a cross-system observable (L1:
  factorization-dependent). The DOF-native scalar is retained only as the own-atom hinge moment `â·N(Bo)`,
  which equals the generalized force merely for tip bodies (L4).
- **CoM vs mobilizer frame is a mean-level distinction, not a fluctuation-level one** (corrected from the
  withdrawn L3(ii) claim): the means differ by `−g^B×⟨F^B⟩`, generically nonzero for constrained TM bodies.
  CoM is preferred purely for factorization-independence (L1). The residual is reported, not neglected.
- **The signal is reaction-side.** The zero-mean free-DOF component is a diagnostic, not the signal; the
  constrained-direction reactions (nonzero mean) carry the allostery information. Any future use as a
  bias/steering force is out of scope (would need a separate detailed-balance spec).
- **Database comparability** still requires the shared-factorization precondition; without it only the
  L1/L2 invariants (`|F|`, `N·F̂`) compare across receptors.

---

## Open questions

- **OQ1 — CLOSED.** L1–L4 are CAS-verified proof-grade this session (Wolfram Language kernel). Wolfram|Alpha
  cloud context was offline; the symbolic kernel executed and closed every identity.
- **OQ2 — Ingestion of the seven prior-art DOIs** (discovery-tier, not authority). None blocks a grounded
  claim (C-FDA rests on eq:4.6 + algebra). scite returned existence/venue only; a citation-status pass
  SHOULD confirm none is contested before publication.
- **OQ3 — Canonical database factorization.** Exact per-helix body boundaries and conserved hinge residues
  (GPCRdb generic numbering) needed for cross-database comparability (PRECONDITION-DB).
- **OQ4 — Reference-dependence of the mean is now a first-class reported quantity, not an assumption.** The
  residual `−(Rᵀ(c−Bo))×⟨F^B⟩` is reported per body; whether it is small for any given system
  (implicit-solvent vs explicit nanodisc) is empirical and SHALL be measured, not assumed.
- **OQ5 (new) — Restrained-ensemble protocol for C-TM6.** The primary activation test requires a
  reference-restrained ensemble (`θ_TM6` near inactive geometry). The restraint form, strength, and the
  correction for restraint bias on the reported mean torque SHALL be specified before the test is run; an
  over-stiff restraint suppresses the very signal it measures.

---

## References

Authority (`references/`, verified equations / codebase symbols): `jain_1993` eq:4.6/4.7 and Appendix A
φ-transform; `spiridon_2020_robosample`; `vaidehi_2015_icmd_gneimo`, `larsen_2014`; `sherman_2011_simbody`;
codebase INV-1 (`bridge/ForceReducer.{hpp,cpp}`), INV-3 units (`Units.hpp`), `world/ReactionReporter.cpp`,
`robot_math.hpp` (`SpatialVec`, `PhiMatrix`), `RobotState::comG/X_GB`, `TestReactionForces`,
`TestForceReducer`.

Discovery (existence/venue verified via scite, not yet authority; cited by title + DOI, author fields not
returned this session):

- *Implementation of force distribution analysis for molecular dynamics simulations.* BMC Bioinformatics,
  12, 101. https://doi.org/10.1186/1471-2105-12-101
- *Mechanical Network in Titin Immunoglobulin from Force Distribution Analysis.* PLoS Computational Biology.
  https://doi.org/10.1371/journal.pcbi.1000306
- *Dynamic Allostery in the Methionine Repressor Revealed by Force Distribution Analysis.* PLoS
  Computational Biology. https://doi.org/10.1371/journal.pcbi.1000574
- *Dynamic Allostery of the Catabolite Activator Protein Revealed by Interatomic Forces.* PLoS Computational
  Biology. https://doi.org/10.1371/journal.pcbi.1004358
- *Change in Allosteric Network Affects Binding Affinities of PDZ Domains: Analysis through Perturbation
  Response Scanning.* PLoS Computational Biology. https://doi.org/10.1371/journal.pcbi.1002154
- *Strain analysis of protein structures and low dimensionality of mechanical allosteric couplings.* PNAS.
  https://doi.org/10.1073/pnas.1609462113
- *Allostery in G protein-coupled receptors investigated by molecular dynamics simulations.* Current Opinion
  in Structural Biology. https://doi.org/10.1016/j.sbi.2019.03.016
