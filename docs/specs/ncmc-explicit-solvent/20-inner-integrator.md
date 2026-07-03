# NCMC explicit-solvent acceptance — Inner integrator

## 1. Scope
The propagator requirement for each construction, classified NORMATIVE vs RECOMMENDED, and the corrector-convergence criterion. Answers mission Q2. Terms from `00-diagnosis-and-scaling.md` §2.

## 2. What each construction requires of the propagator
- **Construction I (endpoint-ΔH):** the step map SHALL be deterministic, volume-preserving, and F-reversible (certified by `checkReversibility`, `RobotIntegrator.hpp:471`). Symplecticity is *not required for exactness* but *is required for acceptance to not scale with `n`* (drift, `00-...` §4.2). So for Construction I, symplectic quality is a RECOMMENDED efficiency property, not a correctness one.
- **Construction II (accept-on-`W`):** each `Φ_λ` SHALL be exactly `π_λ`-invariant. The cheapest exact route is GHMC: one (or few) Verlet steps proposed under `H_λ`, Metropolis accept/reject against `H_λ`, momentum reversal on reject, with periodic (partial or full) momentum refresh. NORMATIVE.
  - **Inherits Construction I's reversibility precondition (F3, NORMATIVE).** GHMC `π_λ`-invariance requires the inner *proposal* map to be volume-preserving and F-reversible (what `checkReversibility` certifies, `:471`). Therefore Construction II inherits — does NOT escape — the Construction-I requirement on the propagator. The "take the step anyway" non-converged corrector path (`:396–407`) yields a state-dependent-iteration map that need not be F-reversible; used as the GHMC proposal it makes the inner kernel non-`π_λ`-invariant and biases the target. Under Construction II a non-converged inner step SHALL be handled as reject / fail-loud / reduce-dt, NEVER silently taken. Guard: bound the `checkReversibility` residual of the inner proposal at the working dt.

## 3. Corrector convergence (source S-b)
The internal velocity corrector is functional iteration of the implicit trapezoid `u1 = u0 + (h/2)(udot0+udot1)`, tol `1e-4` relative, max 10 sweeps, and on non-convergence the step is *taken anyway* with a warning (`RobotIntegrator.hpp:316–407`).

**CLAIM I1.** A non-converged corrector is a genuine, sign-definite shadow-work source (S-b), distinct from and additive to the fundamental symplectic floor (S-a). A step whose velocity is not the converged implicit-trapezoid solution is not the intended (near-symplectic, reversible) map and injects heat that accumulates linearly in `n`.

**CLAIM I1b (severity escalates under Construction II — F3).** Under Construction I a non-converged corrector is *only* an efficiency loss (its heat lands in `ΔH` and is rejected away, self-correcting). Under Construction II the same non-convergence is a *correctness* defect: as the GHMC proposal it can break F-reversibility, so the inner kernel no longer preserves `π_λ` and accept-on-`W` samples a biased density with no visible symptom. The classification below is therefore Construction-dependent: RECOMMENDED-efficiency under I, NORMATIVE-correctness under II.

**Requirements:**
- **NORMATIVE (Construction I efficiency guard):** when the corrector does not converge at the configured `dt`, the move SHALL be diagnosable as corrector-limited (the reproducer counts warnings and reads `checkReversibility` residual). A non-converged step SHOULD be treated as a signal to lower `dt`, not tuned around by raising `ncmc_steps`.
- **RECOMMENDED:** expose corrector `tol` and `max_iters` as knobs; tightening `tol` below `1e-4` reduces S-b at fixed `dt`. Evaluate a RATTLE-consistent velocity solve (Andersen 1983, `references/papers/andersen_1983_rattle`, Eqs. 2.8/A4) so the velocity constraint is satisfied to solver tolerance rather than to a fixed 10-sweep cap.
- **NOTE (loop-closure RATTLE):** `constraints_.enforceVelocityConstraints`/`enforcePositionConstraints` (`RobotIntegrator.hpp:279,409`) are a separate constraint-corrector; for a cyclic Region A their residual is an additional shadow-work source. Acyclic 2ala has none, so this is not the 2ala cause but SHALL be assessed for cyclic receptors.

**CLAIM I2 (leverage).** If the reproducer attributes the reported `ΔH` mainly to S-b, a converged corrector is higher-leverage than reducing `dt` (removes a whole additive term at fixed `dt`). If S-a dominates (clean symplectic floor), only `dt`, a better integrator, or Construction II helps — and only Construction II removes the `N_bath` scaling. So the corrector fix is a *local* win for the 2ala case; it does not substitute for Construction II at receptor scale.

## 4. BAOAB / stochastic inner integrators
BAOAB Langevin (Leimkuhler & Matthews, *AMRX* 2013; *J. Chem. Phys.* 138:174102, 2013) has excellent *configurational* accuracy but is **not exactly `π_λ`-invariant** (`O(dt²)` configurational error). Therefore:
- It SHALL NOT be used as the `Φ_λ` in Construction II *without* a Metropolis correction (MALA/GHMC), because accept-on-`W` assumes exact `π_λ`-invariance.
- It MAY be used as the inner proposal *inside* the GHMC accept/reject (BAOAB-proposed, Metropolized), which restores exact invariance and can improve inner acceptance vs plain Verlet.

## 5. Touch list
- `include/RobotIntegrator.hpp` `verletStep` corrector block (tol/iters knobs; optional tighter/RATTLE solve). Convention at risk: the quaternion exp-map advance (lines 224–258) is deliberately non-Simbody and must stay reversible; do not "fix" it toward linear-Taylor.
- New GHMC wrapper (Construction II) around `stepTo`, owning per-substep `H_λ` accept/reject and reject-flip of `(u,v_s)`.
