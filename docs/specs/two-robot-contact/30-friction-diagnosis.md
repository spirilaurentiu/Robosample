# Two-robot contact world — Is the "friction" physical or numerical? (question 4)

Terms from `00-...` §2. Separates genuine coupling-induced thermalization from the
documented Free-joint KE pump and the non-converged corrector.

## 1. The three candidate sources of effective friction
Effective friction `γ_eff` (the decay rate of `C(t)`, `20-...` §2) decomposes:
- **(γ-phys)** physical coupling of `u_ξ` to intra-R and R–E DOF. In GB/OBC this
  is elastic, recurrent, and WEAK (T1); total energy conserved.
- **(γ-pump)** the Free/Ball root-joint KE pump (`RobotIntegrator.hpp:216-227`),
  ~+1300 kJ/mol/traj at 1 fs from the port's reimplemented Free-joint
  `(N, Ndot, qddot)`; a SECULAR, sign-definite energy SOURCE (negative friction /
  heating) localized at the root body and propagated through the ABA recursion.
  The exp-map quaternion advance (`:224-258`) is the existing workaround.
- **(γ-corr)** the non-converged velocity corrector taken anyway
  (`:400-407`): a sign-definite heat source (positive friction / dissipation),
  additive in `L`.

`γ-pump` and `γ-corr` are numerical; `γ-phys` is real. The testbed (GB/OBC, two
robots) is chosen so `γ-phys` is small, sharpening the separation.

## 2. The diagnostic (NVE, refresh OFF)
Run ONE long NVE trajectory with momentum refresh DISABLED (pure Hamiltonian
flow), and measure:

- **Total-energy drift** `⟨H(t) − H(0)⟩` over `L` steps. A correct symplectic-ish
  map: BOUNDED, oscillatory, `O(dt²)` fluctuation, NO secular term. Secular linear
  growth ⇒ energy pump/sink (γ-pump or γ-corr).
- **Per-DOF kinetic-energy flow**: track KE in (a) `ξ`, (b) other R DOF, (c) E.
  Elastic (γ-phys): KE sloshes back and forth, recurrent, total conserved.
  Thermalizing (γ-phys, chaotic T2): KE relaxes toward `kT/2` per DOF, total
  conserved. Pump (γ-pump): total KE grows secularly, sourced at the ROOT body.
  Corrector (γ-corr): drift concentrated on steps where the corrector warns.
- **Localization**: attribute the energy source by body. Root-body-localized
  secular growth matching the Simbody-vs-port Free-joint discrepancy ⇒ γ-pump.
  Contact-localized `O(dt²)` drift vanishing as `dt→0` at fixed physical time ⇒
  γ-phys discretization friction (benign; Metropolis-corrected under full refresh).

**Tie to the Simbody differential oracle.** The KE pump's root cause is a known
Free-joint kinematics defect (`RobotIntegrator.hpp:227` ROOT-CAUSE TODO). The
diagnostic SHALL cross-check the per-body energy source against the Simbody
differential oracle (`.claude` memory: robotics-oracle-campaign; the
disasm-vs-refactor ABA differential test) on a single free body: if the two-robot
drift localizes to the root and matches the oracle's Free-joint discrepancy, it is
γ-pump and the exp-map workaround is incomplete; if it distributes across contacts
and scales `O(dt²)`, it is γ-phys.

## 3. Pass/fail reasoning (LEMMA-grade bounds)
- **LEMMA DR1 (energy-drift bound).** For the exp-map propagator on the two-robot
  NVE trajectory, `|⟨H(t)−H(0)⟩|` SHALL be bounded and `O(dt²)` — i.e. the drift
  rate `d⟨H⟩/dt → 0` as `dt→0` at fixed physical time, and `⟨H(t)−H(0)⟩` shows no
  linear-in-`t` term over `L` steps. Expected structure: fluctuation amplitude
  `~ √N_dof · dt²` (bounded), NOT secular. A residual secular term that does NOT
  vanish as `dt→0` is γ-pump (integrator bug), not physics.
- **LEMMA DR2 (KE conservation of the total).** Total KE + PE conserved to the
  DR1 bound; any per-body KE source whose integral grows `∝ L` is a defect, not
  elastic exchange (which integrates to ~0 by recurrence).

**Attribution rule.**
- Bounded `O(dt²)` drift, distributed over contacts, vanishing with `dt` → γ-phys.
  The friction is REAL and benign; the fix is transport/policy (`20-...`), not the
  integrator.
- Secular drift, root-localized, `dt`-scaling anomalous, matches Simbody oracle →
  γ-pump. Fix: close the Free-joint `(N, Ndot, qddot)` root cause; until then the
  exp-map is the correct propagator and SHALL NOT be reverted.
- Drift correlated with corrector-warning steps → γ-corr. Fix: lower `dt` / tighten
  corrector (`docs/specs/ncmc-explicit-solvent/20-inner-integrator.md`).

## 4. Touch list (diagnostic only)
- Per-body / per-DOF KE logging over an NVE trajectory (diagnostic, not physics).
- Reuse of `checkReversibility` (INV-REV, `10-...`) and the Simbody differential
  oracle harness; no engine physics change.
