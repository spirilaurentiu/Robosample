# Two-robot contact world — Diagnosis and reformulation

## 1. Problem restatement

**User phrasing.** "A single robot accepts large HMC timesteps. TWO robots in
contact collapse the timestep ~25x due to intermolecular steric clashes, so the
method loses its edge over Cartesian MD/HMC. Trajectories go from ballistic to
diffusive ('friction/dampening'); NCMC has not rescued it. Fix it with a
momentum-preserving mixed torsional/Cartesian contact world."

**Codebase restatement.** Robosample samples the target molecule in the
articulated-body (torsional) generalized coordinates `φ` (`state_.q/u`, mass
metric `M_φ` from the Jain O(n) recursion) via blocked-Gibbs + NVE-Verlet HMC:
draw `u ~ N(0, RT·M_φ⁻¹)` (`reinitialize`, `multiplyBySqrtMInv`,
`World.cpp:1686`), run `L` fixed-`dt` velocity-Verlet steps
(`RobotEngine::stepTo`, `World.cpp:1457`), accept once on
`ΔH = H_end − H_start` with `H = pe + ke + keSolvent + fixman − ½RT ln sin²γ₂
− nmaCorr` (`currentTotalEnergy`, `World.cpp:1804`). Fixman is acceptance-only
(a pure function of `φ`, not a dynamical torque; `calcFixman`, `World.cpp:851`).

The **timestep collapse** is genuine and unavoidable in *any* representation: the
intermolecular LJ `r⁻¹²` wall is a stiff, non-covalent mode present in both
torsional and Cartesian pictures. Torsional coordinates buy a large `dt` by
*constraining out intramolecular* stretch/bend modes; they do **not** constrain
the *inter*molecular contact mode, so a clash raises `ω_max` and the Verlet
stability bound `dt < 2/ω_max` drops `dt` ~25x for both representations. The
torsional method therefore keeps its `dt` edge only *away* from contact. The
question is not "recover the large `dt`" (impossible at a stiff contact) but:
**at the collapsed `dt`, does the torsional-mixed move still decorrelate the slow
interface dihedral in fewer force evaluations than Cartesian MD/HMC?**

The proposed method is a **mixed contact world**: the target robot in `φ`
(ABA), the second robot as a flat-Cartesian **environment** `E`, coupled through
ONE shared OpenMM force evaluation per step, run as a single NVE trajectory of
`L` steps with momenta persisted across the trajectory and one Metropolis test
at the end. The scaffolding exists: `RobotEngine::verletStep`
(`include/RobotIntegrator.hpp:125-276`) already advances a flagged Cartesian atom
subset (`cartSolventAtoms`) in flat velocity-Verlet on the same shared force,
certified symmetric/reversible/volume-preserving; the `L`-then-one-MH loop is
`World.cpp:1457-1485`.

### 1.1 Two readings, kept separate
- **(R-transport)** Even at the collapsed `dt`, torsional generalized-momentum
  transport of the slow dihedral is *ballistic* (`⟨Δξ²⟩ ~ t²`) while Cartesian is
  *diffusive* (`~ t`). If the ballistic mean free path exceeds the barrier width,
  the mixed move wins per force-eval. This is the make-or-break line
  (`20-...`). Decides whether the method *can* win.
- **(R-defect)** The observed "friction" is not (only) physical thermalization
  but a numerical momentum sink — the documented Free/Ball root-joint KE pump
  (`RobotIntegrator.hpp:216-227`, ~+1300 kJ/mol/traj) and/or the non-converged
  velocity corrector (`RobotIntegrator.hpp:400-407`), plus the *policy* artifact
  of full momentum resampling at short `L·dt`. Decides whether the friction is
  fixable (`30-...`).

Both are pursued; they are distinguished by the diagnostics in `30-...`, not
assumed.

## 2. Binding vocabulary (defined once; reused across files)

- **Target robot R** — the molecule integrated in generalized coordinates `φ`
  (kinematic tree of rigid bodies / mobilizers). Speeds `u`; mass metric
  `M_φ(φ) = Jᵀ M J` (Spiridon & Minh 2017, `references/papers/spiridon_2017_cdhmc_gibbs`,
  eq:mass-metric-tensor), assembled implicitly by the ABA recursion.
- **Cartesian environment E** — the second robot, integrated as free Cartesian
  atoms via the existing `cartSolventAtoms` mechanism (`setCartesianSolvent`,
  `World.cpp:1559`). Positions `x_E`, velocities `v_E`, diagonal **constant**
  mass matrix `M_E`. NOTE: `E` reuses the "solvent" plumbing but is a *robot*,
  not water; the term "environment" is used to avoid the solvent connotation.
- **Mixed configuration `(φ, x_E)`** with joint mass metric `M_mixed`
  (derived block-diagonal in `10-...`).
- **Slow collective coordinate `ξ`** — the interface dihedral(s) whose transition
  is the sampling target (e.g. the clashing torsion; the relative-orientation
  angle). Its generalized speed `u_ξ` is a component/projection of `u`
  (`state_.u()`); its value `ξ` is a dihedral of the atom positions. Both are
  engine state, logged per step.
- **`H`** — the acceptance Hamiltonian exactly as `currentTotalEnergy` assembles
  it (`World.cpp:1804`). `ke = ½ uᵀ M_φ u` (`calcKineticEnergy`),
  `keSolvent = ½ Σ_E m|v|²` (`calcSolventKE`, `World.cpp:1594`),
  `fixman = U_F(φ)` (`calcFixman`).
- **`U_F(φ)`** — Fixman compensating potential,
  `½RT·(lnDetM − lnDetZ − lnDetMCartesian_)` (`World.cpp:880`); Spiridon & Minh
  2017 eq:3, `U_F = ½β⁻¹ ln(|M_{N_f}|/|M_{3N}|)`. Pure function of `φ`.
- **Momentum correlation time `τ_p`**, **ballistic mean free path `ℓ`**,
  **MSAD exponent `α`** — transport observables defined in `20-...` §2.

## 3. Scope and non-goals
- IN: the mixed-integrator correctness (`10-...`), the transport observable and
  momentum policy that decide viability (`20-...`), the physical-vs-numerical
  friction diagnostic (`30-...`).
- SECONDARY: local soft-core NCMC as a complementary move (`40-...`).
- OUT: recovering the pre-contact `dt` (physically impossible at a stiff
  contact); the stiff/rigid *Hessian*-determinant correction (Echenique 2006,
  `references/papers/echenique_2006_stiff_rigid_constraints`) — Robosample
  corrects the mass-metric determinant only, a pre-existing single-robot
  approximation the two-robot design neither worsens nor is responsible for.

## 4. Testbed
TWO ROBOTS IN IMPLICIT SOLVENT (GB/OBC). GB/OBC is a conservative potential of
mean force: no friction, no random force, no explicit bath. Therefore coherent
torsional momentum has only three sinks — (i) the finite, low-dimensional set of
intra-R and R–E DOF (elastic, recurrent), (ii) integrator defects (KE pump,
non-converged corrector), (iii) the full-refresh policy. This is what makes the
testbed decisive: it removes real thermalization as a confound, so a *strong*
friction signal is prima facie numerical (R-defect), and it isolates the
transport question (R-transport) cleanly.
