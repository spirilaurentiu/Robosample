# A simple molecular mechanics integrator in mixed rigid body and dihedral angle space

Andreas Vitalis, Rohit V. Pappu. J. Chem. Phys. 141, 034105 (2014). DOI: 10.1063/1.4887339.

## Abstract

We propose a numerical scheme to integrate equations of motion in a mixed space of
rigid-body and dihedral-angle coordinates, targeting biomolecular systems with tree-like
topology. By approximating the effective mass matrix as diagonal and lumping all bias
torques into the time dependencies of the diagonal elements, we exploit the formal
decoupling of the individual equations of motion. Energy conservation is imposed
independently for every degree of freedom, yielding the numerical integration scheme.
All auxiliary operations cost O(N_at). Coupling to two thermostats extends the method to
constant-temperature ensembles. The integrator is stable and free of mass-metric tensor
(MMT) artifacts by construction. Two systems - liquid water and an alpha-helical peptide
in continuum solvent - establish thermodynamic accuracy against reference sampling on
exactly the same degrees of freedom.

## I. Introduction

Molecular mechanics propagates internal/relative molecular motions with force fields,
usually via MD (numerical integration of Newton's equations). Increasing the integration
time step reduces force-evaluation cost but amplifies discretization error and eventually
causes instability, especially for low-mass particles under stiff potentials. Two common
mitigations: redistribute/change atomic masses, or eliminate the fastest motions via
constrained dynamics (Lagrange-multiplier solvers like SHAKE, or explicit generalized
coordinate spaces). This paper concerns the generalized-coordinate approach.

Some important force fields (ECEPP, ABSINTH, ROSETTA) assume dihedral angles as the only
internal DOF (fixed bond lengths and angles). A stable MD engine producing unbiased
equilibrium ensembles for such force fields is of continued interest.

Consider N_at atoms in N_mol molecules with position vectors r_i in a space-fixed global
frame, forming the 3N_at Cartesian state vector r, with separable Hamiltonian and
potential U(r). The Cartesian Hamilton equations are eq:1. Because M is diagonal, the
3N_at velocity equations are uncoupled; any coordinate can be constrained by zeroing its
momentum without changing others' kinetic contributions, so the momentum-integration
factor in the partition function is constant for a given constraint set and temperature.

In generalized coordinates phi (dimension 3N_at), the kinetic energy is
E_k = (1/2) p_phi^T G^{-1} p_phi = (1/2) omega^T G omega, where G = J^T M J is the
mass-metric tensor (MMT) and J is the Cartesian->generalized Jacobian (eq:2, eq:3).
For a separable Hamiltonian the canonical partition function is
Q = int sqrt(det G) exp[-beta U(phi)] dphi; without constraints det G is independent of
phi and Q factorizes into thermal and configurational parts.

The generalized-momenta EOM (eq:2) reveal bias terms from the time derivative of J
(rotating-frame torques). Formally imposing constraints (zeroing a subset of omega)
leaves a sqrt(det G_S) term whose det(G_S) **varies with phi** (unlike det G) - these are
the **MMT artifacts**. Fixman-style compensating potentials containing det(G_S) correct
this but are not routinely used and can be costly. Recursive O(N_at) computation of
det(G_S) or its derivatives exists (treating each molecule as a chain of rigid bodies with
hinges). Elements of J for angular DOF are eq:3.

### The approach of this paper

This method neither ignores MMT artifacts nor corrects them with a Fixman potential. We
never evaluate det(G_S), det(G), or their derivatives. Instead we derive propagators for
phi and omega corresponding to an (artificial) dynamics on a modified constant-energy
hypersurface with Omega_k = I_kk^{1/2} omega_k and phi_k as independent variables, where
I_kk are the diagonal MMT elements and the modified kinetic energy is
E'_k = (1/2) sum_k I_kk omega_k^2. This restores the simple structure of eq:1 while lumping
all bias torques into a time-dependence of effective masses. Given r and grad U(r), the
I_kk and F_phi are computed recursively in linear time. Because the kinetic energy is
modified and dynamic coupling reduced, the dynamics are altered for all but special cases -
this is precisely what avoids MMT artifacts by construction. Design goals: simplicity,
versatility, thermodynamic (not kinetic) correctness; the method is in the spirit of MC
propagators.

## II. Methods

### A. Motivation of approximations

Goals: (1) avoid MMT artifacts (comparability to MC in the same coordinate space);
(2) O(N_at) auxiliary computations, easy implementation; (3) stable numerical integration;
(4) seamless multi-molecule support; (5) thermodynamic (not kinetic) correctness.

The system is represented in mixed rigid-body + dihedral-angle coordinates of dimension
K <= N_at. Dihedrals are a subset of a Z-matrix (M-1 bond lengths, M-2 bond angles, M-3
dihedrals). Bond lengths and angles and some dihedrals are constraints. Rigid-body
coordinates of each molecule are encoded via Cartesian coordinates of its first three
atoms (Appendix A). A reversible mapping is stipulated: eq:4, phi = A(r), r = A^{-1}(phi),
with A^{-1} imposing a tree structure (build from first three atoms using Z-matrix
variables, backward dependency along the main chain and to branch points). Z-matrix
variables are convenient but limit flexible rings (one ring bond is not represented).

The key approximation is a diagonal mass matrix (eq:5): the modified KE uses I_D and
generalized velocities, differing from the standard KE. This removes kinetic coupling,
formally decouples the EOM (structurally like unconstrained Cartesian dynamics), and
restores the equipartition principle to its intuitive form. Units across phi are
heterogeneous. Effective masses are actual masses (translation) and rotational inertias
(rotation), arranged as eq:6. The elements I_kk, F_phi, and estimates for dot(r) come from
the recursions eqs:B1-B3.

### B. Numerical integration scheme on a constant energy hypersurface

Mandate conservation of E'_k(omega) + U(r) over a discrete step t1 -> t2 (eq:7). This
single condition is insufficient; because I_D is diagonal (analogy to the Cartesian case),
we require eq:7 to hold for **each of the K DOF independently**, giving a per-DOF quadratic
whose solution is eq:8. Eqs:7,8 carry a time dependence of I_kk that preserves kinetic
energy for angular variables under the diagonal assumption - fundamentally different from
how true bias torques would preserve angular momentum and total energy with a non-diagonal
mass matrix. The initial choice of omega defines the explored hypersurface.

Eq:8 is time-reversible but not guaranteed to have a real solution. Approximating I_kk(t2)
in the second square-root term by sqrt(I_kk(t1) I_kk(t2)) leaves one meaningful root
(eq:9); eq:9 (not time-reversible) is used as a guide to pick the correct solution of eq:8.
Positions update via eq:10 (simple leapfrog increment); the exception is rigid-body
rotation, handled by a quaternion (eq:12) with no explicit phi_k.

Eq:8 is also implicit (I_kk(t2) needed for new velocities). The assumed time-dependence
I_kk(t) masks the true dependence I_kk({phi_{i!=k}(t)}). Required I_kk(t2) values are
obtained by guessing phi_k(t2), not by extrapolating I_kk. Procedure: obtain a guess of
omega_k(t_{1.5}) using eq:8 with t2 -> t_{1.5} and delta_t/2, increment positions by
delta_t/2 to phi_k(t2), compute I_kk(t2) guesses, then restore the configuration and apply
eq:8 as written. The scheme is force-explicit (no extra force evaluations).

With time-independent masses, eq:8 recovers the standard leapfrog integrator; integrator
behavior depends on the rate of change of I_kk (magnitude of effective bias torques). If
I_kk changes rapidly, velocities are updated iteratively via eq:11: the step is partitioned
into Lambda segments assuming linear evolution of I_kk between t1, t_{1.5}, t2. The force is
held constant (force-explicit); eq:11 is not time-reversible; the eq:9-analog picks each
sub-step's root. Lambda should be a multiple of 2; benefits taper for large Lambda. Eq:11
better represents changing effective masses (e.g., rigid water rotation).

If hidden computations are recursive, complexity is O(N_at), less than the cost of U(r) for
all but trivial cases. Limiting cases (monoatomic gas, molecular liquids, polymer mixtures)
are handled uniformly.

### C. Technical issues

Rigid-body coordinates appear explicitly in the EOM (difference from the spatial-operator
formalism, which maps external motion to the base body). A rule assigns the base for each
flexible molecule (terminus or chain middle tested here). Rigid-body rotation uses the
quaternion update eq:12: omega_{x/y/z} are angular velocities at t2 about fixed lab axes
through the molecular COM, c set by unit-norm. Rigid translation of each COM is
straightforward (constant masses). For flexible molecules the COM is updated after the full
conformational update, so displacement is mismatched relative to the velocity increment;
linear momentum is conserved (no external forces) under eq:5, but sum_k p_k is conserved
only for rigid molecules.

Thermostats (temperature is a function of omega^T I_D omega):
- **Andersen**: stochastic per-DOF coupling; velocities reset from the target Maxwell
  distribution using the current-conformation I_D. Applied right after computing
  F_phi(t_{1.5}), before velocity increments. Asynchronous/independent coupling makes
  equipartition artifacts unlikely but may slow collective dynamics.
- **Velocity rescaling (Bussi et al.)**: global rescaling from a single stochastic process.
  In eqs:8,9,11 replace omega_k(t1) with alpha_T omega_k(t1), where alpha_T is the global
  rescaling factor from instantaneous vs target temperatures.
Both use a coupling time tau_T. Full cycle in Appendix C.

### D. Underlying equations of motion (justification of artifact-freedom)

The EOM integrated by eqs:8/11 are NOT Hamiltonian (eq:2) because of the approximation
eq:5. The thermostats enforce per-DOF equipartition eq:13. The bias torques (visible in
eq:9) conserve omega_k I_kk^{1/2}, so neither omega_k nor omega_k I_kk are independently
distributed at zero potential.

### Derivation (not implemented)

Starting from eq:8 for DOF k (eq:14), substitute Omega_k(t) = omega_k(t) I_kk(t)^{1/2}
(eq:15); in the linear approximation the sums are twice the half-step means and the finite
difference is the rate of change; letting delta_t -> 0 gives eq:16. Excluding the trivial
omega=0 case yields the underlying EOM eq:17, showing bias torques are hidden by treating
Omega as the dynamical variable; the I_kk^{1/2} factors preserve volume in phi. With
Lagrangian L = sum_k (1/2) Omega_k^2 - U, the first line of eq:17 is eq:18. A further
substitution Phi_k = phi_k I_k^{1/2} fails for non-constant masses (then dot(Phi_k) != Omega_k),
so eq:17 is generally artificial dynamics. Writing the canonical partition function in
independent variables Omega, phi (eq:19) and integrating momenta gives eq:20 - the pure
Boltzmann configurational integral times a temperature-only prefactor, i.e. no det(G)
weighting. Eq:20 hides that the method does not transform Cartesian momenta canonically
(eq:5).

### E. Simulation protocols

All simulations used CAMPARI (http://campari.sourceforge.net). See `checks.md` for the
full parameter fixtures. Summary:

1. **Flatness test:** 18-atom PEG-like linear polymer (15 dihedral + 6 rigid-body DOF),
   eq:11 with Lambda=4, dt=5 fs, three mass distributions, base at atoms 8-10 or 1-3,
   Andersen tau_T=1 ps. Reference: Langevin + SHAKE (Skeel-Izaguirre impulse integrator,
   dt=5 fs, friction 1 ps^-1), constraints on bonds+angles (leaving 15 dihedrals) or bonds
   only. 50 runs x 10 ns each.
2. **Integrator stability:** two capped (GS)_50 chains, cubic 200 A box, PBC; U = 12th-power
   repulsion (cutoff 10 A) + amide-planarity dihedral potentials. eq:11, Lambda=4. dt 4-10 fs
   (correct masses) / 10-30 fs (adjusted masses). 20 runs x 1 ns.
3. **Rigid water:** 1095 TIP4P, cubic 32 A box, PBC, 12 A cutoff, reaction-field
   electrostatics, velocity-rescaling tau_T=1 ps; reference Cartesian leapfrog + SETTLE.
4. **FS peptide:** N-Acetyl-A5(AAARA)3A-N'-methylamide, ABSINTH implicit solvent, 40 A
   droplet, ~0.15 M NaCl, half-harmonic boundary (0.05 kcal/mol/A^2), 12 A cutoff,
   backbone-phi blocking potential, Andersen tau_T=10 ps, eq:11 Lambda=4. Reference REMC.
   Helicity from torsional segment statistics: >=2 consecutive alpha-basin residues = a
   segment (N_s); a length-N_alpha segment gives N_alpha-2 H-bonds (N_h); length-1 runs -> N_1.

## III. Results (summary)

- **A. Absence of MMT artifacts:** with U=0, eq:11 gives flat dihedral histograms for all
  15 DOF across three mass distributions and both bases (Fig. 1a-c,f); the canonically
  transformed reference (Fig. 1d) shows the expected artifact. Confirms the Sec. II D
  reasoning; artifacts avoided by construction.
- **B. Integrator stability:** (GS)_50 stable to dt~6 fs (correct masses, ~340 K); serine
  chi_2 dominates error; redistributing hydroxyl masses to 8.5 Da each raises stability to
  dt~16 fs. Integrator error is stochastic (rare-event jumps, not linear drift).
- **C. Liquid water:** matches SETTLE reference in energetics, pair-correlation functions,
  C_v, diffusion, dielectric, rotational correlation (Tables I, II; see `checks.md`). At
  dt=2 fs Lambda in {1,2,4} equivalent; at dt=5 fs stability improves with Lambda. Small
  <T> vs <T_c> mismatch (<=0.2%) is tied to integrator error and should be monitored. Not as
  stable as dedicated rigid-body integrators, but thermodynamically accurate.
- **D. FS peptide:** helix-coil transition (N_h, N_s, N_1 vs T) overlaps the MC/REMC
  reference within error; same melting temperature. Ensembles thermodynamically identical
  across N/M/C-base choices but kinetically distinct (base controls transition rates;
  M-base gives the most natural dynamics). Equipartition issues (flying-ice-cube for weakly
  coupled rigid-body/chi DOF; Cl- vs Na+ imbalance depending on base) shrink as dt
  decreases. <omega^T I_D omega> ~ <p^T M^{-1} p> in the mean (eq:22 condition) but
  fluctuations differ (eq:21 inequality); only the modified KE reproduces the ideal K/2 k_B.

## IV. Discussion and conclusions

Frozen bond lengths/angles bound the I_kk, so bias torques are an unlikely instability
source. The reasoning in Sec. II D requires approximate equipartition (eq:13). The scheme
avoids systematic biases in phi (eq:20) at the cost of not preserving total angular
momentum, introducing artificial dynamics, and not preserving Cartesian-momentum phase-space
volume. Approximate agreement <omega^T I_D omega> ~ <p^T M^{-1} p> is restated as eq:22
(vanishing averaged off-diagonal MMT contribution). Dihedral velocities are expected to be
less cross-correlated than, e.g., interatomic-distance velocities, favoring eq:22.

The method performs MD in mixed dihedral + rigid-body space with a diagonal-mass
approximation (eq:5) giving Cartesian-like structure and O(N_at) recursions (Appendix B).
It is thermodynamically accurate on challenging systems vs holonomic-constrained MD and MC
references, handles any polymer/small-molecule mixture uniformly (except flexible rings),
and adds negligible cost over Cartesian force evaluation. Loss of accurate dynamics is the
main caveat. A prior simplified variant (I_kk lagging by delta_t/2) was combined with MC
splicing; hybrid MC + this dynamics is a promising direction (an MC propagator that "jumps"
in phi would obviate the blocking potential). Ongoing work: a proper Langevin integrator
and a unified internal-coordinate sampling engine.

## Appendix A: Coordinate operations

Cartesian coordinates are built with a backward dependency of r_i on three previously built
reference atoms r_j, r_l, r_m; the first three atoms of each molecule are stored explicitly
(rigid relative orientation). Placement uses eq:A1 (Z-matrix -> Cartesian, NeRF-style),
applied hierarchically to form A^{-1}. If the dihedral is rotatable, its phi_k is a DOF and
a_k is its unit bond vector; r_j, r_l, r_m must be non-collinear.

Different building directions (bases of motion) can use a single A^{-1} plus compensatory
rotations of the first three atoms; the quaternions for these rotations follow from the
eq:10 increments of the dihedrals in question.

J is described differentially; Y is the reduced matrix of covariant base vectors for
flexible DOF (3N_at x K, eq:A2), identical in form to eq:3. Y is adjusted for different
bases by the sign of a_k and which terms vanish. Instantaneous Cartesian velocities from
eq:A3 (sum over DOF base-ward of atom i, incl. parent branches and rigid-body motion) give
the true KE (1/2) p^T M^{-1} p. For rigid rotation the reference frame uses the COM (b_k)
and lab-frame axes (a_k).

## Appendix B: Recursion formulas

I_D and F_phi are computed in the same inward (tip->base) recursion; only the angular case
is given. Projected force: eq:B1 (axis-projected net torque about b_k), sums over all
tip-ward atoms incl. sub-branches; branch values combine at merge points; rigid rotation is
the last step (sum over all atoms). Effective inertia: eq:B2 (last term = Frobenius inner
product of outer-product matrices). Cartesian atomic velocities via outward (base->tip)
recursion: eq:B3, sums over base-ward DOF. All are O(N_at) because the sums contain no
DOF-specific terms and can be accumulated.

## Appendix C: Numerical implementation (integration cycle)

Assume the conformation at t_{1.5} (current), I_D at t_{0.5}, velocities at t_1, and a prior
guess of I_D at t_1. At most one U(r)/grad evaluation per step. Cycle:

1. Store I_D at t_{0.5}; initialize new F_phi, I_D.
2. Compute U(r), grad U(r), F_phi, I_D at t_{1.5} (latter two via recursions).
3. Compute dot(r) using Y at t_{1.5} and omega at t_1.
4. Thermostat: (a) velocity rescaling - infer current T from omega, I_D at t_1; derive
   alpha_T. (b) Andersen - for each DOF draw a uniform [0,1] number, compare to delta_t/tau_T;
   if smaller, reassign omega_k(t1) from a pseudo-Boltzmann distribution using I_kk at
   t_{1.5}; set alpha_T = 1.
5. Back up coordinates at t_{1.5}.
6. Apply eq:8 with omega_k(t1)=alpha_T omega_k(t1), delta_t=delta_t/2, I_kk(t2)=I_kk(t_{1.5});
   compute increments directly as phi_k(t2)-phi_k(t_{1.5}) = 0.5 delta_t omega_k(t_{1.5}).
7. Update coordinates; back up prior I_D guess at t_1; compute I_D guess at t_2.
8. Restore coordinates to t_{1.5}.
9. Velocity iteration eq:11 (omega_k(t1)=alpha_T omega_k(t1)) using I_D at t_1, t_{1.5}, t_2 -> omega at t_2.
10. Increments from eq:10; update coordinates to t_{2.5}.
11. Update each molecule's COM (rigid-translation omega unchanged).
12. Accumulate properties (dynamics at t_2, structure at t_{2.5}).
13. Return to step 1.

Coordinate updates (steps 7, 10) include the pre-rotation of the 3 reference atoms per
molecule when the base of motion differs from the natural A^{-1} structure. Timing (single
core Xeon E5410): auxiliary O(N_at) recursions add negligible cost vs force evaluation (see
`checks.md`).
