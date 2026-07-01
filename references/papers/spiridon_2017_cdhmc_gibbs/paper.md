# Hamiltonian Monte Carlo with Constrained Molecular Dynamics as Gibbs Sampling

Laurentiu Spiridon, David D. L. Minh. J. Chem. Theory Comput. (2017). DOI: 10.1021/acs.jctc.7b00570.

## Abstract

Compared to fully flexible molecular dynamics, simulations of constrained systems can use larger time steps and focus kinetic energy on soft degrees of freedom. Achieving ergodic sampling from the Boltzmann distribution, however, has proven challenging. Using recent generalizations of the equipartition principle and Fixman potential, the authors implement Hamiltonian Monte Carlo based on constrained molecular dynamics as a Gibbs sampling move. By mixing HMC based on fully flexible and torsional dynamics, they reproduce free energy landscapes of simple model systems and enhance sampling of macrocycles.

## 1 Introduction

Molecular dynamics (MD) simulations of nanoscale systems often become trapped in local minima because the timescales of important conformational transitions exceed those accessible with modern computers. Enhanced-sampling algorithms that do not preserve kinetics but sample configurations from a target distribution (e.g. the canonical Boltzmann distribution) are a response to this.

One way to enhance sampling is to impose holonomic constraints on high-frequency motions that limit the MD time step. Constraints allow larger time steps and focus kinetic energy on slower degrees of freedom. Constraining bond lengths involving hydrogen enables a 1 -> 2 fs boost. Constraining bond lengths and (some/all) bond angles, leaving torsions flexible (torsional MD), is numerically stable up to 10 fs time steps and has been applied to NMR/crystallographic structure refinement and protein structure prediction. Constraining even larger subsets (holding protein domains rigid, allowing hinges to flex) enables observation of large conformational changes.

Enabling advances: constrained MD integrators and compensating-potential methods that scale **linearly** rather than cubically in the number of DOF; an equipartition theorem for internal coordinates (allowing correct velocity assignment at a specified temperature); and software (GneimoSim, Simbody, Molmodel).

Constrained dynamics can be done in Cartesian coordinates enforcing constraints each step (Tao et al. algorithm: arbitrary block size, energy-conserving, reversible, symplectic - but inefficient when constrained blocks share atoms, as in torsional MD), or in internal coordinates fixing specific DOF. The internal-coordinate approach needs fewer DOF and eliminates force-field terms like bond length. The first internal-coordinate MD integrator needed an O(n^3) mass-matrix inversion; Jain et al. gave a recursive O(n) algorithm. But because constrained and unconstrained dynamics sample different probability densities, an O(n) integrator alone does not yield correct Boltzmann sampling.

### Mass metric tensor and probability densities

Distinctions between constrained and unconstrained probability densities are encapsulated in the mass metric tensor. For Cartesian coordinates `q` with linear momenta `p_q`, the kinetic energy is `½ p_q^T M^{-1} p_q` with `M ∈ R^{3N x 3N}` a position-independent diagonal matrix of atomic masses. For `N_φ` generalized coordinates `φ = φ(q)` with momenta `p`, the kinetic energy `½ p^T M_φ^{-1} p` is based on the mass metric tensor `M_φ ∈ R^{N_φ x N_φ}`, defined `M_φ = J^T M J` with Jacobian `J_kl = ∂q_k/∂φ_l`. Unlike Cartesian coordinates, `M_φ` is position-dependent.

In the Boltzmann distribution the joint density is Eq. boltzmann-joint, with Hamiltonian Eq. hamiltonian. Thermodynamic predictions ignore momenta; the marginal over coordinates (Eq. 2) carries a `|M_φ|^{1/2}` factor. If unconstrained there are 3N DOF (`M_{3N}`); if partitioned into `N_f` flexible and `3N - N_f` constrained coords, the flexible-block tensor `M_{N_f}` is a sub-block of `M_{3N}`. Because in general `|M_{3N}| ≠ |M_{N_f}|`, the marginal densities of constrained and unconstrained dynamics differ for the same potential.

### Fixman potential

Fixman proposed simulating with modified potential `U'(φ_f) = U(φ_f) + U_F(φ_f)` (Eq. 3) so the constrained-dynamics marginal matches the unconstrained marginal. Its importance grows with chain length; Echenique et al. concluded it should be considered for peptides longer than 2 residues. General determinant evaluation (e.g. Cholesky) is O(n^3), but an O(n) algorithm for the Fixman potential and its derivative (the **Fixman torque**) for serial robotic systems and general branched molecules exists (Jain et al.).

Alternatives to circumvent cubic Fixman cost: specialized internal-coordinate force fields approximating a Cartesian source FF without the compensating potential (Katritch, Chen); a molecular-mechanics integrator on a modified energy hypersurface approximating Boltzmann without explicit Fixman (Vitalis & Pappu) - but these are restricted to a single predetermined constraint set; explicit free-energy cost of imposing constraints (König & Brooks) - requires extra intermediate thermodynamic states.

### Hamiltonian Monte Carlo

The authors use HMC (Hybrid Monte Carlo), an MCMC sampler based on Metropolis-Hastings that uses MD to generate candidates. Given initial config `φ_i`, velocities are randomized and a trajectory propagated for a set number of steps to a trial `φ_f`, accepted/rejected per Eq. 4, where `T(φ_j|φ_k)` is the attempt probability (including velocity-generation probability). Advantages over pure MD: (1) velocity initialization need not obey exact Maxwell-Boltzmann - only the probability ratio of initial/final velocities is needed; (2) the integrator need not preserve the Boltzmann distribution - only the forward/reverse propagation probability ratio matters (unity for a reversible, symplectic, deterministic integrator). This permits large time steps, constrained integrators, and **distinct guidance vs acceptance potentials**: use unmodified `U(φ)` during dynamics and modified `U'(φ)` only in the acceptance probability, drastically reducing Fixman evaluations.

Forrest & Suter (1994) pioneered constrained-dynamics HMC (CDHMC), simulating polymers via torsional MD. They initialized moves with Maxwell-Boltzmann velocities using fictitious position-independent moments of inertia and excluded the Fixman potential from the propagator. They did not rigorously sample Boltzmann because (1) the acceptance Hamiltonian omitted the Fixman potential (no efficient method existed then) and (2) they used only torsional MD, never relaxing bonds/angles.

### Gibbs sampling to restore ergodicity

Fully rigorous ergodic sampling with only one constraint set is impossible: a constrained simulation samples a conditional `ρ(φ_f|φ_c)`, which only approximates the marginal `∫ dφ_c ρ(φ_f, φ_c)` when the conditional and marginal are similar. Matching them (e.g. Kandel et al. selectively relaxing angle constraints) limits efficiency, needs foreknowledge of correlated DOF, and is not rigorous unless constrained DOF are fully independent from flexible ones - and it cannot be ergodic since it never accesses other values of the constrained DOF.

To sample all configuration space, constraints must be **alternated** so every DOF gets a chance to be flexible - e.g. torsional MD alternated with fully flexible Cartesian MD. Such mixed simulations are a special case of **Gibbs sampling** (a MCMC method generating samples from a series of conditional distributions; replica exchange and expanded-ensemble simulations are also Gibbs sampling). This work is stated to be the first use of constrained dynamics in Gibbs sampling moves for rigorous Boltzmann sampling of molecular systems.

## 2 Methods

### 2.1 Software (Gmolmodel)

The square root of the mass matrix and its determinant were tested against the Eigen 3.2 library. Gmolmodel was developed for use in Alchemical Grid Dock (AlGDock), which calculates binding free energies between flexible ligands and rigid protein conformations (useful in implicit ligand theory). Gmolmodel implements a shared-memory model to send configurations to and obtain forces from the Molecular Modeling Toolkit that AlGDock is based on. Boost.Python (Boost 1.55) provides a Python interface to CDHMC within AlGDock. AlGDock is open-source (MIT), on GitHub (CCBatIIT/AlGDock).

### 2.2 Sampling Algorithms

HMC moves are based on either constrained or unconstrained MD.

Constrained moves used two constraint types: (i) **torsional dynamics** - all bonds and angles constrained, torsions flexible (majority of CDHMC sims); (ii) **rigid body** - only the central torsion of butane flexible. For both, velocities were drawn per the internal-coordinate equipartition principle so that initial and final velocity probabilities are given by the Boltzmann-weighted kinetic energy. Dynamics were propagated on the **unmodified** potential with the Simbody Velocity Verlet integrator, which satisfies constraints on positions and velocities and computes velocity-dependent Coriolis/gyroscopic forces within tolerance η = 10^-4. Unless stated, trajectories were 10 MD steps of 4 fs. The Fixman compensating potential was included (or excluded, for testing) in the Metropolis-Hastings acceptance criterion.

Markov-chain step pseudocode (given current `φ_{t-1}`, produce next `φ_t`):

1. Initialize the trial trajectory.
   - (a) Set initial config `φ_0* = φ_{t-1}`.
   - (b) Draw initial momenta `p_0* ~ N(μ=0, Σ = k_B T · M^{-1}(φ^{t-1}))` (Eq. momenta-draw).
2. Propagate by Velocity Verlet with step size ε. While `n < N`:
   - (a) position predictor (Eq. vv-position)
   - (b) project positions (Eq. vv-position-proj)
   - (c) i. compute velocity-dependent forces `f_p`; ii. momentum update (Eq. vv-momentum); iii. project momenta (Eq. vv-momentum-proj); iv. if the relative projection correction exceeds η (Eq. vv-tolerance), repeat.

   `proj(·)` is a projection satisfying the constraints.
3. Accept or reject the final configuration using the modified Hamiltonian (Eq. accept-pseudocode), where `H'(φ,p) = H(φ,p) + U_F(φ)` (Eq. modified-hamiltonian) and `U ~ U(0,1)`.

Unconstrained HMC moves drew velocities from Maxwell-Boltzmann, used Velocity Verlet with 1.5 fs, and Metropolis-Hastings acceptance. Most moves used 100 dynamics steps, but acceptance-rate comparisons used 10 steps. Unless stated, sims were at 300 K.

Four simulation types combined the two move kinds:

1. **UDHMC** - only unconstrained-MD HMC moves.
2. **CDHMC-noFixman** - constrained MD without Fixman torque; Fixman potential NOT used in acceptance.
3. **CDHMC** - constrained MD without Fixman torque; Fixman potential USED in acceptance.
4. **MIXED** - 15 CDHMC moves alternated with 5 UDHMC moves (for alanine dipeptide, 10 and 1). Constrained moves used no Fixman torque for dynamics but the Fixman potential for acceptance.

## 3 Results and Discussion

Simulations tested (1) reproducing the unconstrained Boltzmann distribution and (2) enhancing sampling. Correctness tested on a serial chain, butane, and alanine dipeptide; efficiency on five macrocycle molecules. (Concrete numbers are in `checks.md`.)

### 3.1 Torsion Angle Distribution of an Idealized Serial Chain

A 4-bead idealized chain (C4) with fixed bonds/angles and no nonbonded or torsional terms should be uniform in torsion space, and its Fixman potential has a known closed form. Results verify that including the Fixman potential is necessary to draw from the correct (uniform) Boltzmann distribution: CDHMC-noFixman shows periodic distortion, while CDHMC and MIXED recover uniformity. This also confirms the constrained-dynamics integrator does **not** require the Fixman torque to sample correctly (consistent with using distinct guidance vs acceptance potentials, Duane et al.). Unconstrained moves are unnecessary here because the torsion is independent of other DOF.

### 3.2 Potential Energies of Butane

Shirts' test: log-ratio of energy histograms at two temperatures is linear in U with slope `-(β_2 - β_1)` (Eqs. 5, 6). CDHMC and MIXED (torsional and rigid-body) all pass; fitted slopes are within a standard deviation of the expected `0.13363` (Table 1 in `checks.md`).

### 3.3 Free Energy Landscape of Alanine Dipeptide

Free energy `F(Φ,Ψ) = -k_B T ln[ρ(Φ,Ψ)]` (Eq. free-energy). UDHMC and MIXED landscapes are qualitatively consistent; MIXED shows lower barriers in two transition regions. Free energies of C5, PPII, α_L relative to C7eq agree within one std between methods (numbers in `checks.md`), evidencing MIXED samples the same distribution as UDHMC.

### 3.4 Conformational Transitions in Alanine Dipeptide

The MIXED empirical transition matrix has smaller diagonal and larger off-diagonal elements than CDHMC, i.e. MIXED undergoes conformational transitions more readily per step. Mean first passage times (MFPT, via PySAL) are shorter for MIXED for every state pair, especially to the isolated α_L. Per integrator step, MIXED is more efficient than UDHMC. (Tables 3, 4 in `checks.md`.)

### 3.5 Sampling Efficiency of Macrocycle Simulations

Five macrocycles (1R6, AA0, AB0, ACZ, ADN) sampled. Constrained dynamics tolerates larger time steps: UDHMC acceptance ~35% at 2 fs and exactly 0% at 3 fs, whereas CDHMC stays >0.1 at 10 fs and >0.7 at 4 fs. At 625 K, MIXED samples a similar or broader configuration space than UDHMC in the same MD-step count; hierarchical clustering (cutoff 0.25 Å) shows MIXED conformations are nearly a superset of UDHMC (MIXED ≈ UNION, UDHMC ≈ INTERSECTION). The 4 fs step accesses clusters faster than 1.5 fs. Constrained dynamics may improve sampling by a mechanism similar to self-guided Langevin dynamics (SGLD): suppressing high-frequency motion and focusing kinetic energy on low-frequency DOF - taken to the extreme, plus larger time steps.

## 4 Conclusions and Future Directions

Feasibility of Boltzmann sampling and improved efficiency via CDHMC as a Gibbs move is demonstrated by: (1) uniform torsion marginal in an analytic model, (2) correct energy-histogram ratio across temperatures, (3) 2D free energy landscape matching flexible dynamics. Efficiency shown via transition matrices/MFPT and via more distinct conformations per MC trial (with fewer energy evaluations).

Constrained dynamics for CDHMC does **not** require the Fixman torque; omitting it gives ~25% speedup. But as Fixman importance grows with system size, neglecting the torque for large systems may steer toward low-acceptance conformations.

Future directions: torsional dynamics may be less efficient at low temperature (1-4 contacts, van der Waals repulsions between atoms three bonds apart); modified torsional terms or softened vdW/electrostatics can help; use modified potentials as guidance Hamiltonian while retaining the atomistic FF for acceptance. New constraint types (rigid multimeric-protein domains) could accelerate large conformational transitions. Constrained dynamics could enable higher-temperature replicas in replica exchange (enhancing relative domain motion rather than unfolding). Tuning of CDHMC:UDHMC ratio and time step; possible use of SGLD integrator or advanced acceptance criteria in CDHMC moves.
