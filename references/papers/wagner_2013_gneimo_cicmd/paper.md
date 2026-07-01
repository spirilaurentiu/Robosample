# Advanced Techniques for Constrained Internal Coordinate Molecular Dynamics

Wagner, Balaraman, Niesen, Larsen, Jain, Vaidehi. *J. Comput. Chem.* 2013, 34, 904-914. DOI: 10.1002/jcc.23200

## Abstract

Internal coordinate molecular dynamics (ICMD) describes a protein with bond, angle,
and torsional coordinates rather than Cartesian coordinates. Freezing high-frequency
bonds and angles gives constrained ICMD (CICMD). This paper develops a framework for:
(1) initializing velocities for nonindependent CICMD coordinates, (2) efficient
center-of-mass (CM) velocity computation during CICMD, (3) advanced integrators
(Runge-Kutta, Lobatto, adaptive CVODE), and (4) cancelling the "flying ice cube"
effect that arises in Nose-Hoover dynamics. The Generalized Newton-Euler Inverse Mass
Operator (GNEIMO) method is the CICMD implementation. GNEIMO allows a hierarchy of
coarse-grained models by rigidly constraining any group of atoms, including
adaptive on-the-fly "freezing and thawing" of degrees of freedom guided by secondary
structure, which is shown to fold four proteins to native topologies.

## Introduction

Constraints in MD eliminate high-frequency DOF and allow larger timesteps. All-atom
constraint methods include SHAKE and RATTLE (solve Cartesian equations of motion then
iteratively enforce bond-length constraints). ICMD instead uses bond/angle/torsional
(BAT) coordinates. "Torsional molecular dynamics" is the special case of CICMD that
freezes bond lengths and angles; the molecule becomes a collection of rigid clusters
connected by hinges, each cluster a rigid set of atoms with frozen internal bonds/angles.

CICMD equations of motion are more complex: fewer DOF, but the mass matrix is dense and
configuration-dependent. Naive approaches scale as O(N^3). The authors developed a
spatial operator algebra (SOA)-based Generalized Newton-Euler Inverse Mass Operator
(GNEIMO) method solving the same equations of motion exactly at O(N) cost, N = number
of DOF.

CICMD models fold small proteins faster/more reliably than Cartesian models, refine
homology models more accurately, and reproduce large-scale domain motion. Advantages go
beyond larger timestep: CICMD has ~1 order of magnitude fewer DOF and retains essential
BAT degrees of freedom.

Outstanding technical issues addressed here:

- The traditional equipartition principle does **not** hold for CICMD models (atom
  velocities within a cluster are fully correlated). A companion work derived a rigorous
  equipartition principle using new *modal velocity coordinates*; this paper gives a
  low-cost velocity-initialization procedure based on it.
- MD is initialized with zero total linear/rotational momentum, but drift accumulates;
  a composite-body-inertia technique is developed to efficiently null linear+angular CM
  motion.
- A "flying ice cube" behavior emerges in CICMD (thermal energy bleeds from configuration
  DOF); a mathematical explanation and avoidance technique are given.
- Holonomic constraints make the Hamiltonian **nonseparable** (kinetic energy depends on
  configuration), so Verlet is unsuitable; RK4 and Lobatto integrators are analyzed.
- The GNEIMO platform lets the user freeze/thaw any DOF to run all-atom MD, all-torsion
  MD, or arbitrary coarse-graining, including adaptive on-the-fly clustering.

## Method Development

GNEIMO supports: constant-energy (N,V,E) and Nose-Hoover (N,V,T) ensembles; generalized
Born (GB/SA) implicit solvation; multiple molecules including explicit solvent; the Fixman
correction potential; temperature-REXMD with Metropolis switching; periodic boundary
conditions; LAMMPS architectural integration (ICMD coordinate files, MPI-parallel force
fields); GPU-accelerated OpenMM force fields; and soft atom-pair restraints.

### Initialization of velocities in modal coordinates

Focus: torsional MD for tree topologies with 1-DOF internal hinges. Constrained-dynamics
DOF are coupled, so the Boltzmann distribution cannot directly assign generalized
velocities. Instead, identify independent *modal velocity* DOF.

For temperature T, the system thermal energy with N DOF is given by eq:1. Using the SOA
mass-matrix factorization the kinetic energy is eq:2. Define modal velocity coordinates v
by eq:3 (a reversible velocity transformation); recover generalized velocities by eq:4.
Substituting into eq:2 gives eq:5: the kinetic energy is a plain sum of squares over the
independent modal coordinates, so equipartition holds mode-by-mode.

**Algorithm (velocity initialization):**

1. For desired temperature T, use eq:1 to get target kinetic energy Re.
2. Draw modal velocities v from a zero-mean, unit-variance normal distribution
   (except the six DOF of the base clusters).
3. Compute generalized velocities from v via eq:4 (an O(N) recursive base-to-tips scatter).
4. Reset any nonzero CM velocity (method below); this initializes the base cluster's
   6-DOF hinge velocity.
5. Compute total kinetic energy and rescale all velocities to match the target T.

### Resetting the CM velocity in spatial coordinates

Notation follows Jain (ref 27). Simulations start with zero translational/rotational CM
velocity, but numerical error accumulates nonzero linear/angular CM velocity that must be
periodically reset. Computing V_CM (eq:6) needs the system spatial inertia M_S and system
spatial momentum h_S; both depend on atomic positions, momentum also on velocities. Multiple
chains are treated as one system.

**Computing system spatial inertia M_S:** the 6x6 system spatial inertia referenced to the
base cluster equals the composite-rigid-body inertia R(n) (eq:7), using the identities of
eq:7b. Here phi(j,k) is the 6x6 rigid-body transform between clusters j and k, M(k) is the
spatial inertia of cluster k, R(k) is the composite inertia of cluster k and all its
children, phi-tilde = phi - I, and E is the base pick-off operator. The first moment of
R(n) gives the CM location relative to the base-cluster frame.

**Computing system spatial momentum h_S:** eq:hS gives the base-cluster-referenced momentum;
for the 6-DOF base hinge H*(n)=I and theta-dot(n)=V(n), so it splits as eq:hS2.

**Computing V_CM:** combining, eq:8 restates eq:6 with M_S = R(n); the explicit isolated-molecule
CM spatial velocity is eq:9. Adding delta_V to the base-cluster spatial velocity V(n) adds
R(n)*delta_V of spatial momentum; to zero the momentum, eq:deltaV gives delta_V = -V_CM. So
applying an extra base-cluster spatial velocity of -V_CM nulls and resets the system's spatial
momentum. Tested for single-chain; multi-chain tests are future work.

### Flying ice cube effect in Nose-Hoover dynamics

The *flying ice cube* effect drains energy from high-frequency modes into zero-frequency
(translational/rotational) and low-frequency modes, until the molecule drifts rigidly in a
frozen conformation. Known in velocity-rescaling thermostats (Berendsen); here shown to also
arise in Nose-Hoover. Derivation in Cartesian coordinates:

### Derivation (not implemented)

The mass-weighted Nose-Hoover velocity equation (eq:10) summed over all particles, using
sum_i F_i = 0, yields the CM-momentum evolution eq:11. Integrating gives the CM velocity
eq:12; since the exponent is -ln(s_t), a **negative** ln(s_t) makes V_CM grow, producing the
flying-ice-cube instability. The CM kinetic energy follows as eq:13, with c a fit constant
proportional to initial CM kinetic energy.

Empirically (Figures 1-2): CM kinetic energy growth is lower with GB/SA solvation than in
vacuum, and far larger in the cluster model than in all-atom Cartesian - attributable to the
larger cluster-model timestep (20 fs vs 1-2 fs). At 20 fs, KE_CM grew ~10x vs 1 fs (measured
at 90 ps). Mitigation: periodically reset CM velocity using the method above.

### Coarse-graining methods

Three freeze/thaw schemes:

- **Automated clustering model.** A basic cluster-model file defines fixed DOF. Default: every
  terminal atom is merged into a rigid cluster with its nonterminal neighbor, yielding dynamics
  over all torsions. Peptide-bond dihedrals are not constrained; side chains have all torsions
  free except closed rings.
- **Manual freeze and thaw.** The user freezes/thaws DOF beyond the basic model at any time via
  the Python interface (e.g., treat helices as rigid clusters connected by flexible loops).
- **Dynamic coarse graining ("Dynamic Clustering").** Automated hierarchical clustering during
  MD. At a user-defined frequency a STRIDE scan detects secondary structure; detected helix/strand
  residues are mapped to their backbone-torsion clusters and locked. Tied to GNEIMO-REXMD:
  locking can be gated by REXMD temperature threshold and by an angular-velocity upper threshold;
  stress forces at frozen hinges can be monitored to unlock. Locking method from ref 31.

## Results

### Protein structure prediction using the dynamic clustering algorithm

Fold four proteins from extended structures. Secondary structure predicted by PSIPRED; helical
regions built and treated as clusters. GNEIMO-REXMD: 12 replicas, T 300-1050 K, random switching
every 7.5 ps. Helices/strands detected by STRIDE during exchange are frozen if REXMD T > 400 K.
Adaptive-timestep CVODE (Adams-Moulton) integrator used.

Test proteins (resolved subrange): 1BDD (11-56), 1EON (7-31), 1PRB (11-53), 1UBQ (1-35). Starting
structure contains only predicted secondary structure, otherwise extended. Population peaks in
backbone CRMSD: 1BDD 5-7 Å, 1UBQ 6-10 Å, 1EON 7-8 Å, 1PRB 8-10 Å. Closest-to-crystal backbone
CRMSD: 1BDD 4.007 Å, 1EON 4.198 Å, 1PRB 3.726 Å, 1UBQ 4.325 Å.

12 replicas each. 1PRB, 1UBQ: 3 ns/replica (36 ns total). 1BDD, 1EON: 0.3 ns/replica (3.6 ns
total) because velocity reinitialization at each freeze/thaw exchange caused early loss of
secondary structure with longer runs. 1BDD folding path: 12-16 Å broad interhelical contacts ->
8-11 Å incorrect 3-helix packings -> < 7 Å correct native topology, within ~40 ns total. REXMD +
dynamic clustering is for structure prediction, not folding-pathway analysis.

### MD simulations of crystal structures of proteins

Three high-resolution crystal structures: Crambin (1CRN, 1.50 Å), Defensin (1DFN, 1.9 Å), BPTI
(4PTI, 1.50 Å). GNEIMO (N,V,T) all-torsion dynamics 5 ns at 310 K, after 500 ps equilibration.
Hoover bath relaxation 250 fs. Timesteps 1-30 fs, integrators RK4 and Lobatto.

**Integrator performance.** Lobatto = implicit Lobatto IIIa-b partitioned RK (adaptation of
explicit Stormer-Verlet symplectic method); RK4 = standard 4th-order explicit. Lobatto is 2nd
order and computes thermostat/Coriolis forces twice per step but needs only ONE expensive
position-dependent force-field evaluation per step; RK4 makes FOUR force evaluations per step. To
compare fairly, plot metrics against normalized timestep (fs per force computation). RK4 stable up
to 16 fs raw = 4 fs normalized; Lobatto stable up to 10 fs raw = 10 fs normalized. RK4
temperature-STD flat for normalized timesteps < 4 fs; Lobatto error rises with timestep, failing
above normalized ~9-10 fs. Cartesian all-atom unstable above 2 fs. A slow long-term increase in
bath potential energy was observed (attributed to dissipative nonconserving integrators); future
work: multiple Nose-Hoover chains, energy-conserving integrators, Nose-Poincare.

**Structural properties.** Mean backbone CRMSD < 2.5 Å for most torsional runs, uncorrelated with
timestep for Lobatto and RK4. Cartesian crambin and BPTI drifted from crystal; Cartesian 1DFN
stayed near folded.

**RMSF.** Per-residue RMSF (deviation from time-average position) compared to crystal via
B = (8 pi^2 / 3) RMSF^2 (eq:bfactor), B = crystallographic B-factor. Simulation minus B-factor-derived
RMSF is near zero for good integrator/timestep combinations; few trends correlate with timestep,
implying timestep/integrator choice has little bearing provided no crash. Lobatto is the most
efficient accurate choice.

## Conclusions

CICMD-specific issues addressed: rigorous generalized-coordinate velocity initialization (modal
coords), efficient CM-velocity nulling, and an explanation+solution for the flying-ice-cube effect.
GNEIMO supports CHARMM/AMBER force fields; Nose-Hoover, Berendsen, and rescaling thermostats; GB/SA
and distance-dependent-dielectric solvation. Long-timescale torsional dynamics with large timesteps
suffers the flying-ice-cube effect under Nose-Hoover unless CM kinetic energy is periodically nulled.
RK4, Lobatto, and adaptive CVODE were implemented and tested; both RK4 and Lobatto were stable and
reproduced crystal flexibility for normalized timesteps up to 10 fs (Lobatto) / 4 fs (RK4); within
range, stability metrics were largely timestep-independent. The 2nd-order Lobatto integrator allows
the largest normalized timestep and is the most efficient for GNEIMO.

A dynamic-clustering coarse-graining toolkit was developed; four proteins were folded from extended
structures to molten-globule-like states within 4-5 Å of crystal. Not yet numerically tested:
(1) multiple chains with explicit solvent, (2) the Fixman correction potential to remove the
systematic bias from holonomic constraints in thermodynamic-property calculations (a SOA-based
computational framework for the Fixman compensating potential was derived, ref 39).
