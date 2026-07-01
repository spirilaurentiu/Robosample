# Robosample: A rigid-body molecular simulation program based on robot mechanics

Laurentiu Spiridon, Teodor Asvadur Şulea, David D.L. Minh, Andrei-Jose Petrescu.
BBA - General Subjects (2020). DOI: 10.1016/j.bbagen.2020.129616.

## Abstract

Compared with all-atom molecular dynamics (MD), constrained MD methods allow for
larger time steps, potentially reducing computational cost, motivating continued
interest in improving constrained MD algorithms to increase configuration-space
sampling. Robosample is a software package that implements high-performance
constrained dynamics algorithms, originally developed for robotics, and applies
them to biomolecular simulations. As in the gMolmodel package (Spiridon and Minh,
2017), Robosample uses Constrained Dynamics Hamiltonian Monte Carlo (CDHMC) as a
Gibbs sampling move - a Monte Carlo move where a subset of coordinates is allowed
to change. In addition to Cartesian and torsional dynamics moves, Robosample
implements spherical and cylindrical joints that can be distributed along the
molecule by the user. In alanine dipeptide simulations, the free energy surface is
recovered by mixing fully flexible with torsional, cylindrical, or spherical
dynamics moves. Ramachandran dynamics, where only the two key torsions are mobile,
accelerate the slowest transition by an order of magnitude. Simulations of a
complex glycan cover significantly larger regions of configuration space when
mixed with constrained dynamics.

## 1. Introduction

Structural biology has provided abundant information about free-energy minima of
biomolecules, but processes involving flow over an energy landscape (folding,
large conformational changes) may evade a simple structural description. Many
biomolecules contain innately flexible moieties (intrinsically disordered regions,
glycans) that can only be described statistically. Simulation accuracy is often
limited by conformational sampling.

One approach to enhanced sampling is the application of holonomic constraints.
Standard MD integrators require small time steps to model high-frequency motion;
constrained dynamics eliminates the highest-frequency degrees of freedom, allowing
longer time steps. Early methods include SHAKE and RATTLE (constrain bonds with
hydrogen); torsional dynamics (constrain bond lengths and angles) has been applied
to crystallographic and NMR structure refinement, loop modelling, and rapid
folding.

In a constrained MD simulation, molecules may be represented as multibody chains
where atoms or *groups of atoms* are rigid bodies and joints map onto chemical
bonds. Such representations reflect energy as a function of internal molecular
motion more naturally than Cartesian coordinates, and allow coarse graining that
preserves a fine-grained energy function: a generalized-coordinate formulation may
drop DOFs during dynamics while preserving the energy function (e.g. a phenyl group
as one rigid body, yet all ring atoms contribute to the resulting torque).

A limitation of constrained dynamics is ergodicity. With constraints, samples are
drawn from a conditional rather than joint probability distribution. Vaidehi and
Jain suggested freezing/releasing DOFs during MD; Kandel et al. relaxed torsional
dynamics by allowing flexible bond angles. Spiridon and Minh used constrained MD to
generate candidate configurations for Markov chain Monte Carlo - an example of
Gibbs sampling, in which a Markov chain is propagated by sampling from a conditional
distribution. As long as every DOF can be modified by at least one Monte Carlo
move, ergodicity can in principle be achieved.

Robosample implements Hamiltonian Monte Carlo based on constrained MD and fully
flexible MD using algorithms borrowed from robotics. It was built from the Molmodel
API and its generalization gMolmodel (used for Gibbs sampling in AlGDock).
Robosample refactored and optimized gMolmodel, works as a standalone application
and API, and implements many joint types (Ball, Cylinder, Universal, Slider, etc.).
It emphasizes spherical and cylindrical joints that capture torsion-angle and
torsion-bond couplings. Robosample also calculates the gradient of the mass matrix
determinant, allowing MD trajectories to sample from the conditional probability of
the Boltzmann distribution.

## 1.1. The constrained generalized-coordinates molecular simulation basis

Simulations with `m` constraints use two formulations: the Lagrange multiplier
formulation in Cartesian coordinates (up to 3N DOFs plus m constraint equations),
and the reduced-coordinate formulation based on `N_f = 3N - m` generalized
coordinates obtained from a coordinate transformation. The limiting case of both is
m=0, N_f=3N (0 multiplier equations, maximum reduced coordinates). Starting from
this limiting case with a transformation `{φ} ∈ R^{3N×3N}`, `{φ}` is split into m
coordinates `{φ_v}` that are constrained and the free set `{φ_f}` (the reduced
coordinate set).

In the fully flexible limiting case, dynamics follow Hamilton's equations (see
`equations.md` eq:1), with a Hamiltonian combining kinetic and potential energy
(eq:2), where the mass matrix tensor is `M_tot = J^T M J` with J the Jacobian of the
coordinate transformation and M the diagonal Cartesian mass matrix (eq:2-def). The
Boltzmann probability of a microstate in the canonical ensemble is eq:3.

Most thermodynamic quantities of interest depend on configurations, not momenta.
The marginal configuration probability follows from a Gaussian integral over
momenta (eq:4), producing a `|M_tot|^{-1/2}` prefactor.

A system in reduced coordinates evolves according to Hamiltonian eq:5, with marginal
distribution eq:6. There is a **mismatch between eq:4 and eq:6**: the determinants
differ. To conduct constrained simulations whose unnormalized density matches that
of unconstrained dynamics, one uses the Fixman correcting potential (eq:7). There
are at least two ways to use it so that the marginal probability of flexible
coordinates matches in constrained and fully flexible regimes: (1) include the
Fixman potential in the Metropolis-Hastings acceptance criterion (HMC); (2) compute
its derivative, the Fixman torque, and apply it during MD. With the Fixman
correction the reduced-coordinate Hamiltonian is eq:8.

Generalized-coordinate dynamics is more computationally complex than fully flexible
dynamics: solving for acceleration requires inverting the mass matrix tensor, and
the Fixman torque is required for MD trajectories to sample the Boltzmann
distribution. Direct inversion is O(N^3), but O(N) algorithms exploit the tensor's
special structure. One is the **Spatial Operator Algebra (SOA)** developed by
Rodriguez and Kreutz (1988), a multibody dynamics formulation for internal
coordinates. SOA relies on the equivalence between Kalman filtering /
Bryson-Frazier smoothing equations in signal processing and the dynamics of
kinematic chains: spatial forces per link are input, joint torques are output. SOA
uses linear operators acting on velocities, accelerations, and forces, associated
with recursive algorithms that pass through the chain without matrix multiplication
or inversion. SOA also yields the mass matrix tensor determinant, its gradient,
logarithm, and square root - quantities needed for the statistical mechanics of
constrained systems.

## 2. Methods

### 2.1. Robosample architecture

Robosample builds on the Simbody and Molmodel APIs. Simbody was designed for
general, accurate mechanical-engineering calculations with biomedical emphasis; it
implements SOA and offers integrators for the Hamilton equations. To sample from
the Boltzmann distribution, the authors added libraries for the determinant, square
root, and gradient of the mass matrix tensor. Molmodel represents a chemical object
as a tree and maps Simbody serial multibody subgraphs (dependent on the defined
rigid bodies) onto the chemical tree. Robosample calls Simbody for trajectory
integration while energies and forces are evaluated in Molmodel and can be computed
on GPU via the OpenMM API.

Note (Fig. 1): the main speed bottleneck is where non-bonded energy terms must be
evaluated in Cartesian coordinates. OpenMM provides GPU capability.
Note (Fig. 2): each "world" has different rigid-body definitions and joint types,
mapped onto the same chemical graph; rigid bodies are grouped and joints connect
them.

#### 2.1.1. Gibbs sampling implementation

Robosample uses constrained MD to sample from the Boltzmann distribution following
two conditions:

1. Each sample is drawn from the conditional probability of the Boltzmann
   distribution. Velocities are initialized per a generalized equipartition theorem
   for generalized coordinates. Trajectories are propagated deterministically, and
   the final configuration is accepted/rejected by the Metropolis-Hastings
   criterion.
2. Every DOF can be modified by some Monte Carlo move. This is satisfied by
   alternating different blocks of coordinates (e.g. torsions) that are flexible in
   a given Gibbs move. Each block specifies a "world" (a defined set of rigid bodies
   and joints), implemented as a different Simbody multibody chain mapped onto the
   same Molmodel molecular graph. At the start of each Gibbs move, Cartesian
   coordinates update the generalized coordinates in the relevant world. Including a
   fully flexible world ensures every DOF is accessible.

#### 2.1.2. Parameters setup

Simulating multiple worlds causes a combinatorial increase in parameters: number of
worlds, rigid-body specifications, joint types, number of MD steps per HMC proposal,
and MD time step. General principles: define rigid bodies by secondary structure or
domains; sample informative variables (torsions) more often than bond
lengths/angles; include at least one fully flexible world per cycle for ergodicity.

HMC proposal parameters are optimized via autocorrelation time and acceptance rate.
For MD trajectories, the time step (ε) and trajectory length (L, number of steps)
must be specified; integration time is T = L·ε. Too short T causes a random walk and
long autocorrelation; increasing T reduces autocorrelation; too long T reduces the
effective sample size per L. The acceptance rate is the major consideration for the
time step: large ε reduces steps needed for a given T but increases integrator error
and reduces acceptance; small ε minimizes error but needs more steps. The acceptance
rate optimally balancing these is **0.651**.

Workflow for optimizing parameters:
- a. run short single-world trials of different T with ε = 1 fs;
- b. pick the maximum samples-per-transition over T;
- c. given the chosen T, increase ε in mixed trials until acceptance ~0.651.

A GUI assists with setup: one window for general parameters (temperature, number of
cycles), one tabular window for per-world parameters (MD steps, integration time
step).

#### 2.1.3. Spherical and cylindrical joints

Torsional dynamics is attractive due to the low frequency of torsional motions and
their ability to travel large distances in conformational space, but torsions are
strongly coupled to other DOFs. GneimoSim added angle-bend mobilities; Robosample
introduces a **spherical joint** (three rotational DOFs, expressed as Euler angles
or quaternions) and a **cylindrical joint** (two DOFs: a translation and a rotation
about the translational axis). Both were implemented in Molmodel based on the Simbody
Ball and Cylinder mobilized-body types.

A child rigid body attaches to its parent through a mobilizer (joint), described by
four reference frames: one for each rigid body and two defining the joint (one fixed
on the parent, one mobile on the child). In modelling molecules, the two mobilizer
frames are naturally placed at atom centers oriented along bonds, but in general
robotics they can be placed anywhere within the bodies. Placing joint frames
randomly (not at atom centers) results in poorer transition rates (Suppl. Table S3);
Robosample places mobilizer frames along bonds by default, letting the user specify
the fixed and mobile atoms.

### 2.2. Experimental setups

Two systems tested: alanine dipeptide, and a medium-sized model of the first glycan
of hepatitis C virus protein E2 (E2N1) attached to the E2-derived [−12 + 6] peptide
around its Asn 17 site.

#### 2.2.1. Alanine dipeptide

Input files generated with AmberTools16 using ff14SB. Two dynamic models: (a) 7
mobile joints (all bonds except 2 terminal) - "all-bond dynamics" ('all-'); (b) 2
joints (N-Cα and Cα-C bonds) - "Ramachandran dynamics" ('Rama-'). Each model set for
three move types: torsional (TD), torsion-bond cylindrical (Cyl), and torsion-angle
spherical (Ball) - 6 worlds/regimens mixed with fully flexible worlds for
ergodicity. Separate fully flexible runs performed for reference.

HMC parameters tuned by trials with increasing trajectory lengths and a short fixed
step of 2 fs (1 fs for the flexible world). Autocorrelation time of the end-to-end
distance selects the trajectory length maximizing effective sample size per
picosecond; the time step is then increased until acceptance near 0.651.

Four simulation types compared: (1) fully flexible; (2) mixed with torsional
(mixed-TD); (3) mixed with cylindrical (mixed-CYL); (4) mixed with spherical
(mixed-BALL). Three instances: (a) all-bond, unoptimized; (b) all-bond, optimized;
(c) optimized Ramachandran - 11 experiments total. Each cycle = 1 fully flexible
world + 10 CDHMC worlds. Each experiment repeated 3 times.

#### 2.2.2. E2N1 model

The E2N1 glycan attached to peptide [−12 to +6], peptide kept fixed in the
antibody-recognized conformation. Polypeptide coordinates from an HCV-E2 homology
model (Modeller v9.21); glycan from Carbohydrate Builder v2.1 (glycam.org). Input
prepared with AmberTools16 using ff14SB (polypeptide) and GLYCAM06 (polysaccharide).
Structure minimized with OpenMM v7.4 until 1 kJ/mol.

CDHMC simulations at 300 K with four worlds: (1) fully flexible; (2) rigid
monosaccharide rings; (3) mixed spherical/cylindrical for the polysaccharide; (4)
torsional at the Asn anchor. HMC trajectory lengths × time steps: 108 × 1 fs,
25 × 8 fs, 158 × 0.5 fs, 11 × 20 fs respectively. In worlds 2 and 3, monosaccharide
hydroxyl, carbonyl, secondary amino methyl, and methylene groups were kept rigid;
spherical joints between rigid polysaccharide groups, cylindrical where only two
atoms. Each simulation repeated three times.

Four simulation types: two MD (1a fully flexible MD, 1b rigid body MD / RB-MD) and
two HMC (2a fully flexible, 2b multiple rigid body worlds / CDHMC). Rigid bodies
covered monosaccharides or their lateral groups.

Trajectories analyzed with VMD 1.9.4, NumPy 1.11.0, SciPy 0.17; visualized with VMD;
charts with Matplotlib 1.5.1 and Excel. Conformational-space exploration assessed by
all-vs-all RMSD analysis and cluster count analysis (Ward's method hierarchical
agglomerative clustering, flattened at 2 Å). Runtime measured in seconds. Hardware:
Intel i9-9900K (3.60 GHz), 2× GeForce RTX 2080Ti, 64 GB RAM, Ubuntu 18.04. E2N1
flexible MD (747 DOF) ~0.022 s/step; RB-MD (105 DOF) ~0.008 s/step; CDHMC world 2
(74 DOF) ~0.007 s/step; world 3 identical to RB-MD; world 4 (4 DOF) ~0.001 s/step.

## 3. Results

### 3.1. Mixed simulations reconstruct the free energy surface of alanine dipeptide

Alanine dipeptide (N-acetylalanine-N-methylamide) is small but has a highly
frustrated potential energy surface. Projecting configuration space onto φ and ψ
retains most PES maxima/minima; the potential of mean force over φ,ψ is widely used
to assess sampling methods.

To validate reproduction of prior results, simulations used the same parameters as
Spiridon and Minh 2017; fully flexible and mixed-TD results closely resembled prior
results. Optimized all-bond trajectory lengths and time steps: fully flexible 108 ×
1.87 fs (202.0 fs), torsional 25 × 17.1 fs (427.5 fs), ball 158 × 2.1 fs (331.8 fs),
cylindrical 47 × 6.1 fs (286.7 fs). Ramachandran: torsional 11 × 44.73 fs (492.0 fs),
ball 51 × 8.957 fs (456.8 fs), cylindrical 49 × 12.386 fs (606.9 fs). Ball and
cylinder dynamics use larger time steps than flexible but much smaller than
torsional; larger rigid bodies of Ramachandran dynamics allow larger steps.

All simulations recovered the free energy surface. CDHMC based on the new
cylindrical and spherical joints accurately recovers the surface irrespective of
body size. Simulations at 625 K have basins in the same positions as prior results
at 625 K and at 800 K.

Optimization increased sampling efficiency of rare events. For all types except
all-bond Ball, optimized parameters give faster transitions to/from the αL basin
(crossing φ=0). The all-bond Ball case improved only slightly - the Ball joint is
less subject to diffusive motion.

Mean first passage times between basins (Table 4): fully flexible simulations have
the slowest transitions between C5, PPII, C7eq, and αL. In all-bond dynamics, Ball
and Cylinder joints perform significantly better than flexible Cartesian, torsional
best. All-bond dynamics are outperformed by Ramachandran dynamics. Compared to fully
flexible, MFPT across the φ=0 barrier decreases by an order of magnitude; MFPT
between other basins on the same side decreases by about half. Joint efficiency
increases in order: Torsion, Cylinder, Ball.

### 3.2. Using multiple worlds improves simulation efficiency

Glycans are flexible polysaccharides attached to proteins, rarely present in solved
structures, and thus of interest for conformational analysis. Multiple cyclic
structures make them a good target for rigid-body simulation; Cartesian simulations
are hindered by high-frequency motion within monosaccharides (mostly chair
conformation).

Both MD and HMC simulations were performed using fully flexible vs mixed movements
on the glycan-peptide system.

### 3.3. Robosample increases simulation efficiency

In equal wall-clock time, the mean and fluctuation of RMSD are significantly larger
in mixed-world versus fully flexible HMC simulations. Cluster count analysis
corroborates: HMC simulations involving rigid bodies access more clusters in the
same run time. The contrast is only evident with HMC (multiple rigid-body worlds),
not MD (single rigid-body regime). Multiple-world HMC routinely reaches RMSD up to
30 Å from the initial structure, with distinct structures across independent runs;
fully flexible MD/HMC stay around 15 Å or less.

## 4. Discussion and future directions

Robosample is a good environment for constrained molecular simulations achieving
ergodicity via Gibbs sampling. Beyond torsional and angle/torsion mobilities, it
introduces robotic mobility combinations by exploiting Simbody's mechanical joints,
allowing examination of arbitrary DOF couplings (e.g. weak bond/torsion couplings in
cylindrical joints).

A simple HMC optimization plus different joint types gave enhanced transitions
to/from the αL basin. With all-bond dynamics, all joint types outperformed flexible
dynamics, torsional providing the greatest acceleration. Ramachandran dynamics (only
φ, ψ mobile) was significantly more efficient than fully-flexible or all-bond, with
Ball giving the shortest MFPT for the rarest transition.

The glycan system confirms alternating constraints explore larger configuration-space
volumes than fully flexible simulations (cluster count, RMSD matrix). Robosample may
also aid structure prediction: switching DOF types via joints (cylindrical moves for
side chains, spherical for loops) can sample large configuration-space volumes in
search of the deepest minimum.

Limitation: rigid-body simulation cannot accurately simulate time-dependent
quantities (correlated motion, diffusion coefficients, binding kinetics); coarsening
loses accuracy due to entropic effects in lost DOFs (e.g. rigid monosaccharides
prevent observing chair-envelope-boat transitions).

Future work: optimizing world combinations (flexible HMC followed by multiple
constrained HMC moves); spherical and cylindrical moves may be complementary;
combination with replica exchange; additional Simbody joints (Slider).

## Key references (harvested)
- SHAKE: Ryckaert, Ciccotti, Berendsen 1977. RATTLE: Andersen 1983.
- Fixman potential: Fixman 1974, 1978. Fixman torque: Jain 1997; Jain et al. 2013.
- SOA: Rodriguez, Kreutz, Jain 1989; Rodriguez 1987; Rodriguez, Jain, Kreutz-Delgado 1991.
- Simbody: Sherman, Seth, Delp 2011. Molmodel: Flores et al. 2011. OpenMM 7: Eastman et al. 2017.
- CDHMC/gMolmodel: Spiridon & Minh 2017. AlGDock: Minh 2019.
- HMC: Duane et al. 1987. Metropolis-Hastings: Hastings 1970. Gibbs sampling: Geman & Geman 1984.
- IC-MD equipartition: Jain, Park, Vaidehi 2012. GneimoSim / potential-distortion fix: Kandel et al. 2016. IC-MD foundation: Vaidehi & Jain 2015.
- Optimal HMC tuning / acceptance 0.651: Beskos et al. 2013; Betancourt 2015, 2016.
