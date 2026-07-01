# Torsion Angle Dynamics: Reduced Variable Conformational Sampling Enhances Crystallographic Structure Refinement

Luke M. Rice and Axel T. Brünger (1994), *Proteins: Structure, Function, and Genetics*.

## Abstract

A reduced-variable conformational sampling strategy for macromolecules based on
molecular dynamics in torsion angle space is evaluated using crystallographic
refinement as a prototypical search problem. Bae and Haug's algorithm for
constrained dynamics (originally developed for robotics) is used. Their
formulation solves the equations of motion *exactly* for arbitrary holonomic
constraints, and hence differs from commonly used approximation algorithms
(e.g. SHAKE). It uses gradients calculated in Cartesian coordinates, and thus
also differs from internal-coordinate formulations. Molecular dynamics can be
carried out at significantly higher temperatures because the high-frequency bond
and angle vibrations are eliminated. The sampling strategy combines high-temperature
torsion angle dynamics with repeated trajectories using different initial
velocities. Best solutions are identified by the free R value (or R value if
experimental phase information is included). For a test system with diffraction
data to 2 Angstrom resolution, slow-cooling protocols fail to converge if the
backbone rms coordinate deviation exceeds 1.25 Angstrom, but torsion angle
refinement can correct deviations up to approximately 1.7 Angstrom.

## Introduction (routing summary)

Conformational sampling underlies free-energy perturbation, structure
prediction, and structure determination by X-ray crystallography and NMR.
Approximate global-minimum searches use simulated annealing, generating many
conformations via molecular dynamics or Monte Carlo. Both methods normally
endow a system of N atoms with the full 3N Cartesian degrees of freedom. A
large simplification follows from recognizing that torsion angles (about one
tenth of the total degrees of freedom) encode most of the conformational
variability. A sampling strategy restricted to torsion angles profits from this
~10-fold reduction in the number of variables.

Molecular dynamics under holonomic constraints has two families of methods:
SHAKE (iterative constraint correction each step; an approximation, impractical
at high temperature), internal-coordinate formulations (exact but force
calculation in internal coordinates is complex), and spatial-operator / recursive
formulations that use small (6x6) matrices to convert Cartesian forces into
arbitrary internal coordinates (formally similar to Kalman filtering/smoothing).
The approach chosen here (Bae and Haug) is equivalent to the spatial-operator
one — small matrices map Cartesian forces into internal coordinates — but is
derived using D'Alembert's Principle instead of matrix-inversion techniques,
making the equations more physically transparent.

## Theory

### Crystallographic refinement

Refinement is the global minimization of the objective function

<!-- eq:1 -->
$$ E = E_{\text{chem}} + w_{\text{X-ray}}\, E_{\text{X-ray}} $$

where $E_{\text{chem}}$ contains empirical information about chemical
interactions, $E_{\text{X-ray}}$ describes the difference between observed and
calculated diffraction data, and $w_{\text{X-ray}}$ is a weight balancing the
forces from each term.

$E_{\text{chem}}$ is a function of all atomic positions describing covalent
(bond lengths, bond angles, torsion angles, chiral centers, planarity) and
nonbonded (van der Waals, hydrogen bonding, electrostatic) interactions:

<!-- eq:2 -->
$$ E_{\text{chem}} = \sum_{\text{bonds}} k_{\text{b}}(r - r_0)^2
+ \sum_{\text{angles}} k_{\theta}(\theta - \theta_0)^2
+ \sum_{\text{dihedrals}} k_{\phi} \cos(n\phi + d)
+ \sum_{\text{chiral,planar}} k_{\omega}(\omega - \omega_0)^2
+ \sum_{\text{atom-pairs}} \left( a\,r^{-12} + b\,r^{-6} + c\,r^{-1} \right) $$

In practice, since crystallographic refinement is not very sensitive to the
accuracy of the empirical energy function, a simpler "geometric" energy function
is often substituted for $E_{\text{chem}}$ in which all nonbonded terms are
replaced by a purely repulsive quartic term.

The most common form of $E_{\text{X-ray}}$ is the crystallographic residual, the
sum over squared differences between observed $|F_{\text{obs}}(\mathbf{h})|$ and
calculated $|F_{\text{calc}}(\mathbf{h})|$ structure factor amplitudes:

<!-- eq:3 -->
$$ E_{\text{X-ray}} = \sum_{\mathbf{h}} \left[ |F_{\text{obs}}(\mathbf{h})| - k\,|F_{\text{calc}}(\mathbf{h})| \right]^2 $$

A phase-restraint penalty term based on the difference between observed and
model-calculated phases can be added:

<!-- eq:4 -->
$$ E_{\text{X-ray}} = \sum_{\mathbf{h}} \left[ |F_{\text{obs}}(\mathbf{h})| - k\,|F_{\text{calc}}(\mathbf{h})| \right]^2
+ w_{\text{p}} \sum_{\mathbf{h}} f\!\left[\phi_{\text{obs}}(\mathbf{h}) - \phi_{\text{calc}}(\mathbf{h})\right] $$

where $w_{\text{p}}$ is the weight of the phase restraint and $f$ is a square-well
function of width equal to the arc cosine of the figure of merit $fom(\mathbf{h})$
for each reflection. An alternative that restrains the real and imaginary parts
of the structure factor (the "vector residual"):

<!-- eq:5 -->
$$ E_{\text{X-ray}} = \sum_{\mathbf{h}} fom(\mathbf{h}) \left\{
\left[A_{\text{obs}}(\mathbf{h}) - k A_{\text{calc}}(\mathbf{h})\right]^2
+ \left[B_{\text{obs}}(\mathbf{h}) - k B_{\text{calc}}(\mathbf{h})\right]^2 \right\} $$

where $A$, $B$ denote the real and imaginary parts of the structure factor.

Simulated annealing requires a target function, a mechanism to generate a
Boltzmann distribution of conformations at temperature $T$, and an annealing
schedule $T_1 > T_2 > \cdots > T_l$. Molecular dynamics integrates Newton's
equations of motion:

<!-- eq:6 -->
$$ m_i \frac{\partial^2 \mathbf{r}_{i}}{\partial t^2} = -\frac{\partial E}{\partial \mathbf{r}_{i}} $$

Here $\mathbf{r}_i$ and $m_i$ are the coordinates and mass of atom $i$ and $E$ is
the potential energy (the target comprising $E_{\text{chem}}$ and
$E_{\text{X-ray}}$). Temperature coupling adds pseudo-friction forces
proportional to the atomic velocities $\mathbf{v}_i$:

<!-- eq:7 -->
$$ \mathbf{f}_{i}^{\text{frict}} = \beta_{i}\left(\frac{T_{0}}{T} - 1\right)\mathbf{v}_{i} $$

where $T_0$ is the bath temperature, $\beta_i$ a force constant, and $T$ the
actual (kinetic) temperature.

### Torsion angle dynamics

Unconstrained MD for $N$ particles solves $\ddot{\mathbf{r}}_i = \mathbf{a}_i(\mathbf{r})
= \mathbf{F}_i(\mathbf{r})/m_i$. For holonomic constraints there are two approaches:
(1) switch to generalized internal coordinates $\mathbf{q}_i$ (gradients hard to
compute but accelerations depend only on positions, so Verlet integration
works); or (2) keep Cartesian coordinates (gradients straightforward and
topology-independent) at the cost that the acceleration becomes a function of
both positions and velocities:

<!-- eq:8 -->
$$ \mathbf{a}(\mathbf{r},\dot{\mathbf{r}}) = \mathbf{M}^{-1}(\mathbf{r})\,\mathbf{Q}(\mathbf{r},\dot{\mathbf{r}}) $$

where $\mathbf{a}$ is the system acceleration vector and $\mathbf{M}$, $\mathbf{Q}$
denote the system inertia matrix and generalized force vector. Because the
acceleration depends on velocity, a fourth-order Runge-Kutta integration scheme
is used (not Verlet).

The authors derive the equations of motion in a form specific to torsion angle
dynamics (the full generality is Bae and Haug; the more complete form is in the
Appendix).

Consider two bodies $i$ and $j$ connected by a bond of fixed length
$|\mathbf{h}_{ij}|$, with the only allowable relative motion a rotation about
$\mathbf{h}_{ij}$. Vectors (Fig. 1): $\mathbf{r}_i$, $\mathbf{r}_j$ locate the
centers of mass of bodies $i$ and $j$ w.r.t. an inertial frame; $\mathbf{s}_{ij}$,
$\mathbf{s}_{ji}$ are the attachment points of each body w.r.t. its own center of
mass; $\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i$; and scalar $q_{ij}$ is the
relative rotation angle about the bond.

The single-rotation constraint relates the angular velocities of the two centers
of mass in the lab frame:

<!-- eq:9 -->
$$ \boldsymbol{\omega}_j = \boldsymbol{\omega}_i + \hat{\mathbf{h}}_{ij}\,\dot{q}_{ij} $$

where $\hat{\mathbf{h}}_{ij} = \mathbf{h}_{ij}/|\mathbf{h}_{ij}|$ is the unit vector
along the bond. The position of body $j$:

<!-- eq:10 -->
$$ \mathbf{r}_{j} = \mathbf{r}_{i} + \mathbf{r}_{ij}
= \mathbf{r}_{i} + \mathbf{s}_{ij} + |\mathbf{h}_{ij}|\,\hat{\mathbf{h}}_{ij} - \mathbf{s}_{ji} $$

Differentiating and rearranging gives the center-of-mass velocity of body $j$ in
terms of that of body $i$:

<!-- eq:11 -->
$$ \dot{\mathbf{r}}_{j} = \dot{\mathbf{r}}_{i}
- \mathbf{r}_{ij} \times \boldsymbol{\omega}_{i}
- (\hat{\mathbf{h}}_{ij} \times \mathbf{s}_{ji})\,\dot{q}_{ij} $$

Vector notation to compact the equations:

<!-- eq:12 -->
$$ \mathbf{Y}_i = \begin{bmatrix} \dot{\mathbf{r}}_i \\ \boldsymbol{\omega}_i \end{bmatrix} $$

<!-- eq:13 -->
$$ \mathbf{q}_{ij} = [q_{ij}] $$

<!-- eq:14 -->
$$ \delta \mathbf{Z}_i = \begin{bmatrix} \delta \mathbf{r}_i \\ \delta \boldsymbol{\pi}_i \end{bmatrix} $$

where $\mathbf{Y}$ is the $6\times 1$ vector of translational and angular
velocities, $\mathbf{q}$ the $1\times 1$ relative angular velocity, and
$\delta \mathbf{Z}$ the $6\times 1$ vector of translational ($\delta \mathbf{r}_i$)
and angular ($\delta \boldsymbol{\pi}_i$) virtual displacements.

The velocity relationship in compact form:

<!-- eq:15 -->
$$ \mathbf{Y}_{j} = \begin{bmatrix} \mathbf{I} & -\tilde{\mathbf{r}}_{ij} \\ \mathbf{O} & \mathbf{I} \end{bmatrix} \mathbf{Y}_{i}
+ \begin{bmatrix} -\hat{\mathbf{h}}_{ij} \times \mathbf{s}_{ji} \\ \hat{\mathbf{h}}_{ij} \end{bmatrix} \dot{\mathbf{q}}_{ij}
= \mathbf{B}_{ij}^{(1)} \mathbf{Y}_{i} + \mathbf{B}_{ij}^{(2)} \dot{\mathbf{q}}_{ij} $$

where $\tilde{\mathbf{r}}_{ij}$ is the skew-symmetric (cross-product) matrix of
$\mathbf{r}_{ij}$. The same notation relates virtual displacements:

<!-- eq:16 -->
$$ \delta \mathbf{Z}_{j} = \mathbf{B}_{ij}^{(1)}\, \delta \mathbf{Z}_{i} + \mathbf{B}_{ij}^{(2)}\, \delta \mathbf{q}_{ij} $$

Differentiating the velocity relationship yields the acceleration relationship:

<!-- eq:17 -->
$$ \dot{\mathbf{Y}}_{j} = \mathbf{B}_{ij}^{(1)} \dot{\mathbf{Y}}_{i} + \mathbf{B}_{ij}^{(2)} \ddot{\mathbf{q}}_{ij} + \mathbf{D}_{ij},
\qquad
\mathbf{D}_{ij} = \dot{\mathbf{B}}_{ij}^{(1)} \mathbf{Y}_{i} + \dot{\mathbf{B}}_{ij}^{(2)} \dot{\mathbf{q}}_{ij} $$

The dynamics is solved using D'Alembert's principle: the virtual work of the
applied forces vanishes over displacements that do not violate the constraints:

<!-- eq:18 -->
$$ \sum_{\text{bodies } k} \delta \mathbf{Z}_k^{T}\left(\mathbf{M}_k \dot{\mathbf{Y}}_k - \mathbf{Q}_k\right) = 0 $$

where $\mathbf{M}_k$ and $\mathbf{Q}_k$ are the inertia matrix and generalized
force vector. One then solves for the acceleration of one body and the relative
acceleration between them to obtain $\mathbf{a}(\mathbf{r},\dot{\mathbf{r}})$, and a
Runge-Kutta step integrates Eq. (8).

### Recursive articulated-body algorithm (prose)

This is a recursive algorithm; the two-body equations extend to many. Atoms are
grouped into rigid bodies, allowing only torsion-angle motions between bodies.
The connectivity defines a tree topology with one arbitrarily chosen body as the
base (root). Each MD step:

1. Start with positions, velocities, and forces for all atoms.
2. Compute center-of-mass positions, velocities, and forces for each rigid group.
3. **Inward sweep (reduction):** starting at the tips of the tree, each chain is
   reduced one body at a time by solving for the relative acceleration between a
   tip and its direct inner body; the tip's inertial properties are then
   aggregated into the inner body, shortening the chain by one link. Continue
   until an expression for the base body's acceleration is obtained.
4. Solve for the base's acceleration (requires inversion of only a $6\times6$ matrix).
5. **Outward sweep (expansion):** the acceleration of each body outboard of the
   base is determined by the base's acceleration plus the relative accelerations,
   propagating outward until the whole tree is covered.
6. A Runge-Kutta step updates positions and velocities; new forces are computed;
   repeat.

The formalism handles several tree-like topologies and closed loops (e.g.
disulfide bonds). In this implementation closed loops are treated approximately
by harmonic distance restraints.

## Methods (test case)

Test refinements: alpha-amylase inhibitor 1HOE-467A, 74 amino acids. Extended
version of X-PLOR. Force field PARHCSDX (Engh and Huber). Diffraction data 70%
complete for $|F_{\text{obs}}| > 2\sigma$. Mean figure of merit of MIR phases =
0.69 for 5-3 Angstrom data (max ~0.82 low-res, min ~0.57 high-res). Resolutions:
5-2 Angstrom ("high") and 5-3 Angstrom ("medium").

Protocol per refinement: 100 steps conjugate-gradient minimization against
$E_{\text{chem}}$ (with $C^{\alpha}$ harmonically restrained) to obtain ideal
bond/angle geometry that stays fixed during torsion angle MD; then 4 ps of
constant-temperature (5000 K or 10,000 K) MD (conventional or torsion angle)
against $E^{*} = E^{*}_{\text{chem}} + w_{\text{X-ray}} E_{\text{X-ray}}$ where
$E^{*}_{\text{chem}}$ uses a purely repulsive nonbonded term. $E_{\text{X-ray}}$
= residual [Eq. 3] for high-res, vector residual [Eq. 5] for medium-res. Timestep
0.002 ps for torsion angle dynamics, 0.0005 ps for conventional. Then 0.1 ps
quench by conventional MD at 300 K, then 120 steps of minimization. Refinements
repeated 10 times from the same initial model with different initial velocities.

## Results and Discussion (routing summary)

Scrambled models generated by increasingly long 600 K MD (without X-ray data)
from a good initial model (backbone/non-H rms 0.4 / 0.8 Angstrom). Three test
models: least/medium/most scrambled, with calculated mean phase differences
$\Delta\varphi_{\text{calc}}$ (5-2 Angstrom) of 79, 82, 84 degrees and backbone
rms 1.25, 1.45, 1.63 Angstrom respectively.

Key findings: slow-cooling simulated annealing corrects mean phase errors up to
~80 degrees (~1.3 Angstrom backbone rms). High-temperature torsion angle
dynamics (10,000 K) plus repeated refinements converges for the most scrambled
model where conventional and 5000 K torsion angle dynamics both fail. 10,000 K is
unattainable for unconstrained MD due to numerical instability from
high-frequency vibrations; constraining bond lengths and angles removes these
vibrations. $R_{\text{free}}$ correlates strongly with backbone rms and is the
practical selector of best models. Torsion angle refinement is ~20% faster than
conventional (Table I). Doubling the high-temperature stage (4 -> 8 ps) further
increases the success rate.

## Conclusions

Torsion-angle-constrained MD is a powerful conformational search tool, giving
significantly increased convergence over conventional refinement at high and
medium resolution. Success derives from combining (1) simplification of
conformational space by constraining bond lengths and angles, (2) high-temperature
sampling, and (3) repeated refinements with different initial conditions. Best
solutions are identified by the free R value (or R value with experimental phase
information). The strategy generalizes to further degree-of-freedom reduction
(constraining secondary-structure elements or tertiary domains), and to other
applications (NMR structure determination, de novo prediction).

## Derivation (not implemented) — Appendix: equations of motion under torsion constraints

The vectors $\mathbf{s}_{ij}$, $\mathbf{s}_{ji}$ are fixed in the center-of-mass
frame of each body; the bond vector $\mathbf{h}_{ij}$ is fixed in the frame of
body $i$. Their time derivatives are $\boldsymbol{\omega}_i \times \mathbf{s}_{ij}$,
$\boldsymbol{\omega}_j \times \mathbf{s}_{ji}$, and $\boldsymbol{\omega}_i \times
\mathbf{h}_{ij}$. These give rise to the compact velocity relation (repeat of
Eq. 15):

<!-- eq:19 -->
$$ \mathbf{Y}_{j} = \begin{bmatrix} \mathbf{I} & -\tilde{\mathbf{r}}_{ij} \\ \mathbf{O} & \mathbf{I} \end{bmatrix} \mathbf{Y}_{i}
+ \begin{bmatrix} -\hat{\mathbf{h}}_{ij} \times \mathbf{s}_{ji} \\ \hat{\mathbf{h}}_{ij} \end{bmatrix} \dot{\mathbf{q}}_{ij}
= \mathbf{B}_{ij}^{(1)} \mathbf{Y}_{i} + \mathbf{B}_{ij}^{(2)} \dot{\mathbf{q}}_{ij} $$

D'Alembert's principle: for a system in equilibrium the net virtual work
vanishes, $\sum_i \mathbf{F}_i \cdot \delta \mathbf{r}_i = 0$; for a system not in
equilibrium, $\sum_i (\mathbf{F}_i - \dot{\mathbf{p}}_i) \cdot \delta \mathbf{r}_i = 0$.
Restricting to displacements $\delta \mathbf{Z}_i$ that do not violate the
constraints yields the constrained equations of motion:

<!-- eq:20 -->
$$ \sum_{\text{bodies } k} \delta \mathbf{Z}_k^{T}\left(\mathbf{M}_k \dot{\mathbf{Y}}_k - \mathbf{Q}_k\right) = 0 $$

The body inertia matrix (spatial mass matrix):

<!-- eq:21 -->
$$ \mathbf{M}_{i} = \begin{bmatrix}
m_{i} & 0 & 0 & 0 & 0 & 0 \\
0 & m_{i} & 0 & 0 & 0 & 0 \\
0 & 0 & m_{i} & 0 & 0 & 0 \\
0 & 0 & 0 & I_{xx} & I_{xy} & I_{xz} \\
0 & 0 & 0 & I_{xy} & I_{yy} & I_{yz} \\
0 & 0 & 0 & I_{xz} & I_{yz} & I_{zz}
\end{bmatrix} $$

The generalized force vector:

<!-- eq:22 -->
$$ \mathbf{Q}_{i} = \begin{bmatrix}
F_{i}^{x} \\ F_{i}^{y} \\ F_{i}^{z} \\
N_{i}^{x} - (\boldsymbol{\omega}_{i} \times \mathbf{I}\boldsymbol{\omega}_{i})_{x} \\
N_{i}^{y} - (\boldsymbol{\omega}_{i} \times \mathbf{I}\boldsymbol{\omega}_{i})_{y} \\
N_{i}^{z} - (\boldsymbol{\omega}_{i} \times \mathbf{I}\boldsymbol{\omega}_{i})_{z}
\end{bmatrix} $$

where $\mathbf{I}$ is the $3\times3$ inertia tensor and $N_i^{x,y,z}$ are the
lab-frame net torque components about the center of mass. The Euler term
$\boldsymbol{\omega}_i \times \mathbf{I}\boldsymbol{\omega}_i$ appears because the
net torque is $\mathbf{N}_i = \mathbf{I}\dot{\boldsymbol{\omega}}_i +
\boldsymbol{\omega}_i \times \mathbf{I}\boldsymbol{\omega}_i$ while only
$\mathbf{I}\dot{\boldsymbol{\omega}}_i$ appears in $\mathbf{M}_i\dot{\mathbf{Y}}_i$.

Expanding the D'Alembert sum for two bodies:

<!-- eq:23 -->
$$ \delta \mathbf{Z}_{i}^{T}\left(\mathbf{M}_{i}\dot{\mathbf{Y}}_{i} - \mathbf{Q}_{i}\right)
+ \delta \mathbf{Z}_{j}^{T}\left(\mathbf{M}_{j}\dot{\mathbf{Y}}_{j} - \mathbf{Q}_{j}\right) = 0 $$

Substituting the constraint relations (Eqs. 16, 15, 17, repeated as Eq. 24):

<!-- eq:24 -->
$$ \delta \mathbf{Z}_{j} = \mathbf{B}_{ij}^{(1)}\, \delta \mathbf{Z}_{i} + \mathbf{B}_{ij}^{(2)}\, \delta \mathbf{q}_{ij},
\quad
\mathbf{Y}_{j} = \mathbf{B}_{ij}^{(1)}\, \mathbf{Y}_{i} + \mathbf{B}_{ij}^{(2)}\, \dot{\mathbf{q}}_{ij},
\quad
\dot{\mathbf{Y}}_{j} = \mathbf{B}_{ij}^{(1)}\, \dot{\mathbf{Y}}_{i} + \mathbf{B}_{ij}^{(2)}\, \ddot{\mathbf{q}}_{ij} + \mathbf{D}_{ij} $$

yields (Eq. 25):

<!-- eq:25 -->
$$ 0 = \delta \mathbf{Z}_{i}^{T}\left(\mathbf{M}_{i}\dot{\mathbf{Y}}_{i} - \mathbf{Q}_{i}\right)
+ \left(\delta \mathbf{q}_{ij}^{T}\mathbf{B}_{ij}^{(2)T} + \delta \mathbf{Z}_{i}^{T}\mathbf{B}_{ij}^{(1)T}\right)
\left\{\mathbf{M}_{j}\left[\mathbf{B}_{ij}^{(1)}\dot{\mathbf{Y}}_{i} + \mathbf{B}_{ij}^{(2)}\ddot{\mathbf{q}}_{ij} + \mathbf{D}_{ij}\right] - \mathbf{Q}_{j}\right\} $$

Factoring terms in $\delta \mathbf{Z}_i$ and $\delta \mathbf{q}_{ij}$ (Eq. 26):

<!-- eq:26 -->
$$ \delta \mathbf{Z}_{i}^{T}\left\{
\mathbf{M}_{i}\dot{\mathbf{Y}}_{i}
+ \mathbf{B}_{ij}^{(1)T}\mathbf{M}_{j}\mathbf{B}_{ij}^{(1)}\dot{\mathbf{Y}}_{i}
+ \mathbf{B}_{ij}^{(1)T}\mathbf{M}_{j}\mathbf{B}_{ij}^{(2)}\ddot{\mathbf{q}}_{ij}
+ \mathbf{B}_{ij}^{(1)T}\mathbf{M}_{j}\mathbf{D}_{ij}
- \mathbf{Q}_{i}
- \mathbf{B}_{ij}^{(1)T}\mathbf{Q}_{j}\right\}
+ \delta \mathbf{q}_{ij}^{T}\left\{
\mathbf{B}_{ij}^{(2)T}\mathbf{M}_{j}\mathbf{B}_{ij}^{(1)}\dot{\mathbf{Y}}_{i}
+ \mathbf{B}_{ij}^{(2)T}\mathbf{M}_{j}\mathbf{B}_{ij}^{(2)}\ddot{\mathbf{q}}_{ij}
+ \mathbf{B}_{ij}^{(2)T}\mathbf{M}_{j}\mathbf{D}_{ij}
- \mathbf{B}_{ij}^{(2)T}\mathbf{Q}_{j}\right\} = 0 $$

Because $\delta \mathbf{Z}_i$ and $\delta \mathbf{q}_{ij}$ are linearly
independent by construction, both coefficients must vanish. Setting the
$\delta \mathbf{q}_{ij}$ coefficient to zero and solving for the relative
angular acceleration (Eq. 27):

<!-- eq:27 -->
$$ \ddot{\mathbf{q}}_{ij} = -\left(\mathbf{B}_{ij}^{(2)T}\mathbf{M}_{j}\mathbf{B}_{ij}^{(2)}\right)^{-1}
\left\{\mathbf{B}_{ij}^{(2)T}\mathbf{M}_{j}\mathbf{B}_{ij}^{(1)}\dot{\mathbf{Y}}_{i}
+ \mathbf{B}_{ij}^{(2)T}\left[\mathbf{M}_{j}\mathbf{D}_{ij} - \mathbf{Q}_{j}\right]\right\} $$

This expression for $\ddot{\mathbf{q}}_{ij}$ is substituted into the
$\delta \mathbf{Z}_i$ coefficient, giving an expression for $\dot{\mathbf{Y}}_i$
(the acceleration of body $i$), which in turn determines $\ddot{\mathbf{q}}_{ij}$
and hence the system accelerations. For a branch-free chain of $n$ bodies the
same reduction proceeds inward from the tip (each $\delta \dot{\mathbf{q}}$ is
linearly independent of all other virtual displacements), shortening the chain by
one body at a time until the first body's acceleration is obtained, then
propagating relative accelerations back outward. Branched chains and closed loops
are handled analogously (Bae and Haug).
