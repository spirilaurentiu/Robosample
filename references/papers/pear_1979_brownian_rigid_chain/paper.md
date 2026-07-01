# Brownian dynamics study of a polymer chain of linked rigid bodies

M. R. Pear and J. H. Weiner, Department of Physics and Division of Engineering, Brown University.
*J. Chem. Phys.* **71**, 212 (1979). doi:10.1063/1.438119

## Abstract

A Brownian dynamics model for the backbone chain of a macromolecule is developed
as a system of linked rigid bodies so that constraints on valence angles and bond
lengths are satisfied exactly. For comparison, a corresponding flexible model is
developed in which bond lengths and valence angles are held nearly constant by strong
harmonic potentials. Equilibrium properties and barrier crossing rates are examined
theoretically and by computer simulation of both models, with differences arising due
to the presence of constraints in the rigid case. A compensating potential based on
the metric determinant of unconstrained coordinates in the rigid model is found to
eliminate the effect of constraints. Barrier crossing rates in the transition state
approximation are studied when a force fixed in space is applied to the end atoms of
the three-bond chain. An exact transition state rate formula developed for this case
predicts curved Arrhenius plots of barrier crossing rates; this result is confirmed by
computer simulation.

## I. Introduction

Long-chain molecules transition readily between conformations because rotation about
backbone bonds is opposed only by low energy barriers. The kinetic aspects (rates of
transitions vs. chain characteristics and barrier heights) are modeled here with the
Langevin equation. In the Langevin equation a systematic damping force and a random
fluctuating force are added to the equations of motion, so that for a particle of mass
$m$ in a potential $U$:

<!-- eq:1.1 -->
$$ m\ddot{\mathbf{x}} = -\nabla U - \eta \, \dot{\mathbf{x}} + \mathbf{L}(t) $$

where $\eta$ is the viscosity and $\mathbf{L}(t)$ is the fluctuating Langevin force,
satisfying the statistical relations of Eq. (1.2). "Brownian dynamics" denotes
simulations based on the Langevin equation, as opposed to molecular dynamics where the
heat-bath degrees of freedom are simulated explicitly.

This paper presents a computer program for simulating a long-chain molecule with fixed
bond lengths and valence angles. Distinctive features: (a) inertia effects are included
(applicable to both low- and high-viscosity regimes); (b) geometric constraints of
fixed bond lengths and valence angles are satisfied *exactly* by treating the chain as
a collection of linked rigid bodies, using the Wittenburg formalism for systems of
rigid bodies; (c) a simplified treatment of the Langevin forces is used (requires a
smaller time step than more elaborate schemes for comparable accuracy).

The general procedure is specialized to a three-bond chain with 90-degree valence
angles to examine: (a) equilibrium statistics of rigid vs. flexible models; (b) barrier
crossing rates; (c) effect of applied stress on barrier crossing rates.

**Rigid vs. flexible equilibrium statistics.** The rigid model constrains bond lengths
and valence angles as geometric constraints (only dihedral rotation remains). The
flexible model maintains bond lengths and valence angles nearly constant via harmonic
potentials with large spring constants. The classical equilibrium statistics of the
two models differ: a *metric determinant* appears in the configuration-space
distribution for rigid models but is absent for flexible models, no matter how large
the spring constants. Go and Scheraga conclude that the flexible model is more accurate
for equilibrium properties (the metric determinant should not appear), and that beyond
reaching that conclusion the flexible model may be treated as rigid for computing
equilibrium properties. For the freely rotating three-bond chain with 90-degree valence
angle, the dihedral distribution is uniform for the flexible model and nonuniform
(proportional to $\sqrt{g(\phi_3)}$) for the rigid model. Fixman's compensating
potential $U(\phi) = k_B T \ln\sqrt{g(\phi)}$ makes the rigid model reproduce the
flexible model's uniform distribution.

## II. The Model

### A. Linked rigid body formulation

We apply the Wittenburg formalism, a general systematic numerical method for solving
the motion of a system of connected rigid bodies, starting from Newton's law and the
law of moment of momentum. Definition requires: (1) connectivity, (2) description of
each rigid body, (3) kinematics of the connection points (hinges).

**Geometry (Fig. 1).** Carbon atoms $C_{i-2}, C_{i-1}, C_i, C_{i+1}$. Bond vector
$\mathbf{l}_i$ is segment $C_{i-1}C_i$; valence angle $\theta_i$ is the angle between
$\mathbf{l}_i$ and $\mathbf{l}_{i-1}$; dihedral angle $\phi_i$ measures rotation of
$\mathbf{l}_i$ about the axis along bond $\mathbf{l}_{i-1}$. The angle $\phi_i = 0$
corresponds to a *trans* conformation.

The rigid chain with $N$ bonds constrains valence angles $\theta_i$ ($i=2,\dots,N$) and
bond lengths $|\mathbf{l}_i|$ ($i=1,\dots,N$) to constant values $\theta$ and $l$. Only
the $N-2$ dihedral angles $\phi_i$ ($i=3,\dots,N$) vary.

**Mechanical model (Fig. 2).** Each link consists of one backbone atom with two
massless connecting rods of length $l/2$ forming angle $\theta$. Hinges limit rotation
to the axis along the bond connecting two backbone atoms. A torque internal to each
hinge produces a rotational potential governing that hinge's dihedral angle. In addition
to internal degrees of freedom there are six extra: three for translation of the chain
center of mass (eliminated from the formulation) and three Bryant angles
$\phi_0, \phi_1, \phi_2$ (like Euler angles) for the orientation of the first body. The
state is described by $N+1$ angles $\phi_0,\dots,\phi_N$: $\phi_0,\phi_1,\phi_2$ orient
the first body relative to a fixed frame; $\phi_3,\dots,\phi_N$ are the internal dihedral
angles.

**Concepts (specialized from Wittenburg for a linear chain):**

(1) *Connectivity.* A tree structure; only one path between any two atoms. The
connectivity matrix $\underline{T}$ has $T_{ij} = -1$ if $C_j$ is on the path between
$C_i$ and the body containing $C_{N-1}, C_N$. See Eq. (2.1). The hinge between $C_j$ and
$C_{j+1}$ is the upper hinge of $C_j$ and the lower hinge of $C_{j+1}$.

(2) *Augmented bodies.* The augmented body $B_j$ for body $C_j$ collapses atoms
$C_0,\dots,C_{j-1}$ to a point mass at the lower hinge and $C_{j+1},\dots,C_N$ to the
upper hinge. $B_j$ has a new center of mass (barycenter); its moment of inertia tensor
about the barycenter is $K_j^*$. Three body vectors $\mathbf{b}_{j0}, \mathbf{b}_{jj},
\mathbf{b}_{jN}$ give the positions of the lower hinge, the center of mass of body $C_j$,
and the upper hinge, relative to the barycenter of $B_j$. For end groups $\mathbf{b}_{10}$
and $\mathbf{b}_{N-1,N}$ are zero.

(3) *Kinematics of hinges.* The relative angular velocity of body $C_i$ w.r.t. $C_{i-1}$
is Eq. (2.2); extended to the first body via Eq. (2.3); the absolute angular velocity is
Eq. (2.4).

(4) *Internal hinge torques.* Torques internal to hinges incorporate the rotational
potential $V_\phi$; the torque on body $C_{i-1}$ through the lower hinge along
$\mathbf{p}_i$ has magnitude Eq. (2.5), with $\tau_i = 0$ for $i=0,1,2$.

The Wittenburg method yields the matrix equation Eq. (2.6), $\underline{A}\ddot{\phi} =
\underline{B}$, with $\underline{A}$ and $\underline{B}$ given by Eqs. (2.7)-(2.8). For a
chain with $N$ bonds, $\underline{\phi}$ is a column of $N+1$ elements, $\underline{A}$ is
a symmetric positive-definite $(N+1)$-dimensional matrix, and $\underline{B}$ a column of
$N+1$ elements. The matrices $\underline{\mathbf{p}}$ (Eq. 2.9) and $\underline{K}$
(Eq. 2.10) are defined below; $\underline{\mathbf{f}}$ (Eq. 2.11), $\underline{\mathbf{M}}'$
(Eq. 2.12), and the Langevin/damping terms $\underline{\mathbf{L}}$ (Eqs. 2.13-2.15) and
$\underline{\mathbf{M}}_\eta$ (Eqs. 2.16-2.20) follow.

All tensorial quantities are first defined in a body-fixed coordinate system, then
transformed to a common reference frame (successive transformations similar to Flory's
chain-statistics transformations; Appendix A) to compute $\underline{A}$ and
$\underline{B}$.

### B. Langevin force simulation

The stochastic Langevin force $\mathbf{L}(t)$ is described via its impulse
$\mathbf{B}(\delta t)$ (Eq. 2.21), whose Cartesian component satisfies the Gaussian
probability distribution Eq. (2.22). The solution method has two steps per time step:
(1) a fourth-order Runge-Kutta integration from $t$ to $t+\delta t$ *without* the
Langevin force; (2) increment the velocities with a random impulse having the statistics
of Eq. (2.22). For the flexible (Cartesian) case, the velocity update is Eqs. (2.23)-(2.24).
For the rigid model the Cartesian Langevin impulses (Eq. 2.24) are transformed to
generalized-coordinate velocity increments (Eq. 2.25) by substituting the impulses in
place of the Langevin forces in Eqs. (2.13)-(2.15) to form an impulse matrix
$\hat{\underline{L}}$.

This Langevin treatment is elementary compared to the iterative scheme of Fixman or the
series method of Helfand (which contains this approximation as one term). Accuracy is
controlled by varying the time step $\delta t$.

### C. The three-bond chain

For $N=3$ no interior links appear in the rigid model; the chain has a single internal
degree of freedom ($\phi_3$). The end atoms $C_0, C_3$ may have greater mass than the
interior atoms $C_1, C_2$, with parameter $\alpha$ = ratio of mass of $C_0$ to $C_1$ (and
of $C_3$ to $C_2$), crudely modeling the rest of a chain acting on an interior three-bond
section. In the flexible model, reduction to three bonds gives one dihedral $\phi_3$ plus
five additional internal degrees of freedom (bond lengths $l_1,l_2,l_3$ and valence angles
$\theta_2,\theta_3$).

## III. Equilibrium statistics in flexible and rigid models

### A. Theory

As the harmonic forces in the flexible model become stronger, it approaches the rigid
chain as a limiting case, but the equilibrium statistics differ: taking thermal averages
and letting spring constants go to infinity do not commute. Comparison uses the three-bond
chain with 90-degree valence angles and no internal rotational potential.

For the flexible model, the internal-coordinate equilibrium distribution is Eq. (3.1);
integrating over the other coordinates gives the *uniform* dihedral distribution
Eq. (3.2), even in the infinite-spring-constant limit. For the rigid model, the metric
determinant produces the nonuniform distribution Eq. (3.3) with normalization Eq. (3.4).
Fixman's theorem computes the metric determinant $g$ from the matrix $H_{ij}$ (Eq. 3.5),
which for the three-bond 90-degree chain yields the closed form Eq. (3.6) with
coefficients $a,b,c$. Adding the compensating Fixman potential Eq. (3.7) eliminates the
metric determinant and restores the uniform distribution Eq. (3.8).

### B. Computer simulation results

Comparing flexible (F), rigid without compensation (R), and rigid with Fixman potential
(RF), with $\alpha=10$ to emphasize the metric-determinant effect: $\rho_F$ is uniform;
$\rho_R$ agrees with the theoretical Eq. (3.3) using the determinant Eq. (3.6); with the
Fixman potential $\rho_{RF}$ is again uniform.

## IV. Barrier crossing rates

### A. Theory

A single-barrier potential on $[-\pi,\pi]$ is used, Eq. (4.1), with $E_b = k/4$ so that
the piecewise-quadratic potential is continuous with continuous derivatives. Parameter
$\phi_b$ sets the barrier-peak position. Two rates are considered: (1) the transition-state
rate (all crossings of the barrier peak counted); (2) the effective rate (crossing complete
only after the entire barrier region is traversed). The transition-state rate approximates
the effective rate well at low viscosity but increasingly overestimates it as viscosity
$\eta$ grows.

The transition-state rate for the rigid model is computed exactly from Eq. (4.2), with the
distribution Eq. (4.3), kinetic energy Eq. (4.4), covariant metric tensor Eq. (4.5), and
$g(\phi_3) = |G_{ij}|$. Using the Gaussian integral formula Eq. (4.6), the normalization
constant becomes Eqs. (4.7)-(4.8); with the velocity integral Eq. (4.9), Eq. (4.2) yields
the transition-state rate Eq. (4.10). $G^{33}$ is the indicated component of the
contravariant metric tensor $G^{ij}$.

Eq. (4.10) predicts a splitting of the crossing rate depending on the barrier position
$\phi_b$. Using the compensating Fixman potential, $g(\phi_b)$ drops out and the rate
becomes Eq. (4.11), independent of barrier position (for $\phi_b=0$ and $\phi_b=\pi$,
where $G^{33}$ has the same value) and higher than the uncompensated rates.

### B. Computer simulation results

The splitting predicted by Eq. (4.10) is verified for $\alpha=10$. For $\delta t=0.1/\omega$
(with $\omega=1.33\times10^{13}\ \mathrm{sec}^{-1}$), agreement is poor at low temperature
due to the approximate thermal simulation; halving to $\delta t=0.05/\omega$ improves
agreement. The Fixman potential coalesces the $\phi_b=0$ and $\phi_b=\pi$ curves and gives
higher rates. Effective rates behave qualitatively the same. The flexible model shows no
metric determinant, so its effective rate is independent of barrier-peak position.

## V. Applied stress

### A. Theory

Forces $\sigma$ and $-\sigma$ are applied to end atoms $C_0$ and $C_3$. The
nonuniform-$\phi_3$ effect is removed with the compensating Fixman potential. The total
potential (excluding the Fixman potential) is Eq. (5.1); with the direction of $\sigma$ as
a fixed orientation vector, Eq. (5.2), where $R(\phi) = l(3+2\cos\phi)^{1/2}$ is the
end-to-end distance and $\Theta$ is the angle between $\mathbf{R}$ and $\sigma$.

Because the form of Eq. (5.2) leaves the minimum unchanged, the saddle point at
$\phi_3=\phi_b,\ \Theta=0$, and a *planar* critical hypersurface, the transition-state rate
can be computed exactly (a two-dimensional rate theory). The orientation angles must be
explicitly considered (Eq. 5.3), with the canonical distribution Eq. (5.4). Using the
coordinates $\Theta$ (angle between $\mathbf{R}$ and $\sigma$), azimuthal $\Phi$, and
$\Psi$ (rotation about $\mathbf{R}$), the rate is Eq. (5.5), reduced to Eq. (5.6) with
normalization Eqs. (5.7)-(5.8). Performing the $\Theta$ integration gives the final form
Eq. (5.9) with $I_c(\sigma)$ from Eq. (5.10), evaluated numerically with
$R(\phi)=l(3+2\cos\phi)^{1/2}$. In the small-$\sigma$ limit, Eq. (5.9) reduces to the
zero-stress result Eq. (4.11).

Eq. (5.9) predicts a separation of the rate curves for $\phi_b=\pi$ and $\phi_b=0$, now due
to applied stress (not the metric determinant, which the Fixman potential removed). The
Arrhenius plots are curved, with curvature depending on $\phi_b$, and splitting increasing
at low temperature.

### B. Computer simulation results

Including $\sigma$ requires adding the moments Eq. (5.11) to $\underline{\mathbf{M}}'$ for
the three-bond chain. Simulations with $\alpha=1$, $\sigma l/E_b=0.5$ support the predicted
curvature and splitting (Fig. 10); with $\sigma l/E_b=0.25$ the splitting is reduced,
especially at high temperature (Fig. 11).

## VI. Conclusions

A Brownian-dynamics model of a macromolecule backbone with exactly satisfied bond-length
and valence-angle constraints was developed using the Wittenburg linked-rigid-body
formalism. This formalism generalizes readily to side-group inertia, branched chains, and
chains with more internal flexibility, making it more valuable here than usual Lagrangian
techniques. A corresponding flexible model (bond lengths / valence angles held nearly
constant by strong springs) was also developed. Both were specialized to the three-bond
90-degree chain.

For no internal rotational potential, the flexible model gives a uniform dihedral
distribution and the rigid model a nonuniform one proportional to $\sqrt{g}$; adding the
Fixman potential $U = k_B T \ln\sqrt{g(\phi)}$ restores uniformity. An exact
transition-state rate formula showed a dependence of rate on barrier-peak position for the
uncompensated rigid model, removed by the Fixman potential (and absent in the flexible
model). Effective rates behaved qualitatively similarly. Because the analytic
transition-state formulas are exact, they serve as a check on numerical accuracy;
sufficiently small $\delta t$ was required.

Applied-stress study produced an exact two-dimensional transition-state formula predicting
inherent Arrhenius-plot curvature, traced to the temperature-dependent average
misorientation $\langle\cos\Theta\rangle$ between $\sigma$ and the end-to-end vector
$\mathbf{R}$. Curvature directions are opposite for $\phi_b=0$ and $\phi_b=\pi$. Simulation
and theory agree well.

## Appendix A: Details of linked rigid body formalism

### 1. Coordinate transformations

A body-fixed coordinate system is defined for each of the $N-1$ rigid bodies (following
Flory). For body $C_i$ ($i=1,\dots,N-1$): the $x$ axis is along the bond from $C_{i-1}$ to
$C_i$; the $y$ axis is perpendicular to $x$ and in the plane of $C_{i-1}, C_i, C_{i+1}$; the
$z$ axis completes a right-handed system. A vector $\mathbf{V}$ with components in the frame
of body $C_i$ (written $\mathbf{V}^{(i)}$) is expressed in the frame of $C_{i-1}$ by
Eq. (A1) using the matrix $T(\phi)$ of Eq. (A2). The transformation from the first body
frame to the reference frame uses the matrix $a^{01}$ of Eq. (A3) with $c_i=\cos\phi_i,\
s_i=\sin\phi_i$.

### 2. Description of individual bodies

For interior bodies $B_2,\dots,B_{N-2}$ ($i=2,\dots,N-2$) the augmented body vectors are
Eq. (A4), with total mass $M$, $\beta_i = N+\alpha-1-i$, and $\alpha$ = ratio of the mass of
$C_0$ to $C_1$. After transforming to the reference frame, the augmented body tensor $K_i^*$
is Eq. (A5). The first body's vectors are Eq. (A6) and its augmented inertia tensor is
Eq. (A7) with $a,b,c$ from Eq. (A8). The other terminating body is obtained by transforming
with $\mathbf{R}$ of Eq. (A9): $\mathbf{b}_{N-1,N}=0$, $\mathbf{b}_{NN}=\mathbf{R}\mathbf{b}_{11}$,
$\mathbf{b}_{N0}=\mathbf{R}\mathbf{b}_{12}$, $K_N^* = \mathbf{R}K_N^*\mathbf{R}^T$.

### 3. Hinge vectors

For $i=2,\dots,N-1$ the vectors $\mathbf{p}_{i+1}$ are simply unit vectors along the $x$ axis
of the frame fixed to body $C_i$ (the internal-hinge rotation axis is along the bond). For
the frame of $C_1$, the vectors $\mathbf{p}_0, \mathbf{p}_1, \mathbf{p}_2$ are Eq. (A10);
their partial derivatives (needed in Eq. 2.11) then follow.

## Appendix B: Details of flexible chain formalism

If $\mathbf{x}_i=(x_1^{(i)},x_2^{(i)},x_3^{(i)})$ is the Cartesian position of backbone
atom $C_i$ in a chain of $N+1$ atoms, the generalized coordinates $(l_i,\theta_i,\phi_i)$
are Eqs. (B1a)-(B1c). The $3N+3$ equations of motion are Eq. (B2), with gradient operator,
bond potential $V_l$ (Eq. B3), valence-angle potential $V_\theta$ (Eq. B4), and rotational
potential $V_\phi$ (Eq. B5). Substituting Eqs. (B1) and differentiating w.r.t. Cartesian
coordinates gives the equations of motion Eq. (B6), where $\mathbf{F}_{\alpha_k}^j$ is the
force on atom $C_j$ due to the presence of $\mathbf{x}_j$ in generalized coordinate
$\alpha_k$ (zero if $\alpha_k$ is undefined for that $k$). The bond-stretching forces are
Eq. (B7), the valence-angle forces Eq. (B8) with $K_{\theta_j}$, $\mathbf{e}_j$ (Eq. B9),
$d_j^k$ (Eq. B10), and the rotational forces Eq. (B13) using $a_j^k$ (Eq. B11) and
$^iA_b^j$ (Eq. B12).

## Units

Reduced units: the mass $m$ of a carbon atom, the length $l$ of a C-C bond, and a reference
temperature $T_{\mathrm{ref}} = 600\ \mathrm{K}$ are all unity. A time step $\Delta t = 1$
in this system corresponds to $\Delta t = 7.51\times10^{-14}\ \mathrm{sec}$ in cgs units.
