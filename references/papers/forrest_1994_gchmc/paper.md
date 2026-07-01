# Generalized coordinate hybrid Monte Carlo

B. M. Forrest and U. W. Suter, Institut für Polymere, ETH Zürich.
Molecular Physics, 1994, 82:2, 393-410. DOI: 10.1080/00268979400100304.

## Abstract

A novel hybrid Monte Carlo (HMC) algorithm for the off-lattice simulation of
dense, atomistically detailed polymer systems. Bond lengths and bond angles are
constrained to their equilibrium values and a generalized-coordinate description
using the torsional degrees of freedom is employed. Simple, decoupled,
Cartesian-like equations of motion are obtained by introducing **fictitious
angular momenta** conjugate to the torsional degrees of freedom. A unique time
scale, and hence a mapping onto real time, is obtained by allocating
time-independent, **effective moments of inertia** to each angular momentum
variable.

## 1. Introduction

Efficient thermalization requires generating statistically independent
configurations at constant temperature. For dense polymer systems, sample
configurations are hard to alter because of chain connectivity and lack of free
volume, causing strong correlations between successive samples.

MD is limited primarily by high-frequency bond-length and bond-angle stretching
motions, forcing a time-step of at most a few femtoseconds. Constraining bond
lengths/angles to equilibrium values removes these motions but introduces
constraint-equation convergence problems. An alternative is a
generalized-coordinate description in the torsional degrees of freedom: a more
natural coordinate system for fixed bond lengths/angles that gives at least a
threefold (usually order-of-magnitude) reduction in degrees of freedom, but at
the cost of much more complicated equations of motion requiring matrix inversion.

Conventional MC alters only a few degrees of freedom per step to keep acceptance
rates appreciable, giving long autocorrelation times at high density. This work
presents an HMC algorithm that is a global-update method (like MD) but is exact
(like MC), providing the configuration-space partition function and free of
numerical instabilities and systematic discretization errors. The approach uses
*fixed* bond lengths and bond angles within a *generalized-coordinate* framework
while avoiding the usual computational complications.

## 2. The hybrid Monte Carlo algorithm

HMC was originally proposed for lattice gauge theories and later shown suitable
for condensed-matter simulation. Global changes are made between MC acceptance
decisions by running an MD integration for a specified number of steps. Momenta
conjugate to the coordinate variables are introduced; Hamilton's equations are
discretized and integrated, propelling the system through configuration space.
Since momenta feel the forces, coordinates tend to evolve toward decreasing
energy (higher acceptance), but random initial momenta also allow uphill moves.

The momenta are **refreshed** (redrawn from a Maxwell-Boltzmann distribution)
after **every** MC acceptance decision, whether or not the move was accepted -
otherwise detailed balance cannot be satisfied. The momenta need not be the
actual momenta conjugate to positions: in general one may introduce *fictitious*
momenta that merely propagate the coordinates. This more general approach is used
here.

HMC retains global updates through MD but without compromising the numerical
stability and exactness of MC. Unlike MD or Langevin, HMC is both numerically
stable and exact regardless of the integration time-step size: it reproduces
unbiased canonical averages, whereas MD/Langevin produce step-size-dependent
observables. The discretization error of the underlying MD is properly accounted
for in the MC acceptance criterion and affects only the acceptance rate. The
one-step version of HMC is equivalent to exact Langevin integration for
equilibrium observables.

## 3. The polybead model

The 'polybead' model of polyethylene: each polymer is a skeletal chain of `Nb`
spherical 'monomers', each a CH2 group (with two terminal CH3 groups) - also
called 'polyethylene in the unified-atom approximation'. Highest-frequency motions
are eliminated by fixing bond lengths and bond angles at equilibrium values
(`lb = 1.53 Å`, bond angle `θ = 112°`).

Each backbone of `Nb` beads is described by `Nb - 3` torsional angles plus three
Euler angles (orientation of the plane of the first three monomers relative to
the lab frame) plus three Cartesian coordinates for the first monomer. The number
of coordinates is `Nb + 3` per chain versus `3 Nb` for a full Cartesian
description. The three Euler angles are termed `φ1, φ2, φ3`; `φ4` is the first
torsional angle (conformation of the first three skeletal bonds), and `φ_Nb` is
the last torsional angle.

The total potential energy of the model (`Nc` chains, each `Nb` beads):

<!-- eq:1 -->
$$\mathscr{E}(\{\boldsymbol{r}_{1}^{(c)}, \{\phi_{k}^{(c)}\}\}) = \sum_{c=1}^{N_{c}} \sum_{k=4}^{N_{b}} V_{\phi}(\phi_{k}^{(c)}) + \sum_{c=1}^{N_{c}} \sum_{i=5}^{N_{b}} \sum_{j=1}^{i-4} V_{LJ}(|\boldsymbol{r}_{i}^{(c)} - \boldsymbol{r}_{j}^{(c)}|) + \sum_{c=2}^{N_{c}} \sum_{i=1}^{N_{b}} \sum_{c'=1}^{c-1} \sum_{j=1}^{N_{b}} V_{LJ}(|\boldsymbol{r}_{i}^{(c)} - \boldsymbol{r}_{j}^{(c')}|)$$

Here `r_1^(c)` denotes the coordinates of the origin of chain `c` and `φ_k^(c)`
is the `k`-th Euler angle (`k = 1, 2, 3`) or `(k-3)`-th torsional angle
(`4 ≤ k ≤ Nb`) of chain `c`. For `i ≥ 4` the position `r_i^(c)` of each bead is a
function of `r_1^(c)` and the angles `φ_1^(c)` through `φ_i^(c)`. The second and
third beads depend only on `r_1^(c)` and the three Euler angles.

The first ('torsional') term reflects the tendency of dihedral angles to adopt
*trans* or *gauche* rotational states, in the Ryckaert-Bellemans form:

<!-- eq:2 -->
$$V_{\phi}(\phi) = C \sum_{n=0}^{5} a_n \cos^n(\phi)$$

with `C = 9.0 kJ/mol`, `a0 = 1`, `a1 = 1.31`, `a2 = -1.414`, `a3 = -0.3297`,
`a4 = 2.828`, `a5 = -3.3943`.

The second and third terms represent intramolecular and intermolecular
'non-bonded' internal energy, modelled by a 12-6 Lennard-Jones potential:

<!-- eq:3 -->
$$V_{\rm LJ}(r_{ij}) = 4\varepsilon \left[ \left( \frac{\sigma}{r_{ij}} \right)^{12} - \left( \frac{\sigma}{r_{ij}} \right)^{6} \right]$$

with `ε = 410 J/mol` and `σ = 3.94 Å` for all beads (no distinction between
internal and terminal beads). The intramolecular component acts only between
beads separated by at least four bonds, since the torsional potential (2)
accounts for the interaction between four successive beads.

In the fully flexible polybead model, with `3 Nc Nb` Cartesian coordinates `{x}`
and `3 Nc Nb` conjugate momenta `{p}`, the canonical average of an observable
`O({x})` at temperature `T` is

<!-- eq:4 -->
$$\langle O \rangle \equiv \langle O(\{x\}) \rangle = \frac{1}{C} \int dp \exp\left[-\frac{p^2}{2m_b k_B T}\right] \int dx \, O(\{x\}) \exp\left[-\frac{\mathscr{E}(\{x\})}{k_B T}\right]$$

where `mb` is the bead mass, `kB` the Boltzmann constant, `C` a normalization
constant. The Gaussian integral over the Cartesian momenta factors out, leaving a
configurational average:

<!-- eq:5 -->
$$\langle O(\{x\}) \rangle = \frac{1}{Z} \int dx \ O(\{x\}) \exp \left[ -\frac{\mathscr{E}(\{x\})}{k_{\rm B} T} \right]$$

with configurational partition function `Z = ∫ dx exp[-E({x})/kB T]`.
Transforming from the full unconstrained Cartesian description `{x}` to the
generalized-coordinate description `{q}` of the corresponding rigid model
(`{q}` is shorthand for `{r_1^(c), φ_k^(c); 1 ≤ c ≤ Nc, 1 ≤ k ≤ Nb}`) gives

<!-- eq:6 -->
$$\langle O \rangle = \frac{1}{Z_q} \int dq \ D(\{q\})O(\{q\}) \exp \left[ -\frac{\mathscr{E}(\{q\})}{k_B T} \right]$$

with `Zq = ∫ dq D({q}) exp[-E({q})/kB T]`. The factor `D({q})` has two sources:
(1) the Jacobian of the Cartesian-to-generalized-coordinate transformation
(still within the unconstrained model), and (2) the corrective **Fixman term**,
which accounts for the 'freezing' of the degrees of freedom associated with
constraining bond lengths and bond angles in the rigid model.

An MC algorithm explores `{q}` space to generate configurations distributed
according to the Boltzmann probability `PB ∝ D({q}) exp[-E({q})/kB T]` in eq (6).
A move is suggested from the old configuration `{q0}` to a new candidate `{qn}`
with proposal probability `PS[{q0} → {qn}]` and accepted with probability

<!-- eq:7 -->
$$P_{\text{acc}} = \min \left\{ 1, \frac{P_{\text{S}}[\{q_n\} \to \{q_0\}] D(\{q_n\}) \exp(-\mathscr{E}(\{q_n\})/k_{\text{B}}T)}{P_{\text{S}}[\{q_0\} \to \{q_n\}] D(\{q_0\}) \exp(-\mathscr{E}(\{q_0\})/k_{\text{B}}T)} \right\}$$

For any proposal scheme, choice (7) satisfies detailed balance and ensures
convergence to the Boltzmann distribution of eq (6).

The aim is to sample *only* the configurational integral (6), not the full
phase-space integral (4). Following Duane et al., moves through `{q}` space are
suggested by introducing artificial dynamics using the actual internal energy
`E({q})` and a 'kinetic energy' involving fictitious Cartesian-like momenta.
Even though the fictitious momenta do not appear in eq (6), their change must be
included (in addition to the change in `D({q}) exp[-E({q})/kB T]`) because they
appear in the ratio `PS[{qn} → {q0}] / PS[{q0} → {qn}]` in eq (7). The internal
energy function used in the equations of motion can itself be fictitious, but the
actual internal energy must be used in the acceptance factor.

Following common practice in off-lattice polymer simulation, the corrective
**Fixman pseudopotential** is ignored. In the acceptance function (7) only the
contribution from the Cartesian-to-Euler-angle transformation is explicitly
included (present in all rigid-body simulations).

All results are for a melt of 20 polybead chains, each 24 CH2 units, in a box of
length 25 Å with periodic boundary conditions, at `T = 480 K`. This gives density
≈ 0.68 g/cm3, roughly polyethylene at this temperature at 1 atm.

## 4. The equations of motion

A tilde denotes *fictitious* variables. A Hamiltonian `H̃ = E + K̃` is constructed
by augmenting the potential energy (1) with the fictitious kinetic energy:

<!-- eq:8 -->
$$\tilde{K}(\{p_1^{(c)}\}, \{\tilde{\pi}_k^{(c)}\}) = \sum_{c=1}^{N_c} \frac{(p_1^{(c)})^2}{2m_c} + \sum_{c=1}^{N_c} \sum_{k=1}^{N_b} \frac{(\tilde{\pi}_k^{(c)})^2}{2\tilde{I}_k^{(c)}}$$

The first term uses the usual Cartesian momenta `p_1^(c)` conjugate to the
chain-origin coordinates (`mc` = mass of chain `c`), while the fictitious momenta
`π̃_k` are conjugate to the angular variables `φ_k`. The constants `Ĩ_k^(c)` can
be chosen freely and are identified as **effective moments of inertia**.

Had the actual kinetic energy in generalized coordinates been used, the momenta
would be angular momenta and the kinetic energy would contain cross-terms coupled
by the inertia tensor. The simple Cartesian-like form (8) gives **decoupled**
Hamiltonian equations of motion:

<!-- eq:9 -->
$$\frac{\mathrm{d}\tilde{\pi}_{k}^{(c)}}{\mathrm{d}\tilde{t}} = -\frac{\partial \mathscr{E}}{\partial \phi_{k}^{(c)}}, \qquad \frac{\mathrm{d}\phi_{k}^{(c)}}{\mathrm{d}\tilde{t}} = \frac{\tilde{\pi}_{k}^{(c)}}{\tilde{I}_{k}^{(c)}}, \qquad \frac{\mathrm{d}\boldsymbol{p}_{1}^{(c)}}{\mathrm{d}\tilde{t}} = -\frac{\partial \mathscr{E}}{\partial \boldsymbol{r}_{1}^{(c)}}, \qquad \frac{\mathrm{d}\boldsymbol{r}_{1}^{(c)}}{\mathrm{d}\tilde{t}} = \frac{\boldsymbol{p}_{1}^{(c)}}{m_{c}}$$

From the second equation, `φ̇_k = π̃_k / Ĩ_k` is the rate of change of angle `φ_k`
with respect to (fictitious/computer) time `t̃`. The first equation gives the rate
of change of the angular velocity, `dφ̇_k/dt̃ = -(∂E/∂φ_k)/Ĩ_k`. Hence the `Ĩ_k`
are effective moments of inertia and `∂E/∂φ_k` is an effective torque on `φ_k`.

Because equations (9) are Cartesian-like, they are integrated with the **leap-frog**
scheme, which is time-reversible and area-preserving in phase space (ensuring
detailed balance).

First, for each chain, momenta `p_1^(c)` and `{π̃_k^(c)}` are generated from a
Maxwell-Boltzmann distribution at temperature `T`:

<!-- eq:10 -->
$$\Pr\left(\left\{\boldsymbol{p}_{1}, \tilde{\boldsymbol{\pi}}_{k}\right\}\right) \propto \prod_{c=1}^{N_{c}} \left[\exp\left(-\frac{\beta(\boldsymbol{p}_{1}^{(c)})^{2}}{2m_{c}}\right) \prod_{k=1}^{N_{b}} \exp\left(-\frac{\beta(\tilde{\boldsymbol{\pi}}_{k}^{(c)})^{2}}{2\tilde{I}_{k}}\right)\right]$$

where `β = 1/kB T`. Then eq (9) is integrated for `N_MD` time-steps of size
`δt̃_MD` using leap-frog, effecting a proposed move `q(t̃) → q(t̃ + Δt̃)` with
`Δt̃ = N_MD δt̃_MD`. The move is accepted with probability

<!-- eq:11 -->
$$P_{\text{acc}} = \Pr\left(\boldsymbol{q}(\tilde{t}) \to \boldsymbol{q}(\tilde{t} + \Delta \tilde{t})\right) = \min\left\{1, \left[\prod_{c=1}^{N_c} \frac{\sin\left(\phi_2^{(c)}(\tilde{t} + \Delta \tilde{t})\right)}{\sin(\phi_2^{(c)}(\tilde{t}))}\right] \exp\left(-\beta \Delta \tilde{\mathcal{H}}\right)\right\}$$

where `ΔH̃ = E(t̃ + Δt̃) - E(t̃) + K̃(t̃ + Δt̃) - K̃(t̃)` is the discretization error
after time integration. The product over `sin(φ2)` (the second Euler angle of each
chain) must be included to reproduce configuration-space ensemble averages: it
arises from the factor `D({qn})/D({q0})` in eq (7). The change in the *total*
Hamiltonian `H̃` appears because the ratio
`PS[{qn} → {q0}] / PS[{q0} → {qn}]` in eq (7) involves `exp(-β ΔH̃)`. If the move
is rejected, `q(t̃ + Δt̃) = q(t̃)` is reset before continuing.

One MC step (MCS) = momentum generation (10) + leap-frog MD of `N_MD` steps +
acceptance decision (11). Regardless of accept/reject, the momenta are refreshed
via (10) and the procedure repeats, generating a Markov chain
`{q(t̃), q(t̃ + Δt̃), q(t̃ + 2Δt̃), ...}`.

### Derivation (not implemented): contrast with actual generalized-coordinate MD

If the computational angular velocities `φ̇_k` are chosen to be the actual angular
velocities `ω_k`, the MD trajectory is that of generalized-coordinate molecular
dynamics with angular kinetic energy

<!-- eq:12 -->
$$K_{\rm ang}(\{\phi_k^{(c)}\}, \{\omega_k^{(c)}\}) = \frac{1}{2} \sum_{i,j} \sum_{b=1}^{N_b} m_b \left(\frac{\partial \mathbf{r}_b}{\partial \phi_i}\right) \cdot \left(\frac{\partial \mathbf{r}_b}{\partial \phi_j}\right) \omega_i \omega_j \equiv \frac{1}{2} \sum_{i,j} I_{ij}(\{\phi\}) \omega_i \omega_j$$

where `I_ij` is the inertia tensor and `mb` the mass of bead `b`. The angular
equations of motion are then coupled:

<!-- eq:13 -->
$$\sum_{j} I_{ij} \frac{\mathrm{d}\omega_{j}}{\mathrm{d}t} + \sum_{j,k} \frac{\partial I_{ik}}{\partial \phi_{j}} \omega_{k} \frac{\mathrm{d}\phi_{j}}{\mathrm{d}t} = -\frac{\partial \mathscr{E}}{\partial \phi_{i}} + \frac{1}{2} \sum_{j,k} \frac{\partial I_{jk}}{\partial \phi_{i}} \omega_{j} \omega_{k}$$

This requires matrix inversion for `ω̇` and `ω` as functions of time; `I_ij` and
`∂I_ij/∂φ` are time-dependent. Thus (13) is more computationally demanding than
(9) and does not admit a straightforward time-reversible, area-preserving
discretization such as leap-frog. This contrast motivates the fictitious-momentum
construction (8)-(9).

The most expensive part, computing the generalized forces (torques), can be
reduced from `O((Nb Nc)^3)` to `O((Nb Nc)^2)` operations using an iterative scheme
that computes derivatives recursively from the outermost torsional angle inward to
the Euler angles (a generalization of the Abe et al. method to many chains with
periodic boundary conditions).

## 5. Extracting a unique time scale

For any choice of time-independent constants `Ĩ_k^(c)`, the algorithm is
canonical. To reconcile the fictitious algorithmic time scale with real time, an
**instantaneous moment of inertia** `I_k` is associated with angular variable
`φ_k`, imagined by freezing the chain conformation, holding all other angles
fixed, and considering only rotations in `φ_k`:

<!-- eq:14 -->
$$\mathscr{I}_k \equiv \sum_{b \in M_k} m_b [e_k \times (r_b - P_k)]^2$$

Here `M_k` is the set of beads moved by angle `φ_k`, `mb` the mass of bead `b`;
the rotation axis for `φ_k` has unit direction `e_k` and passes through `P_k`.
This is conformation-dependent and time-varying, but in equilibrium fluctuates
around a time-independent average (borne out by figure 1).

For the first three (Euler) angles, each rotates about an axial vector through the
chain origin `r1`, so `e_k` is that axis and `P_k = r1`. For torsional angles 4-24
(`φ4` innermost, `φ24` outermost), `e_k = (r_{k-1} - r_{k-2}) / ||r_{k-1} - r_{k-2}||`
and `P_k = r_{k-1}`. The outermost angle `φ24` has a constant moment of inertia
`mb lb^2 sin^2 θ` since it rotates only the last bead.

The mean instantaneous moments of inertia `⟨I_k⟩` measure the average difficulty
of rotating each angle, and are natural choices for the time-independent constants
`Ĩ_k`. Being equilibrium values, they can be estimated from any canonical
simulation. In the absence of a-priori information one may set `Ĩ_k = constant`
for all `k` (still canonical).

### Derivation (not implemented): effective time-step estimate

Each `φ_k` has an associated mean moment of inertia `⟨I_k⟩` roughly constant in
equilibrium, so by equipartition there is energy `≈ ½ kB T` associated with its
motion: `½ ⟨I_k⟩ ⟨φ̇_k^2⟩ ≈ ½ kB T`, where `φ̇_k` is the mean rate of change of
`φ_k` with respect to real time. If each MC step spans real time interval `Δt`,
then `⟨φ̇_k^2⟩ = ⟨(Δφ_k/Δt)^2⟩ = ⟨(Δφ_k)^2⟩/(Δt)^2`, where `Δφ_k` is the change in
angle `φ_k` during one MC step (monitorable). Solving gives the effective
time-step:

<!-- eq:15 -->
$$\Delta t \equiv N_{\rm MD} \delta t_{\rm MD} \approx \left( \frac{\langle \mathscr{I}_k \rangle \langle (\Delta \phi_k)^2 \rangle}{k_{\rm B} T} \right)^{1/2}$$

Using `Ĩ_k = I_24` for all `k` (with `mb = 14 g/mol`), the angular variables move
on different time scales: the outer torsional angles are penalized (assigned too
large an inertia) and move slower; only `φ24` gets an appropriate value, matching
the three translational chain-origin coordinates. Using `Ĩ_k = ⟨I_k⟩`, all angular
variables move on approximately the same effective time scale, mapping onto real
time: each MCS is an MD trajectory of roughly 150 fs, i.e. each MD step ≈ 15 fs
(mean 14.9 fs over angles; 14.7 fs from chain-origin coordinates).

The choice `Ĩ_k = ⟨I_k⟩` also enhances sampling. The end-to-end vector
autocorrelation function per chain:

<!-- eq:16 -->
$$f_{\text{EEV}}(t) = \left\langle \frac{(\mathbf{r}_{N_b}(t) - \mathbf{r}_1(t)) \cdot (\mathbf{r}_{N_b}(0) - \mathbf{r}_1(0))}{\|\mathbf{r}_{N_b}(t) - \mathbf{r}_1(t)\| \|\mathbf{r}_{N_b}(0) - \mathbf{r}_1(0)\|} \right\rangle$$

The first- and second-degree bond orientational autocorrelation functions:

<!-- eq:17 -->
$$f_{\text{BCF; 1}}(t) = \langle \mathbf{u}_i(t) \cdot \mathbf{u}_i(0) \rangle, \qquad f_{\text{BCF; 2}}(t) = \frac{3}{2} \langle (\mathbf{u}_i(t) \cdot \mathbf{u}_i(0))^2 \rangle - \frac{1}{2}$$

where `u_i` is the unit vector along the `i`-th bond, averaged over all bonds. The
assignment `Ĩ_k = ⟨I_k⟩` gives faster decorrelation and more efficient sampling.

Because of the product over second Euler angles in acceptance (11), the
mean-acceptance identity `⟨exp(-ΔH̃/kB T)⟩ = 1` becomes

<!-- eq:18 -->
$$\left\langle \left[ \prod_{c=1}^{N_c} \frac{\sin\left(\phi_2^{(c)}(\tilde{t} + \Delta \tilde{t})\right)}{\sin\left(\phi_2^{(c)}(\tilde{t})\right)} \right] \exp\left(-\frac{\Delta \tilde{\mathcal{H}}}{k_B T}\right) \right\rangle = 1$$

The original identity is recovered by defining `ΔH̃* = ΔH̃ - kB T ln(Πc)`, where
`Πc` is the product over the second Euler angles, whence `⟨exp(-ΔH̃*/kB T)⟩ = 1`.
Similarly the Gupta et al. relation

<!-- eq:19 -->
$$\langle P_{\rm acc} \rangle \approx \operatorname{erfc} \left( \frac{1}{2} \langle \Delta \widetilde{\mathcal{H}} / k_{\rm B} T \rangle^{1/2} \right)$$

becomes `⟨Pacc⟩ ≈ erfc(½ ⟨ΔH̃*/kB T⟩^{1/2})` in the present case.

## 6. Optimizing the algorithm

Sampling efficiency is optimized by varying the number `N_MD` of integration steps
between acceptance decisions and the time-step `δt̃_MD`. Small `N_MD`, `δt̃_MD` give
small discretization error and high acceptance, but short trajectories keep
consecutive configurations highly correlated (large autocorrelation times). Large
`N_MD`, `δt̃_MD` give large discretization error, small acceptance, and also large
autocorrelation times. An optimal combination minimizes autocorrelation time -
measured not in MC steps but in **MD steps**, since the cost of each MD step is
independent of `N_MD` and `δt̃_MD`.

Efficiency is measured two ways: the mean time `τ_s` for centres of mass to
diffuse a distance equal to the average radius of gyration (`⟨s^2⟩ ≈ 40 Å^2`), and
the mean time `τ_0.75` for the end-to-end vector autocorrelation (16) to decay to
75% of its initial value.

For fixed `δt̃_MD = 1.5×10^-3`, varying `N_MD` shows an optimal trajectory length
`N_MD ≈ 75` for both measures (with a flat minimum). The optimal time-step for
`N_MD = 10` gives ≈ 65% acceptance; optimal values for `N_MD ≥ 100` give more
efficient sampling at lower acceptance. In all cases optimal performance has
`0.3 < ⟨Pacc⟩ < 0.7`. A self-regulating time-step could maintain `⟨Pacc⟩` in this
range.

For a Lennard-Jones fluid the optimal trajectory length was `N_MD ≈ 10`, much
shorter than here. Longer trajectories help the polymer system because when
velocities are refreshed after each MC step, the interdependence of torsional
velocities of beads in the same chain is lost; this correlation is rebuilt during
MD integration, so longer trajectories allow it to develop. Such inter-particle
velocity correlation matters much less for Lennard-Jones particles.

## 7. Conclusion

A novel HMC algorithm for off-lattice simulation of dense, atomistically detailed
polymer systems, using fixed bond lengths and bond angles in a generalized
coordinate framework. Coupled equations of motion are circumvented by auxiliary
angular momenta conjugate to the torsional degrees of freedom, via a fictitious
kinetic energy that is a simple Cartesian-like quadratic sum with constant
coefficients interpreted as effective moments of inertia. Instantaneous moments of
inertia associated with each angular degree of freedom have time-independent
average values in equilibrium; choosing these means for the effective moments of
inertia recovers a unique time scale, mapping computational time to real time.

## Appendix: A brief restatement of the hybrid Monte Carlo method

Given any potential energy `E({q})`, where `{q}` is an arbitrary coordinate system
(need not be Cartesian), introduce a fictitious kinetic energy diagonal in
fictitious momenta `{π̃}` and carry out:

1. Define a Hamiltonian in an extended `({π̃}, {q})` phase space:

<!-- eq:A1 -->
$$\tilde{\mathscr{H}}(\{\tilde{\pmb{\pi}}\},\{q\}) = \mathscr{E}(\{q\}) + \sum_{i} \tilde{\pmb{\pi}}_{i}^{2}/2m_{i}$$

(constant coefficients `m_i`; in Duane et al. taken as unity).

2. From current configuration `{q0}`, generate initial momenta `{π̃0}` according to

<!-- eq:A2 -->
$$P_{\rm G}(\{\tilde{\pi}_0\})[\mathrm{d}\tilde{\pi}_0] \propto \exp\left(-\sum_i (\tilde{\pi}_0)_i^2/2m_ik_{\rm B}T\right)[\mathrm{d}\tilde{\pi}_0]$$

3. Integrate the Hamiltonian equations of motion (e.g. leap-frog) for `N_MD`
time-steps of length `δt_MD`, producing candidate `{qn}` and momenta `{π̃n}`.

4. Accept `{qn}` with probability

<!-- eq:A3 -->
$$P_{\mathbf{A}} = \min\{1, \exp(-\beta \Delta \widetilde{\mathscr{H}})\}$$

where `β = 1/kB T` and
`ΔH̃ = H̃({π̃n},{qn}) - H̃({π̃0},{q0})`.

5. Go to 2.

This generates the distribution

<!-- eq:A4 -->
$$P[dq] \propto \exp(-\beta \mathscr{E}(\{q\}))[dq]$$

in `{q}` space. Since we want `P[dx] ∝ exp(-βE({x}))[dx]` in Cartesian space -
equivalent to `P[dq] ∝ D({q}) exp(-βE({q}))[dq]` in `{q}` space with `D({q})` from
eq (6) - the acceptance probability is corrected:

<!-- eq:A5 -->
$$P_{\mathbf{A}} = \min \left\{ 1, \exp\left(-\beta \Delta \widetilde{\mathcal{H}}\right) D(\left\{q_{n}\right\}) / D(\left\{q_{0}\right\}) \right\}$$

The factor `D({qn})/D({q0})` has nothing to do with HMC but is present for any MC
algorithm executing changes in `{q}` space.

The only additional role the fictitious momenta can play in the acceptance
criterion is altering the phase-space volume `[dq][dπ̃]` during step 3. The
leap-frog scheme conserves this volume: `[dq0][dπ̃0] = [dqn][dπ̃n]` (not true in
general for other discretizations). The procedure relies on a kinetic energy
*diagonal* in the fictitious momenta; a more general form (e.g.
`Σ_ij c_ij π̃_i π̃_j`) may alter the phase-space volume, requiring an additional
Jacobian factor in the acceptance probability.
