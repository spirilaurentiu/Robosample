# Self-Guided Langevin Dynamics via Generalized Langevin Equation

**Authors:** Xiongwu Wu, Bernard R. Brooks, Eric Vanden-Eijnden
**Venue:** J Comput Chem 2016, 37(6), 595–601. DOI: 10.1002/jcc.24015

## Abstract

Self-guided Langevin dynamics (SGLD) is a molecular simulation method that
enhances conformational search and sampling via acceleration of the low
frequency motions of the system. This acceleration is produced via introduction
of a guiding force which breaks down the detailed-balance property of the
dynamics, implying that some reweighting is necessary to perform equilibrium
sampling. Here we eliminate the need of reweighting and show that the *NVT* and
*NPT* ensembles are sampled exactly by a new version of self-guided motion
involving a generalized Langevin equation (GLE) in which the random force is
modified so as to restore detailed-balance. Through the examples of alanine
dipeptide and argon liquid, we show that this SGLD-GLE method has enhanced
conformational sampling capabilities compared to regular Langevin dynamics (LD)
while being of comparable computational complexity. In particular, SGLD-GLE is
fully size extensive and can be used in arbitrarily large systems, making it an
appealing alternative to LD.

**Keywords:** Self-guided Langevin dynamics; generalized Langevin equation;
molecular simulation; conformational sampling; canonical ensemble

## Introduction

Self-guided Langevin dynamics (SGLD) was developed for efficient conformational
exploration and sampling via selective acceleration of the low frequency modes
of the molecules. SGLD is unique in that this acceleration is achieved without
modifying the energy surface or raising temperature. SGLD has been applied to
many studies of events arising on long time scales (protein conformational
reorganization, rotameric substates, conformational transitions, denatured
states).

Because the guiding force involves a running average over the momentum at
previous times, the detailed-balance property of the dynamics is broken in
SGLD. In other words, SGLD is a non-equilibrium sampling method. The stationary
distribution of the method (the "SGLD ensemble") is not known explicitly. In
particular, to perform canonical (*NVT*) or isothermal-isobaric (*NPT*) sampling
via SGLD simulations requires reweighting.

The main result of the present work is to modify the SGLD equations of motion in
a way that restores the detailed-balance property of the dynamics via a suitable
modification of the random force. This new SGLD-like equation is a generalized
(non-Markovian) Langevin equation (GLE) which can be shown to sample exactly the
*NVT* and *NPT* ensembles. This method, referred to as SGLD-GLE, alleviates the
need for reweighting and is therefore fully size extensive, as are molecular
dynamics (MD) and Langevin dynamics (LD).

The idea of using a GLE to selectively accelerate conformational exploration is
not new: it is at the core of the method of Ceriotti et al. One contribution of
the present work is to make the connection between these GLE-based methods and
the general framework of SGLD. SGLD-GLE is more efficient than standard LD and
with a comparable computational cost.

## Theory and methods

We first recall the main equations of SGLD and then derive those of SGLD-GLE.
For brevity, the derivations focus on systems in the canonical (*NVT*) ensemble.
Similar developments hold for the *NPT* ensemble.

### Self-guided Langevin dynamics

The equation of motion for SGLD has the form of eq. (1), where $\dot{\mathbf{p}}_i$
and $\mathbf{f}_i$ are the time derivative of momentum and the interaction force
of particle *i*. $\mathbf{R}_i$ is a zero-mean white-in-time Gaussian random
force whose covariance (eq. 2) involves the mass $m_i$, the collision frequency
$\gamma$, and the simulation temperature $T$.

Compared to standard Langevin dynamics (LD), eq. (1) has an extra term called
the guiding force $\mathbf{g}_i$ given by eq. (3). A "$\sim$" cap denotes the low
frequency portion, given by a local exponential average (eq. 4) which can be
calculated efficiently via a simple update of the current value. The local
averaging acts like a low-frequency filter that reduces the high frequency
components of the motion while keeping its low frequency contributions.

In eq. (3), the parameter $\lambda$ is the guiding factor that controls the
strength of the guiding force and the parameter $\xi$ is an energy conservation
factor used to cancel any energy input from the guiding force (eq. 5). Solving
eq. (5) gives $\xi$ (eq. 6). The guiding forces defined by eqs. (3) and (6)
produce no net work on a simulation system since they act in a direction
orthogonal to the momentum. The $\lambda\gamma\tilde{\mathbf{p}}_i(t)$ term
accelerates the low frequency motions and the $-\xi\gamma\mathbf{p}_i(t)$ term
damps the high frequency ones.

### SGLD-GLE method

Next we modify the random force in the SGLD eq. (1) to put it in the form of a
generalized Langevin equation (GLE), eq. (7). Here $K(t-\tau)$ is the memory
kernel and $\eta_i(t)$ is a zero-mean Gaussian noise whose covariance is related
to the kernel by the fluctuation-dissipation theorem (eq. 8).

If we neglect the energy conservation term proportional to $\xi$ in the SGLD
guiding force (eq. 3), the dissipation term in the GLE (eq. 7) reduces to that in
SGLD (eq. 1) for the specific kernel choice eq. (9), with the convention eq. (10)
that the half-sided delta integrates to $\tfrac{1}{2}P(t)$.

Consistently, this kernel implies that the noise term should be eq. (11), where
$\mu \ge 0$ is a parameter related to the guiding factor $\lambda$. Substituting
into the kernel definition gives eq. (12), which reduces to eq. (9) when eq. (13)
holds: $\lambda = \mu(2-\mu)$. Eq. (13) has two roots which lead to noises that
are statistically equivalent and can be used interchangeably.

Substituting eqs. (9), (11), and (13) into eq. (7), we rewrite the equation of
motion as eq. (14), where the guiding force now also contains a random component
(eq. 15). Eq. (14) is the equation of motion used in SGLD-GLE.

Because SGLD-GLE satisfies detailed-balance, it exactly preserves the ensemble
distribution. This can be checked by noting that the guiding force satisfies the
Markovian auxiliary equation (eq. 16). Eqs. (14) and (16) form a closed system of
Markovian equations, and the Fokker-Planck equation (FPE) for the equilibrium
density $\rho(\{\mathbf{r}_i\},\{\mathbf{p}_i\},\{\mathbf{g}_i\})$ reads eq. (17).

### Derivation (not implemented): stationary distribution

One can check by direct substitution that the solution to the FPE (eq. 17) is
the extended distribution eq. (18), where $C$ is a normalization constant and
$E_p$ is the potential energy. Marginalizing over the auxiliary guiding-force
variables $\{\mathbf{g}_i\}$ gives eq. (19), the canonical Maxwell-Boltzmann
distribution. This shows that eq. (14) samples the canonical distribution
exactly.

### Derivation (not implemented): FPE proof (Appendix)

To prove that eq. (18) solves the FPE (eq. 17), denote the state vector eq. (A1)
and the extended energy $w(z)$ (eq. A2). Define the diffusion structure via
$\sigma$ and $\omega = \sigma\sigma^T/kT$ (eqs. A3, A4) and an antisymmetric
matrix $\kappa = -\kappa^T$ (eq. A5). In these terms the GLE (eq. 7) becomes the
compact form eq. (A6), with $\eta(t)$ a vectorial white noise. The associated FPE
is eq. (A7). To check that $\rho = e^{-w/kT}$ solves it, note two facts:

$$\omega \nabla_z w\,\rho + kT\,\omega \nabla_z \rho = \omega\nabla_z w\,\rho + kT\,\omega\left(-\tfrac{1}{kT}\nabla_z w\,\rho\right) = 0 \quad (A8)$$

$$\nabla_z \cdot \kappa \nabla_z w\,\rho = \kappa \nabla_z\nabla_z w\,\rho - \tfrac{\kappa}{kT}\nabla_z w\,\nabla_z w\,\rho = 0 \quad (A9)$$

the second vanishing because $\kappa$ is antisymmetric.

## Illustrative simulations

### Alanine dipeptide

Alanine dipeptide's conformation is mainly characterized by two dihedral angles,
$\phi$: CT-N-C$\alpha$-C and $\psi$: N-C$\alpha$-C-NT. High frequency motions
(bond stretching) coexist with low frequency motions ($\phi,\psi$ dihedral
changes). The CHARMM all-atom force field was used. A distance-dependent
dielectric constant of $4r$ represented solvent screening; non-bonded
interactions were calculated without a cutoff.

Simulations used a time step of 2 fs and SHAKE to fix bond lengths. Each
simulation lasted 20 ns; conformations every 2 ps were saved. All simulations
were at 300 K except high-temperature LD. A collision frequency of 10/ps was used
for all simulations. SGLD and SGLD-GLE used a local average time $t_L = 0.2$ ps.

The $\phi$-$\psi$ distribution from LD at 300 K shows two major peaks:
I:$(-90°,-70°)$ and II:$(-90°, 160°)$. The frequency of transitions between the
two peaks measures conformational sampling efficiency. The transition between I
and II has a small energy barrier.

Elevated temperature shifts conformational sampling toward high energy
conformations. In SGLD (without removing the net guiding force, for comparison),
the guiding effect shifts sampling toward lower energy conformations because more
kinetic energy is distributed to the center of mass. Compared with
high-temperature LD, SGLD increases the number of transitions much more
effectively while causing much less deviation in average potential energy. For
SGLD-GLE, with $\lambda$ increasing from 0 to 1, the number of transitions
increases from 353 to 515 while the potential energy remains almost constant.
SGLD-GLE at $\lambda=1$ does much better than LD and even much better than
high-temperature LD.

High-temperature LD and SGLD have flattened $\psi$-angle distributions, while
SGLD-GLE produces almost exactly the same distributions as LD. The potential
energy distribution of SGLD-GLE almost overlaps with the LD result, confirming
that SGLD-GLE samples the canonical distribution exactly.

### Liquid argon

An argon fluid tested SGLD-GLE in the *NPT* ensemble. Argon atoms are described
by the Lennard-Jones 6-12 potential with $\varepsilon = 119.8$ K and
$\sigma = 3.405$ Å. 500 argon atoms were placed in a cubic periodic box
($28.53 \times 28.53 \times 28.53$ Å$^3$). A collision frequency of 10 ps$^{-1}$
and a time step of 1 fs were used. SGLD and SGLD-GLE simulations were at 100 K.
The target pressure was 1 atm. Each simulation lasted 10 ns; coordinates and
velocities every 0.05 ps were stored.

High-temperature LD significantly shifts sampling to high energy conformations,
while SGLD-GLE has an energy distribution almost identical to LD. For SGLD, at
high friction constant and guiding factor, the conformational distribution also
shifts toward high energy.

Accelerated conformational sampling corresponds to an increase in the diffusion
constant. Temperature elevation caused fastest energy and volume increases; SGLD
caused smaller changes; SGLD-GLE resulted in almost identical average energies
and volumes to LD while achieving significantly larger diffusion constants.

The dynamic property is understood from the spectrum of the velocity
autocorrelation function, computed via eq. (20). As temperature increases from
100 K to 130 K, the spectrum shifts upward at all frequencies. SGLD and SGLD-GLE
at the same temperature (100 K) enhance the slow motions (low frequencies) and
reduce fast motions (high frequencies). A longer averaging time (larger $t_L$)
further enhances slow motions. The guiding factor determines the increase in the
diffusion constant, while the averaging time determines the portion of motion to
be enhanced.

## Concluding remarks

SGLD-GLE is a modification of SGLD in which a generalized Langevin equation is
used to restore the detailed-balance property of the dynamics. This permits
sampling the *NVT* and *NPT* ensembles exactly, without the reweighting needed in
SGLD. Since the guiding force (eq. 15) is proportional to the friction constant
$\gamma$, the improvement over LD is more pronounced when $\gamma$ is large;
SGLD-GLE becomes comparable to LD in the limit $\gamma \to 0$.

The modification restoring detailed balance typically reduces the overall
acceleration compared to SGLD. Thus SGLD may be preferred if conformational
exploration is the main objective, while SGLD-GLE is the appealing alternative to
LD if unbiased conformational sampling is the goal.

SGLD-GLE has several unique characteristics. First, it is size extensive: the
guiding force is calculated from momentum and random forces, independent of
system size; reweighting is no longer needed, making it suitable for million-atom
systems. Second, the method can be applied to a part of a simulation system (the
guiding factor can be defined per-atom, like the friction constant). Third, the
motion mode to be enhanced can be controlled by the local averaging time.
