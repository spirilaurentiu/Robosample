# Alchemical Grid Dock (AlGDock): Binding Free Energy Calculations between Flexible Ligands and Rigid Receptors

**Author:** David D. L. Minh (Illinois Institute of Technology)
**Source:** arXiv 1507.03703v2 (2019). Open-source, MIT license: https://github.com/ccbatiit/algdock/

## Abstract

Alchemical Grid Dock (AlGDock) is open-source software designed to compute the binding potential of mean force (BPMF) - the binding free energy between a flexible ligand and a rigid receptor for a small organic ligand and a biological macromolecule. Multiple BPMFs can be used to rigorously compute binding affinities between flexible partners. AlGDock uses replica exchange between thermodynamic states at different temperatures and receptor-ligand interaction strengths. Receptor-ligand interaction energies are represented by interpolating precomputed grids. Thermodynamic states are adaptively initialized and adjusted on-the-fly to maintain replica exchange rates. In demonstrative calculations, when the bound ligand is treated as fully solvated, AlGDock estimates BPMFs with a precision within 4 kT in 65% and within 8 kT for 91% of systems. It correctly identifies the native binding pose in 83% of simulations. Performance is sometimes limited by subtle differences in the important configuration space of sampled and targeted thermodynamic states.

## Introduction

AlGDock computes the binding potential of mean force (BPMF) - the binding free energy between a flexible ligand and a rigid receptor - between a small organic ligand and a biological macromolecule.

The BPMF is defined as a ratio of configurational integrals,

<!-- eq:1 -->
$$B(r_R) = -\beta^{-1} \ln \left( \frac{\int I(\xi) J(\xi) e^{-\beta U(r_{RL})} \, dr_L \, d\xi}{\int I(\xi) J(\xi) e^{-\beta [U(r_L) + U(r_R)]} \, dr_L \, d\xi} \right).$$

The internal coordinates (excluding translation and rotation) of a receptor-ligand complex, $r_{RL}$, are partitioned into the receptor, $r_R$, the ligand, $r_L$, and the relative translation and rotation of the species, $\xi$. $\beta = (k_B T)^{-1}$ is the inverse of Boltzmann's constant times the temperature. $I(\xi)$ is an indicator function that specifies whether the receptor and ligand are bound (1) or not (0). $J(\xi)$ is the Jacobian for transforming Cartesian coordinates into the coordinate system used for $r_L$ and $\xi$. $U(\cdot)$ is the potential energy of a species in solvent.

According to implicit ligand theory (ILT), the standard binding free energy can be computed from BPMFs between a ligand and multiple receptor conformations. ILT also allows BPMFs to reweight receptor conformations from the apo (ligand-free) to the holo (ligand-bound) ensemble.

The main methodological distinctions of AlGDock are (1) the use of precomputed nonbonded interaction grids for receptor-ligand interactions, and (2) the adaptive initialization and on-the-fly adjustment of thermodynamic states. Grids remove the $O(N^2)$ scaling of nonbonded evaluation with receptor atom count $N$; once the grid is computed, calculation time does not depend on $N$. Adaptive states improve free energy precision by ensuring sufficient configuration space overlap between adjacent thermodynamic states along the alchemical protocol.

AlGDock is a python module based on the Molecular Modeling Toolkit (MMTK) 2.7.8.

## Methodology

### Thermodynamic Cycle

BPMFs are calculated based on a thermodynamic cycle whose milestone thermodynamic states are labeled A to E. The states between and including milestones X and Y are referred to as states XY. Over the cycle, the receptor-ligand interaction strength is scaled and the temperature is varied (high-temperature states enhance transitions between local energetic minima).

Because states are at different temperatures, the paper works with the **reduced potential energy** $u$, a log probability density that incorporates the inverse temperature $\beta = (k_B T)^{-1}$. The reduced free energy difference between milestones X and Y is denoted $f_{XY}$. Converting reduced to standard quantities divides by $\beta$; units of reduced quantities are $k_B T$.

Milestone temperatures: target $\beta_T^{-1} = k_B(300\,\text{K})$, high $\beta_H^{-1} = k_B(600\,\text{K})$. $U_T(\cdot)$ and $U_S(\cdot)$ are potential energies for the target and sampling force fields, respectively (including MM terms and implicit solvent). $\Psi_g(\cdot)$ is the potential energy due to receptor-ligand interaction grids.

In all simulated states, the ligand is confined to the binding site using a flat-bottom harmonic potential,

<!-- eq:2 -->
$$u_I(d) = \begin{cases} 0 & \text{if } d \le d_0 \\ \frac{1}{2}\beta k(d - d_0)^2 & \text{if } d > d_0 \end{cases},$$

where $k = 10000$ kJ/(mol nm$^2$) is the spring constant, $d$ is the distance between the ligand center of mass and the center of the binding site, and $d_0 = 6.0$ Å is the binding site radius. There is no restriction on ligand rotation.

**Sampling vs. target force fields.** Both use AMBER ff14SB for proteins/ions and GAFF2 with AM1BCC charges for other molecules. Differences: solvation model and whether receptor-ligand interactions are evaluated directly (target, via OpenMM) or by grid interpolation (sampling, via MMTK). The sampling force field uses the generalized Born/surface area model II (OBC) from Onufriev et al., adapted from OpenMM, with only the ligand assumed solvated. Two solvation pathways: **Desolvated** (implicit solvent scaled down to zero at milestone C; no solvent for CD) and **Full** (implicit solvent at full strength for BD). Both should agree in the limit of asymptotic sampling.

**Grid interaction energy** (one electrostatic + two van der Waals grids):

<!-- eq:3 -->
$$\Psi_g(r_{RL}) = \Psi_{PBSA}(r_{RL}) + \Psi_{vdW}(r_{RL}).$$

Electrostatic $\Psi_{PBSA}$: atomic partial charges times the electrostatic potential obtained by trilinear interpolation of a precomputed grid, produced by solving the linear Poisson-Boltzmann equation around the minimized receptor with APBS 1.4 (sequential focusing). Fine grid spacing 0.5 Å; protein dielectric 2.0, solvent dielectric 80.0, solvent radius 1.4 Å, temperature 300.0 K.

vdW $\Psi_{vdW}$: grid-based with a nonlinear transformation, trilinear interpolation, inverse transformation. Inverse-transformation power 4 for the repulsive potential; no transformation for the attractive potential.

**Soft grids and end-point catastrophe.** Between milestones C and D a set of soft Lennard-Jones repulsive and electrostatic grids is introduced. The original grid value $v_o$ is replaced with $v_{max}\tanh(v_o/v_{max})$. For the soft LJ repulsive grid $v_{max} = 10.0$ kJ mol$^{-1/2}$. To avoid electrostatic pinning, $v_{max}$ for the electrostatic grid is set to 10 times the minimum ratio of Lennard-Jones and electrostatic scaling factors.

The reduced potential energy is switched according to the protocol (progress variable $\alpha$, consistent with milestone C at $\alpha=0$ and milestone D at $\alpha=1$),

<!-- eq:4 -->
$$u_{\alpha}(r_{RL}) = \frac{1}{k_B T(\alpha)} \left[ U_s(r_L) + U_s(r_R) + \alpha_{sg}(\alpha)\,\Psi_{sg}(r_{RL}) + \alpha_g(\alpha)\,\Psi_g(r_{RL}) \right]$$

$$\alpha_{sg}(\alpha) = -(2\alpha - 1)^2 + 1$$

$$\alpha_g(\alpha) = \frac{(2\alpha - 1)^2}{1 + \exp\left[-1000(\alpha - \tfrac{1}{2})\right]}$$

$$T(\alpha) = (T_T - T_H)\alpha + T_H$$

This turns on the soft grids first, then the unperturbed grids.

### Sampling

For states BD, ligand conformational sampling combines Hamiltonian Monte Carlo (HMC), external coordinate Markov chain Monte Carlo (MCMC) moves, and Hamiltonian replica exchange.

**HMC:** trial moves based on 50 steps of velocity Verlet molecular dynamics with an adaptive time step between 0.1 and 5.0 fs.

**External coordinate MCMC** (states CD, only when $\alpha < 0.01$): (1) random rotation - a random quaternion converted to a rotation matrix applied about the center of mass; (2) random translation - each dimension drawn from a Gaussian with standard deviation 0.6 Å. Accepted/rejected by the Metropolis criterion. Not attempted for $\alpha > 0.01$ due to low acceptance.

**Hamiltonian replica exchange** (states BC and CD during production): a MCMC move that swaps configurations of a pair of simulations at different states $a$ and $b$ with reduced energies $u_a$ and $u_b$. If $x$ is the original configuration in state $a$ and $y$ in state $b$, the acceptance probability

<!-- eq:5 -->
$$p_{acc} = \min \left[ 1, \; e^{-u_a(y) - u_b(x) + u_a(x) + u_b(y)} \right]$$

preserves the Boltzmann distribution in both states. As replica exchange is a type of Gibbs sampling, arbitrary pairs may be attempted. Each sweep attempts swaps between pairs $1, 2, \ldots, \min(5, K)$ states apart, where $K$ is the total number of states in the direction.

### Stages

1. **Ligand preparation:** 5000 steepest descent minimization steps; temperature ramped 20 K -> 300 K over 30 geometrically spaced simulations of 2500 steps each.
2. **Initialization:** from 50 seed configurations, 2000-step simulations initialize each thermodynamic state for states BD.
3. **Equilibration and production:** conducted for states BC then states CD.
4. **Postprocessing:** milestone B and D samples postprocessed with the target force field.
5. **Estimation:** free energy differences summing to the BPMF are estimated.

### Initialization

Goal: establish a protocol with reasonable time steps and mean replica exchange rate $\langle p_{acc}\rangle$ between all neighboring states. Low exchange probabilities indicate poor configuration space overlap, limiting free energy convergence.

For states BC, the first state ($k=0$) is at 300 K; the ligand is warmed to 600 K. For states CD, the first state depends on whether a fully-bound pose is available: if so, first state is fully bound at 300 K, otherwise fully unbound at 600 K (50 configurations from milestone C placed at 453 random center-of-mass positions - 0.5 positions per Å$^3$ - and 100 random orientations).

After state $k$ is initialized, state $k+1$ is initialized:

**1. Parameter selection** using thermodynamic length. Thermodynamic length is a metric of distance on the manifold of thermodynamic states. Statistical error is minimized (and exchange frequency nearly maximized) when intermediate states are equidistant in thermodynamic length. With parameters $\lambda$ (components $\lambda^i$) and $\gamma \equiv \gamma(\alpha)$ describing the dependence of $\lambda$ on $\alpha$ (so $\gamma(0)$ initial, $\gamma(1)$ final), the thermodynamic length is the path integral

<!-- eq:6 -->
$$\mathcal{L} \equiv \int_0^1 \sqrt{\sum_{i,j} \frac{\partial \gamma^i}{\partial \alpha} \, g(\gamma)_{ij} \, \frac{\partial \gamma^j}{\partial \alpha}} \; d\alpha.$$

The reduced potential energy is $u_\lambda(x) = U_\lambda(x)/(k_B T_\lambda)$. The normalized log probability of observing $x$ is $l_\lambda(x) = -u_\lambda(x) - \ln Z_\lambda$, where $Z_\lambda = \int e^{-u_\lambda(x)}dx$ is the partition function. Elements of the Fisher information matrix are

<!-- eq:7 -->
$$g(\gamma)_{ij} \equiv \sigma_\lambda^2 \left[ \partial^i l_\lambda, \; \partial^j l_\lambda \right],$$

where $\sigma_\lambda^2$ is the covariance in state $\lambda$ and $\partial^i$ denotes a partial derivative with respect to $\lambda^i$. For a protocol in which only one parameter $\lambda^i$ varies with $\alpha$, the length is $\mathcal{L} = \int_0^1 \frac{\partial \gamma^i}{\partial \alpha}\,\sigma_\lambda\left[\partial^i l_\lambda\right]\,dt$.

A simple single-parameter approximation is $\mathcal{L} = \Delta\lambda^i\,\sigma_0[\partial^i l_\lambda]$, with $\sigma_0$ a standard deviation in the initial state. To keep $\mathcal{L}$ approximately constant between stages, the parameter change should be inversely proportional to $\sigma_0[\partial^i l_\lambda]$,

<!-- eq:8 -->
$$\Delta\lambda^i = \frac{s}{\sigma_0\left[\partial^i l_\lambda\right]},$$

where $s$ is the adjustable **thermodynamic speed**.

For the Full solvation pathway, between milestones B and C the varying parameter is temperature $T$. With $l_\lambda = -\frac{U_S(r_L)}{k_B T} - \ln Z_\lambda$, $T$ is incremented by

<!-- eq:9 -->
$$\Delta\lambda^i = -\frac{s_{bc}\,k T^2}{\sigma_\lambda[U_S]},$$

where $s_{bc} = 20.0$. For states CD, with $l_\lambda = -u_\alpha(r_{RL}) - \ln Z_\lambda$, $\alpha$ is incremented by

<!-- eq:10 -->
$$\Delta\lambda^i = s_{cd} \left[ \left| \frac{d\alpha_{sg}}{d\alpha} \right| \frac{\sigma_\lambda[\Psi_{sg}]}{k_B T(\alpha)} + \left| \frac{d\alpha_g}{d\alpha} \right| \frac{\sigma_\lambda[\Psi_g]}{k_B T(\alpha)} + |T_T - T_H| \frac{\sigma_\lambda[u_\alpha(r_{RL})]}{T(\alpha)} \right]^{-1},$$

with $s_{cd} = 0.2$. If the targeted parameter value is exceeded (e.g. $T$ above 600 K), the targeted value is used.

**2. Seed selection.** 50 configurations from state $k$ are resampled as starting seeds for state $k+1$ with weights proportional to $\exp[u_k(x_i) - u_{k+1}(x_i)]$, where $u_k(x) = U_k(x)/(k_B T_k)$. This is sampling importance resampling; in the limit of infinite sampling of state $k$, resampled configurations are Boltzmann-distributed in state $k+1$.

**3. Sampling and adaptation.** 2000-step simulations from each seed; the time step is adapted: if MC acceptance rate > 0.8, increase step by 0.125 fs; if < 0.4, reduce by 0.25 fs; if < 0.1, reduce by 0.5 fs. Repeat until acceptance is between 0.4 and 0.8.

**4. Verification.** The mean replica exchange probability $\langle p_{acc}\rangle$ (sample mean of $p_{acc}$, Eq. 5, over all pairs of initial samples at the same time index from states $k$ and $k+1$) verifies states are neither too distinct nor too similar. If too low (< 0.4), state $k+1$ parameters are reselected with a smaller increment (thermodynamic speed multiplied by 4/5) and simulations repeated. If too high (> 0.99), state $k$ is removed.

### Equilibration and production

Broken into cycles. Each cycle: 1000 iterations of {one HMC move + 20 external coordinate MCMC moves (if $\alpha < 0.01$) per state}, then 25 sweeps of replica exchange. 50 snapshots saved per replica exchange cycle. Demonstrative runs: 8 cycles for states BC, 15 cycles for states CD.

Between cycles, thermodynamic states are inserted if the average replica exchange acceptance rate between any neighboring pair falls below 0.4; the new state is populated by sampling importance resampling.

Equilibration/production separation follows Chodera: integrated autocorrelation time and statistical inefficiency are estimated from the mean potential energy of configurations from the last $c \in \{1, 2, \ldots, C\}$ cycles; the number of independent samples is the number of snapshots in $c$ cycles divided by the statistical inefficiency; equilibration is set by the $c$ giving the largest number of independent samples.

### Estimation

BPMFs are estimated as $\beta_T B(r_R) = f_{AB} + f_{BC,L} + f'_{CD} + f_{DE}$. $f_{BC,L}$ is the free energy of warming the ligand from $T_T = 300$ K to $T_H = 600$ K, and

<!-- eq:11 -->
$$f'_{CD} = -\ln \frac{\int I(\xi)J(\xi)e^{-\beta_T[U(r_L)+\Psi_g(r_{RL})]}\,dr_L\,d\xi}{\int I(\xi)J(\xi)e^{-\beta_H U(r_L)}\,dr_L\,d\xi}.$$

$f_{BC,L} + f'_{CD}$ is used instead of $f_{BC} + f_{CD}$ because it avoids determining the receptor internal energy $U(r_R)$. The receptor desolvation free energy is $f_{AB,R} = \beta_T(U(r_R) - U(r_R))$. $f_{AB,R}$ and $f_{DE}$ are estimated by free energy perturbation (Zwanzig) using configurations from milestones A and E. $f_{BC,L}$ and $f'_{CD}$ are estimated by the multistate Bennett acceptance ratio (MBAR), which uses potential energies from every replica.

### Pose prediction

Configurations from milestone D are clustered by hierarchical clustering with complete linkage (scipy.cluster.hierarchy.linkage). Distances use the Hungarian symmetry-corrected heavy-atom RMSD (as in UCSF DOCK 6). Cluster separation threshold 1.0 Å. Each cluster's probability is obtained by reweighting configurations by

<!-- eq:12 -->
$$w_c = \exp\left[-\beta_T \left(U_T(r_{RL}) - U_S(r_L) - \Psi_g(r_L)\right)\right],$$

or, assuming interaction energies are the only terms that change between milestones D and E,

<!-- eq:13 -->
$$w_c = \exp\left[-\beta_T \left(U_T(r_{RL}) - U_T(r_L) - \Psi_g(r_L)\right)\right].$$

Reduced free energy of each pose $p$ from the cumulative cluster weight:

<!-- eq:14 -->
$$f_{EE_p} = -\ln \frac{\sum_c w_c}{\sum_p \sum_c w_c},$$

where $\sum_c$ is over configurations in the cluster and $\sum_p$ over poses. The pose-specific BPMF is

<!-- eq:15 -->
$$f_{AE,p} = f_{AE} + f_{EE_p}.$$

### Astex diverse set BPMF calculations

For each system, 11 independent simulations were performed with the Desolvated and Full solvation pathways. AMBER input files (ff14SB + GAFF2 + AM1BCC) reused from a previous study. Simulations started from crystallographic and docked poses (UCSF DOCK 6, minimum anchor size 5). Starting poses minimized for 1000 conjugate gradient steps for milestone D. Benchmark on XSEDE Comet (Intel Xeon E5-2680v3, single-core jobs).

## Results

### Thermodynamic state initialization is system-specific and robust

The number of thermodynamic states $N_{states}$ varies widely (evidence of system specificity). For states BC: Desolvated 67-182, Full 51-111 (more states for Desolvated due to solvent removal at milestone C). For states CD, comparable between pathways. Robustness: $\sigma[N_{states}]$ small relative to $\bar{N}_{states}$; for states BC $\sigma[N_{states}] < 2$ for all systems. Time steps converged to 2.75-3.75 fs in the vast majority of initializations.

### Replica exchange acceptance probabilities are reasonable

For states BC, $\bar{p}_{acc}$ estimated during replica exchange are high with low variance; all $\ge 0.92$. For states CD, high but larger variance (most between 0.7 and 1.0), with occasional dips around $\alpha = 0.2$ or $\alpha = 0.8$, none low enough to be a bottleneck.

### Simulation time

Desolvated: 4-44 hours; Full: 6-117 hours. Roughly exponential dependence on ligand size; roughly linear up to 40 atoms. Most time spent in equilibration and production, largest fraction in states CD.

### Convergence of free energy differences

BPMFs estimated within chemical precision (1 kcal/mol = 1.68 RT) for only 28.2% of systems; within 4 $k_B T$ for 75.3% (Desolvated) and 74.1% (Full); within 8 $k_B T$ for 87.1% (Desolvated) and 94.1% (Full). Largest source of imprecision is $f_{DE}$ (transferring the complex between sampling and target force fields).

- $f_{AB}$: converges quickly; RMSE < 1 $k_B T$ after first cycle for all systems; < 0.5 $k_B T$ after 8 cycles for all but 1p62 Desolvated (0.544 $k_B T$).
- $f_{BC}$: Desolvated RMSE < 2 $k_B T$ after 1 cycle except 1jje (3.2 $k_B T$), < 1 $k_B T$ after 8; Full within 1 $k_B T$ after 1 cycle, within 0.15 $k_B T$ after 8.
- $f_{CD}$: slower; Desolvated RMSE > 8 $k_B T$ for 23 systems after 1 cycle, < 2.5 $k_B T$ after 15 for all but 1t40 (2.65 $k_B T$); Full > 8 $k_B T$ for 28 systems after 1 cycle, < 2.5 $k_B T$ after 15 except 1l7f (3.41 $k_B T$).
- $f_{DE}$: Desolvated RMSE > 8 $k_B T$ for 26 systems after 1 cycle, still > 8 $k_B T$ for 9 systems after 15; Full > 8 $k_B T$ for 23 after 1 cycle, 5 after 15.

### False convergence

Comparing Desolvated and Full options: 63 systems (74.1%) agree within error. Consensus BPMF = lower of the two mean BPMFs. Full option is less susceptible to false convergence. Convergence limited by differences in the important configuration space of milestones D and E.

### Native pose identification

Native pose (RMSD < 2 Å from crystal) identification for the demonstrative set; the force field at milestone E performs best. GOLD 80.5%, GLIDE 82%, ICM 91% on the same Astex set; milestone-E free energy ranking here: 75.8% (Desolvated), 83.9% (Full). Full solvation option outperforms Desolvated (lower BPMFs and better pose prediction).

## Conclusions

A robust method to estimate BPMFs for protein-ligand systems. Largest sources of imprecision are configuration space overlap between representations of the complex. The single-parameter (thermodynamic speed) state initialization scheme and the automatic state-insertion rule (when exchange rate < 40%) may be useful for other classes of simulations.
