# Enhanced sampling methods

1. **Tempering / generalized-ensemble methods (change the temperature or ensemble)**

    Intuition: heat the system (or several copies of it) so it has enough thermal energy to jump barriers, then recover the correct room-temperature statistics afterward. These are CV-free methods since you don't have to know in advance what motion matters.

    - **Parallel tempering / temperature replica exchange (T-REMD)**: run many copies at a ladder of temperatures and periodically swap configurations between neighbors via a Metropolis criterion. Hot copies cross barriers; cold copies stay physical; swapping shuttles good configurations down to the temperature you care about.
    - **Simulated tempering:** like the above but a single copy random-walks up and down the temperature ladder itself.
    - **Hamiltonian replica exchange / REST & REST2 (replica exchange with solute tempering):** instead of heating everything, you only *heat* the part you care about (e.g., the solute), which makes the method scale to large solvated systems far better. Solute tempering methods exist because pure temperature methods scale badly with system size (you need more and more replicas).
    - **Integrated tempering sampling (ITS):** sums Boltzmann factors over a range of temperatures into one effective potential.

2. **Collective-variable-based biasing (change the energy surface along a chosen coordinate)**

    Intuition: you guess one or a few collective variables (CVs) e.g. reduced coordinates like a distance, a dihedral, or a fancy combination that capture the slow motion, then add a bias that pushes the system along them or flattens the barrier. Afterward you *reweight* to remove the bias and recover the true free-energy landscape. A CV is a function of input coordinates that offers a simplified description of the system's structure, and there are typically far fewer of them than total degrees of freedom. A *bad* reaction coordinate hides the real barrier and gives wrong or unconverged results.

    - **Umbrella sampling:** place a series of harmonic "windows" along the CV so the system is forced to visit each region, then stitch the windows back together (usually with WHAM).
    - **Metadynamics:** periodically deposit Gaussian *hills* of bias at the system's current location, progressively filling up the basin it's in until it's forced to explore elsewhere. Well-tempered metadynamics makes the hill height shrink over time so the bias converges cleanly.
    - **Adaptive biasing force (ABF):** estimates the mean force along the CV on the fly and cancels it, so the system diffuses freely along that coordinate.
    - **Variationally enhanced sampling (VES):** poses the optimal bias as a variational minimization problem.
    - **Steered MD and temperature-accelerated MD (TAMD):** pull the CV mechanically, or couple it to a fictitious high-temperature variable, respectively.

3. **CV-free boosting / surface-flattening (change the whole energy surface, no CV needed)**

    Intuition: mechanically raise the low-energy regions of the entire potential so barriers shrink everywhere, without committing to any particular coordinate. No reaction coordinate is required, but the boost is *dumb* (it doesn't know where the interesting barrier is), so reweighting can be noisy for high barriers.

    - **Accelerated MD (aMD):** adds a boost potential that lifts the surface wherever the energy is below a threshold.
    - **Gaussian accelerated MD (GaMD):** A boost potential with a nearly Gaussian distribution is applied whenever the system's potential energy falls below a predefined threshold, and a cumulant expansion is used to reconstruct unbiased thermodynamic averages. The Gaussian shape is what makes the reweighting tractable, which was aMD's main pain point.

4. **Path / trajectory-based methods (change the trajectories themselves)**

    Intuition: don't touch the physics at all. Instead, run many trajectories and cleverly manage them so computer time is spent on the transition rather than the waiting. As one path-sampling review frames it, you focus the computing effort on functional transitions rather than stable states. This family is special because it preserves unbiased kinetics (rate constants), which the bias-based families generally distort.

    - **Transition path sampling (TPS):** Monte Carlo in the space of trajectories. It's essentially a Markov chain Monte Carlo algorithm in trajectory space that harvests reactive pathways by perturbing existing ones (throwing ropes over rough mountain passes in the dark).
    - **Transition interface sampling (TIS) and forward flux sampling (FFS):** place interfaces between reactant and product and compute the crossing probability stage by stage. FFS only integrates forward in time, so unlike TPS it does not rely on reversibility of the dynamics and handles non-equilibrium systems.
    - **Milestoning:** put "milestones" across the landscape and compute transition statistics between adjacent ones.
    - **Weighted ensemble (WE):** runs many weighted parallel walkers, periodically resampling (splitting walkers in under-explored regions, merging in crowded ones) to push probability into rare regions while keeping statistics exact. It can achieve superlinear scaling: unbiased estimation of observables such as rate constants and equilibrium populations to greater precision than ordinary parallel simulation. **Adaptive multilevel splitting (AMS)** is a close mathematical cousin.
