Contents lists available at [ScienceDirect](http://www.sciencedirect.com/science/journal/03044165)

# BBA - General Subjects

journal homepage: [www.elsevier.com/locate/bbagen](https://www.elsevier.com/locate/bbagen)

![](_page_0_Picture_5.jpeg)

# Robosample: A rigid-body molecular simulation program based on robot mechanics

![](_page_0_Picture_7.jpeg)

Laurentiu Spiridon[a,](#page-0-0)[⁎](#page-0-1) , Teodor Asvadur Şulea[a](#page-0-0) , David D.L. Minh[b,](#page-0-2)[⁎](#page-0-1) , Andrei-Jose Petrescu[a,](#page-0-0)[⁎](#page-0-1)

- <span id="page-0-0"></span><sup>a</sup> Department of Bioinformatics and Structural Biochemistry, Institute of Biochemistry of the Romanian Academy, Splaiul Independentei 296, Bucharest 060031, Romania
- <span id="page-0-2"></span><sup>b</sup> Department of Chemistry, Illinois Institute of Technology, Chicago, IL 60616, USA

## ARTICLE INFO

Keywords: Molecular dynamics Hamiltonian Monte Carlo Gibbs sampling Robotics Multibody dynamics

#### ABSTRACT

Background: Compared with all-atom molecular dynamics (MD), constrained MD methods allow for larger time steps, potentially reducing computational cost. For this reason, there has been continued interest in improving constrained MD algorithms to increase configuration space sampling in molecular simulations.

Methods: Here, we introduce Robosample, a software package that implements high-performance constrained dynamics algorithms, originally developed for robotics, and applies them to simulations of biomolecular systems. As in the gMolmodel package developed by Spiridon and Minh in 2017, Robosample uses Constrained Dynamics Hamiltonian Monte Carlo (CDHMC) as a Gibbs sampling move - a type of Monte Carlo move where a subset of coordinates is allowed to change. In addition to the previously described Cartesian and torsional dynamics moves, Robosample implements spherical and cylindrical joints that can be distributed along the molecule by the user.

Results: In alanine dipeptide simulations, the free energy surface is recovered by mixing fully flexible with torsional, cylindrical, or spherical dynamics moves. Ramachandran dynamics, where only the two key torsions are mobile, accelerate the slowest transition by an order of magnitude. We also show that simulations of a complex glycan cover significantly larger regions of the configuration space when mixed with constrained dynamics.

Major conclusions: Robosample is a tool of choice for efficient conformational sampling of large biomolecules. General significance: Robosample is intended as a reliable and user-friendly simulation package for fast biomolecular sampling that does not require extensive expertise in mechanical engineering or in the statistical mechanics of reduced coordinates.

# 1. Introduction

Structural biology techniques including crystallography and cryoelectron microscopy have provided us with a wealth of information about the free energy minima of biomolecular systems, but this information may not adequately describe a complex energy landscape. Processes involving flow over an energy landscape, such as folding or large conformational changes, may evade a simple structural description. Moreover, not all biomolecules are completely folded. Many contain moieties, such as intrinsically disordered regions in proteins or glycans, that are innately flexible and may adopt a continuum of configurations that can only be described statistically [\[1](#page-9-0)–4].

Over the past decade, interest in investigating highly flexible molecular species has sharply increased. It has become clear that flexibility plays a fundamental role in higher organisms where, for example, over 50% of the transcriptome do not fold in a unique native conformation and over 60% of the secreted proteins are glycosylated [\[3,](#page-9-1)[5](#page-9-2)]. Flexibility enables the structural plasticity of such molecules, allowing them to play a fundamental role in interaction with structurally diverse partners. For example, intrinsically disordered regions (IDR) may function as chaperones, effectors, assemblers, or scavengers [[6](#page-9-3)]. In their absence, complex multidomain proteins lose their ability to interact and perform their function [[7](#page-9-4)]. Similarly, glycans can mold to protein interfaces. In this way, they can modulate glycoprotein folding or important protein-protein interactions [8–[10](#page-9-5)].

While simulations can potentially describe the entirety of biomolecular energy landscapes, their accuracy is often limited by conformational sampling. Simulations have often used data from bulk experiments such as Small Angle Scattering (SAS) or Nuclear Magnetic Resonance (NMR) as restraints for modelling statistical ensembles with

E-mail addresses: [spiridon.laurentiu@biochim.ro](mailto:spiridon.laurentiu@biochim.ro) (L. Spiridon), [dminh@iit.edu](mailto:dminh@iit.edu) (D.D.L. Minh), [andrei.petrescu@biochim.ro](mailto:andrei.petrescu@biochim.ro) (A.-J. Petrescu).

0304-4165/ © 2020 Published by Elsevier B.V.

<span id="page-0-1"></span><sup>⁎</sup> Corresponding authors.

Molecular Dynamics or Monte Carlo [11–16]. With improved force fields, recent reports have shown that MD simulations are able to produce accurate physical models of an IDRs consistent with SAS and NMR measurements [17]. These reports suggest that, given enough computational power and improved sampling algorithms, the energy landscapes of flexible molecules can be accurately modeled.

One approach to enhanced sampling that has received ongoing attention is the application of holonomic constraints. A general problem with molecular dynamics integrators is that they require small time steps to model high-frequency motion. Constrained dynamics simulations present the opportunity to eliminate degrees of freedom with the highest frequency, leading to longer time steps. Some of the earliest methods, which are still widely used to constrain bonds that include hydrogen atoms, are the SHAKE and RATTLE algorithms [18,19]. Torsional dynamics simulations, which require constraints on bond lengths and bond angles, have been applied to crystallographic and NMR structure refinement [20,21]. Other applications of constrained MD include loop modelling and the rapid folding of small proteins [22].

In a constrained MD simulation, molecules may be represented as multibody chains where atoms or *groups of atoms* are rigid bodies and joints are mapped onto chemical bonds. Such representations more naturally reflect the energy as a function of internal molecular motion than Cartesian coordinates. Representations where groups of atoms are rigid bodies allow for a type of coarse graining that preserves a finegrained energy function. Whereas traditional coarse graining reduces the degrees of freedom (DOFs) in both the model and the dynamics, leading to a perturbed equilibrium distribution, a generalized coordinate formulation may drop DOFs during dynamics while preserving the energy function. For example, the phenyl group can be regarded as one rigid body. Nevertheless, all the atoms in the ring are accounted for when calculating the energy and each one of them ends up contributing to the resulting torque.

One limitation of constrained dynamics is ergodicity, the ability to access all regions of the configuration space. If simulations are performed with constraints, then samples are drawn from a conditional opposed to a joint probability distribution. To address this limitation, Vaidehi and Jain suggested freezing and releasing DOFs during MD simulations [23]. To reduce the difference between the conditional and joint probability distribution, Kandel et al. [24] relaxed the requirement of torsional dynamics by allowing bond angles to be flexible. More recently, we pursued the former direction, using constrained MD to generate candidate configurations for Markov chain Monte Carlo [25]. Such moves are an example of Gibbs sampling [26], in which a Markov chain is propagated by sampling from a conditional probability distribution. As long as all DOFs can be modified by at least one of the Monte Carlo moves, then ergodicity can, in principle, be achieved.

Here, we present Robosample, a user-friendly application that implements Hamiltonian Monte Carlo [27] based on constrained MD and fully flexible MD implemented using algorithms borrowed from robotics. Robosample is freely available at <a href="https://github.com/spirilaurentiu/Robosample">https://github.com/spirilaurentiu/Robosample</a>. It was built starting from the Molmodel API [28] and its generalization gMolmodel API [25] designed to perform Gibbs sampling in AlGDock [29] based on torsional and fully-flexible constrained dynamics.

Robosample refactored and optimized the gMolmodel code and was written to work as a standalone application, as well as an API. In addition to torsional and fully flexible dynamics, Robosample implements many joint types such as: Ball, Cylinder, Universal, Slider etc., allowing to perform various types of rigid body movements (mobilities). Of these we emphasize here the spherical and cylindrical joints that capture torsion-angle and torsion-bond couplings. Moreover, Robosample calculates the gradient of the mass matrix determinant, which allows MD trajectories to sample from the conditional probability of the Boltzmann distribution.

In this paper we present the statistical mechanics basis of Robosample, its architecture and its Gibbs sampling implementation. In addition, Robosample is tested herein on alanine dipeptide and a complex glycan with relevance to hepatitis C virus (HCV) immunity.

# 1.1. The constrained generalized coordinates molecular simulation basis of Robosample

Simulations with m constraints resort on two types of dynamics formulations: the Lagrange multiplier formulation in Cartesian coordinates, with a maximum number of 3 N DOFs and m additional m constraint equations; and the reduced coordinate formulation based on  $N_f=3$  N - m generalized coordinates obtained from a coordinate transformation. The limiting case of both formulations is when m=0 and  $N_f=3$  N. Here the number of needed Lagrange multiplier equations is 0 and the number of reduced coordinates is the maximum number of DOFs. A straightforward approach to imposing constraints in generalized coordinates starts with this limiting case, with  $\{\phi\} \in \mathbb{R}^{3N\times 3N}$  obtained as a result of some sort of a coordinate transformation. Subsequently  $\{\phi\}$  is divided into m coordinates  $\{\phi_v\} \in \mathbb{R}^{3N-m} \times {}^{3N-m}$  that will be constrained and  $\{\phi_f\} \in \mathbb{R}^{N\times N}$  that will remain free, namely the reduced coordinate set.

In the fully flexible limiting case, the dynamics of a system can be described by Hamilton's equations,

$$\frac{d\left[\phi_{f},\phi_{v}\right]}{dt} = \frac{\partial \mathcal{H}}{\partial\left[p_{f},p_{v}\right]}$$

$$\frac{d\left[p_{f},p_{v}\right]}{dt} = -\frac{\partial \mathcal{H}}{\partial\left[\phi_{f},\phi_{v}\right]}$$
(1)

where the Hamiltonian includes the kinetic and potential energy terms:

$$H(\phi_f, \phi_v, p_f, p_v) = [p_f, p_v]^T M_{tot}(\phi_f, p_v) [p_f, p_v] + U((\phi_f, \phi_v)$$
(2)

and  $M_{\rm tot}$  is the mass matrix tensor defined by  $M_{\rm tot} = JMJ$  where J is the Jacobian  $\{dx_{ii}/d\varphi_f\}$  of the transformation and M is the diagonal Cartesian mass matrix. The Boltzmann probability of a microstate  $(\varphi_f, \varphi_v, p_f, p_v)$  in canonical ensemble:

$$\rho(\phi_f, \phi_v, p_f, p_v) \propto e^{-\beta H \left(\phi_f, \phi_v, p_f, p_v\right)}.$$
(3)

<span id="page-1-0"></span>Most thermodynamic quantities of interest are based on configurations, not momenta. The marginal probability of a configuration can be determined by a Gaussian integral over the momenta to yield:

$$\rho(\phi_f, \phi_v) \propto |\mathcal{M}_{\text{tot}}(\phi)|^{-\frac{1}{2}} e^{-\beta U(\phi_f, \phi_v)}$$
(4)

In a similar fashion, a system represented in reduced coordinates will evolve according to the following Hamiltonian:

$$H(\phi_f, p_f) = p_f^T \mathcal{M}(\phi_f) p_f + U(\phi_f)$$
(5)

<span id="page-1-1"></span>and subsequently, marginal distribution of a configuration will be:

$$\rho(\phi_f) \propto |\mathcal{M}(\phi_f)|^{-\frac{1}{2}} e^{-\beta U(\phi_f)} \tag{6}$$

There is a mismatch between Eq. (4) and Eq. (6). We would like to conduct constrained simulations in which the unnormalized density matches the unnormalized density of unconstrained dynamics. To do so, we may use a correcting potential known as the Fixman potential [30,31],

$$U'(\phi_f) = kT \ln \left( \frac{|\mathcal{M}_{\text{tot}}(\phi)|}{|\mathcal{M}(\phi_f)|} \right)^{\frac{1}{2}}$$
(7)

There are at least two ways to use the Fixman potential such that the marginal probability of flexible coordinates  $\phi_f$  in the constrained and fully flexible regimes match. When used for Hamiltonian Monte Carlo [25], the Fixman potential can be included in the Metropolis-Hastings acceptance-rejection criterion [32]. If the derivative of the Fixman

potential, the Fixman torque [\[33](#page-10-16)], is computed, it can be used during molecular dynamics.

With the addition of Fixman correcting potential, the Hamiltonian becomes:

$$H(\phi_f, p_f) = p_f^T \mathcal{M}(\phi_f) p_f + U(\phi_f) + U'(\phi_f)$$
(8)

Compared to fully flexible dynamics, generalized coordinate dynamics is more computationally complex. To solve for acceleration inverting the mass matrix tensor is needed; and calculating the Fixman torque [\[34](#page-10-17)] is required in order for the MD trajectories to sample from the Boltzmann distribution. Directly inverting the mass matrix tensor has a complexity of O(N<sup>3</sup> ). However, efficient algorithms with linear O (N) complexity have been developed by exploiting the special structure of the tensor.

One of the O(N) algorithms is the spatial operator algebra (SOA) developed by Rodriguez and Kreutz in 1988 [[35\]](#page-10-18). SOA is a multibody dynamics formulation for internal coordinates initially used to describe robotic structure motion. SOA relies on the equivalence between the Kalman filtering / Bryson-Frazier smoothing equations used in signal processing [[36,](#page-10-19)[37\]](#page-10-20) and the dynamics of kinematic chains [\[38](#page-10-21)]. The filtering/smoothing equations take a time series of measurements as input and provides a smoothed version of their estimates as an output. In the dynamics of kinematic chains, spatial forces as a function of the link (body) are taken as input and joint torques as output. SOA uses linear operators which act on velocities, accelerations and forces of a linked multibody chain. The operators are associated with recursive algorithms that pass through the chain back and forth without requesting any matrix operations such as multiplication or inversion [[35](#page-10-18)[,39](#page-10-22)[,40](#page-10-23)]. Beside the kinematics and dynamics solutions, the SOA formulation allows for a straightforward derivation of mass matrix tensor-related quantities needed in the statistical mechanics of constrained systems: the mass matrix tensor determinant, its gradient, logarithm, and square root.

# 2. Methods

#### 2.1. Robosample architecture

Robosample further builds on the Simbody and the Molmodel APIs. Simbody API [[41\]](#page-10-24) was designed for general and accurate mechanical engineering calculations with an emphasis on biomedical applications. Simbody implements the SOA formulation and offers a range of integrators to solve the Hamilton equations. To facilitate sampling from the Boltzmann distribution, we added libraries for determinant, square root, and gradient of the mass matrix tensor computation. Molmodel [[28\]](#page-10-11) is a library designed for simulating chemical objects. It represents a chemical object as a tree and maps Simbody serial multibody subgraphs - which are dependent on the defined rigid bodies - onto the chemical tree [\(Fig. 2](#page-3-0)).

Robosample calls Simbody for trajectory integration while energies and forces are evaluated in Molmodel and can be computed via graphical processing units through the OpenMM API [[42\]](#page-10-25).

The program workflow is represented in [Fig. 1.](#page-3-1)

#### 2.1.1. Robosample implementation of Gibbs sampling

Robosample is designed to use constrained molecular dynamics to sample from the Boltzmann distribution following two conditions outlined by [[25\]](#page-10-8).

First, each sample is drawn from the conditional probability of the Boltzmann distribution. Velocities are initialized in accordance to a generalized equipartition theorem for generalized coordinates [\[43](#page-10-26)]. Trajectories are propagated deterministically, and the final configuration of each trajectory is accepted or rejected based on the Metropolis-Hastings criterion.

The second condition is that every DOF can be modified by some

Monte Carlo move. In Robosample, this condition is satisfied by alternating different blocks of coordinates, e.g. torsions, that are allowed to be flexible in a given Gibbs sampling move. Each block of coordinates specifies a "world", a defined set of rigid bodies and joints. Each world is implemented as a different Simbody multibody chain mapped onto the same Molmodel molecular graph [\(Fig. 2](#page-3-0)).

At the beginning of each Gibbs sampling move, Cartesian coordinates are used to update the generalized coordinates in the relevant world. An easy way to ensure that every DOF can be accessed is to include a world where the system is fully flexible.

#### 2.1.2. Robosample parameters setup

Simulating multiple worlds leads to a combinatorial increase in the number of parameters. In addition to selecting the number of worlds, for each world one needs to choose rigid body specifications, the types of joints, the number of MD steps in each HMC proposal, and the MD time step.

Although the development and testing of worlds that lead to efficient sampling of biomolecular systems is an unexplored area some general principles may apply. For instance, in proteins it could prove helpful to define rigid bodies based on secondary structure or protein domains. As Jensen et al. [[44\]](#page-10-27) suggested optimal/informative variables should be sampled more often, as for example torsions in contrast to bond lengths or angles due to that torsions play a larger role in distinguishing molecular conformations. Finally, it is known that including at least one fully flexible world per cycle is an easy way to ensure ergodicity.

HMC proposal parameters may be optimized based on the autocorrelation time and acceptance rate [\[45](#page-10-28)] [\[46](#page-10-29)]. HMC uses MD simulations to generate candidate configurations that are accepted or rejected according to the Metropolis criterion. For these MD trajectories, the time step (ε) and trajectory length (L), or number of integration steps, must be specified. If the integration time (T = L\*ε) is too short, then the system undergoes a random walk in configuration space and has a long autocorrelation time. Increasing the integration time reduces the autocorrelation time of the Markov chain. However, if the integration time is too long, then the effective sample size as a function of L is reduced [\[47](#page-10-30)]. The acceptance rate is the major consideration in selecting the time step. A large time step ε will reduce the number of steps required to attain a desired integration time T. However, a large time step comes at the cost of larger integrator error and a reduced acceptance rate. On the other hand, a small time step minimizes integrator error but increases the number of steps required to attain a desired integration time T. The acceptance rate that optimally balances these considerations is 0.651 [[46,](#page-10-29)[48\]](#page-10-31).

In summary, our workflow for optimizing parameters for Robosample consists of the following three steps:

- a. run short single world trial simulations of different Ts with ε = 1 fs;
- b. pick the maximum the sample-per-transition over T;
- c. given the chosen T (under b.) increase ε in mixed trial simulations until an optimal acceptance rate of ~0.651 is achieved.

To help keep track of all Robosample parameters, we have developed a GUI that assists users in setting up the simulations. The GUI consists of two main windows: one for the general simulation parameters such as the temperature or the number of cycles and one for every world specific parameters such as the number of molecular dynamics steps or integration time step. The second window is displayed in a tabular fashion so the user can have an overview of all the worlds ([Fig. 4\)](#page-4-0).

# 2.1.3. Spherical and cylindrical joints introduced in Robosample

Historically, the most attractive choice for constrained dynamics was torsional dynamics. Torsional dynamics is an intuitive choice due to the low frequency of torsional motions and their ability to travel

<span id="page-3-1"></span>![](_page_3_Figure_2.jpeg)

Fig. 1. Sequential UML diagram for main loop in Robosample that implements Gibbs sampling. The red rectangle is the main speed bottleneck in Robosample where the non-bonded energy terms must be evaluated in Cartesian coordinates. B. Module communication in Robosample. OpenMM, represented in orange has GPU capability. (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

<span id="page-3-0"></span>![](_page_3_Picture_4.jpeg)

Fig. 2. Schematic representation of two possible Simbody kinematic trees. B and C are mapped onto the chemical graph A. Each "world" has different rigid body definitions and joint types. B includes cylindrical and C includes spherical joints. Rigid bodies are colored grey, joints are colored red. (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

large distances in conformational space. Unfortunately, torsions are strongly coupled to other DOFs [[49\]](#page-10-32). To overcome this issue, the developers of GneimoSim added angle bend mobilities to molecular dynamics [[24\]](#page-10-7). Along the same vein, we have introduced, among others, a spherical joint in Robosample. A spherical joint has three rotational DOFs that are expressed in Simbody as Euler angles or quaternions. We have also introduced a cylindrical joint with two DOFs: a translation and a rotation around the translational axis ([Fig. 5](#page-4-1)). These two joints were implemented in Molmodel based on the Simbody Ball and Cylinder mobilized body types.

In Robosample a given 'child' rigid body is attached to its 'parent' rigid body through a 'mobilizer', i.e. a joint. Mathematically this scaffold is described by four reference frames - one for each of the two rigid bodies, and two defining the mobilizer/joint: one fixed on the parent and the second, mobile, on the child (see [Fig. 2,](#page-3-0) in ref. [\[50](#page-10-33)]).

In modelling molecules a user would naturally assume that the two mobilizer frames are placed at given atom centers and oriented with one of their axis along bonds - but in general robotics terms this is not a necessity, as the two joint-frames can be placed in principle anywhere within the parent and child rigid bodies. Preliminary experiments performed during Robosample optimization clearly indicate that placing the joint frames randomly, not in atoms centers result in poorer transition rates as shown in Supplementary Table S3. Consequently, Robosample places the mobilizer frames along bonds by default, letting the user specify the atom to be treated as fixed and the one that is treated as mobile.

2.2. Experimental setups for testing Robosample: defining the simulation parameters

To assess Robosample efficacy, simulations were performed herein on two systems: alanine dipeptide and a medium sized molecular model consisting of the first glycan of hepatitis C virus protein E2 (E2N1) attached to the E2 derived [−12 + 6] peptide around its Asn 17 site. The glycan likely offers protection or delays the recognition of host anti-E2 antibodies [\[51](#page-10-34)].

#### 2.2.1. Alanine dipeptide model and simulation parameters

Input files of alanine dipeptide were generated with AmberTools16 [[53\]](#page-10-35) using the ff14SB force field. Next two dynamic models were set in Robosample: (a) the first with 7 mobile joints corresponding to all bonds excepting the 2 terminal ones, and (b) the second with only 2 joints corresponding to the N-Cα and Cα-C bonds. Experiments performed with the first model were called 'all-bond dynamics', labeled 'all-' and the later Ramachandran dynamics, labeled "Rama-". Each of these two dynamics models was set for three types of moves: torsional (TD), torsion-bond cylindrical (Cyl), and torsion-angle spherical (Ball) resulting in a total of 6 worlds/regimens that were mixed with fully flexible worlds to ensure ergodicity (as in [Fig. 3\)](#page-4-2). In addition, separate fully flexible runs were performed for reference.

For each of the above 6 worlds, HMC parameters were tuned by first running a series of trial simulations with increasing HMC trajectory lengths and a short, fixed step of 2 fs (1 fs for the flexible world). For each HMC trajectory length - the autocorrelation time of the end-to-end distance were calculated and used to select the HMC trajectory length

<span id="page-4-2"></span>Fig. 3. Schematic representation of Robosample operation flow. Simulations are based on a sequence of alternating worlds, each with its own set of parameters including the number of MD steps and MD time step.

that maximizes the effective sample size per picosecond that was further used in all subsequent production runs. For the optimal HMC integration time, the time step was then increase until an acceptance rate near 0.651 was obtained [\(Fig. 6\)](#page-5-0).

We compared the efficiency of four types of simulations: 1) fully flexible (fully-Flex); 2) fully flexible mixed with torsional moves (mixed-TD); 3) fully flexible mixed with torsion-bond cylindrical moves (mixed-CYL); and 4) fully flexible mixed with torsion-angle spherical moves (mixed-BALL). These were performed in three instances: (a) allbond dynamics with unoptimized parameters, (b) all-bond dynamics with optimized parameters and (c) optimized Ramachandran dynamics - resulting in a total of 11 experiments as shown in [Table 1.](#page-8-0) In each type of simulation, the sequence of worlds in each cycle ([Fig. 3\)](#page-4-2) consisted of 1 fully flexible world followed by 10 CDHMC worlds (TD, Cyl, Ball). Each experiment was repeated 3 times.

#### 2.2.2. The E2N1 model and its parameters

The E2N1 glycan was attached to its corresponding peptide [−12 to +6] and the peptide was kept fixed in the conformation known to be recognized by the antibody [\[52](#page-10-36)]. Initial coordinates of the polypeptide

<span id="page-4-1"></span>![](_page_4_Picture_8.jpeg)

Fig. 5. Spherical and cylinder joint schematic representation. The mobile body is colored in blue and fixed frame is left blank. Although there are no rules on where to place the fixed and mobile frames inside the connected bodies, Robosample places the cylinder between two atoms at the default bond length. The spherical joint is placed by default at the beginning of a chemical bond and subsequently the bond acts like a rigid rod that rotates around a fixed atom. (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

were generated by building a HCV-E2 homology model with Modeller v9.21 [[53\]](#page-10-35) and the glycan structure was generated with Carbohydrate Builder v2.1 [[http://glycam.org\]](http://glycam.org). Simulation input files were prepared

<span id="page-4-0"></span>![](_page_4_Figure_11.jpeg)

| 🔞 🖯 🕕 Worlds Ove      | erview                            |                      |                      |  |
|-----------------------|-----------------------------------|----------------------|----------------------|--|
| Experiment 0          |                                   |                      |                      |  |
| World Type:           | World Type: IC                    |                      | TD CT                |  |
| Run Type:             | Normal                            | Normal               | Normal               |  |
| Flex Type:            | From File                         | From File            | From File            |  |
| Flex File:            | <pre><param not="" set=""/></pre> | g_a_GUI/peptide.flex | g_a_GUI/peptide.flex |  |
| Rigid Bodies File:    | <pre><param not="" set=""/></pre> | ing_a_GUI/peptide.rb | ing_a_GUI/peptide.rb |  |
| Reblock Frequency:    | 0                                 | 0                    | 0                    |  |
| Roots:                | 0                                 | 0                    | 0                    |  |
| Sampler:              | HMC                               | HMC                  | HMC                  |  |
| Timestep:             | 0,002                             | 0.05                 | 0.05                 |  |
| MD Steps:             | 0                                 | 100                  | 100                  |  |
| Boost MD Steps:       | 0                                 | 10                   | 10                   |  |
| Samples Per Round:    | 0                                 | 10                   | 10                   |  |
| Seed:                 | 0                                 | 1                    | 2                    |  |
| Use Fixman Potential: | True                              | True                 | False                |  |
| Use Fixman Torque:    | True                              | True                 | False                |  |
| Root Mobility:        | Root Mobility: Weld               |                      | Weld                 |  |
| Experiment 1          |                                   | IC                   |                      |  |
| World Type:           | World Type: IC                    |                      | IC                   |  |
| Run Type:             | Normal                            | Normal               | Normal               |  |
| Flex Type:            | From File                         | From File            | From File            |  |
| Flex File:            | g_a_GUI/peptide.flex              | g_a_GUI/peptide.flex | g_a_GUI/peptide.flex |  |
| Rigid Bodies File:    | ing_a_GUI/peptide.rb              | ing_a_GUI/peptide.rb | ing_a_GUI/peptide.rb |  |
| Reblock Frequency:    | 0                                 | 0                    | 0                    |  |
| Roots:                | Roots: 0                          |                      | 0                    |  |
| Sampler:              | Sampler: HMC                      |                      | HMC                  |  |
| Timestep:             | 0,002                             | 0.002                | 0,002                |  |
| MD Steps:             | MD Steps: 50                      |                      | 50                   |  |
| Boost MD Steps:       | Boost MD Steps: 10                |                      | 25                   |  |
| Samples Per Round:    | 10                                | 10                   | 25                   |  |
| Seed:                 | Seed: 0                           |                      | 2                    |  |
| Use Fixman Potential: | Jse Fixman Potential: False       |                      | False                |  |
| Use Fixman Torque:    | Use Fixman Torque: True           |                      | True                 |  |
| Root Mobility:        | Cartesian                         | Free                 | Free                 |  |

Fig. 4. Robosample GUI assists the user with the simulation parameters. The window on the left is used to set overall simulation parameters while the one on the right gives an overview of the different worlds' specific parameters.

<span id="page-5-0"></span>![](_page_5_Figure_2.jpeg)

![](_page_5_Figure_3.jpeg)

![](_page_5_Figure_4.jpeg)

![](_page_5_Figure_5.jpeg)

Fig. 6. Curves used for finding the optimal HMC trajectory length for Ramachandran dynamics: A) fully flexible, B) mixed-RamaTD, C) mixed-RamaBall and D) mixed-RamaCyl. The blue line traces the independent samples per move. The orange one represents the number of samples per integration time and it was used to find the optimal HMC trajectory length. (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

with AmberTools16 [\[54](#page-10-37)] using the ff14SB force field for polypeptides and the GLYCAM06 force field for the polysaccharide moiety [[55,](#page-10-38)[56](#page-10-39)]. The structure was minimized with the same force field using OpenMM v7.4 until it reached 1 kJ/mol.

For the E2N1 glycan, CDHMC simulations were performed at 300 K with four worlds: −1) fully flexible; −2) with rigid monosaccharide rings; −3) with mixed spherical and cylindrical moves for the polysaccharide structure; and − 4) torsional mobilities at the Asn anchor. HMC dynamics trajectory lengths and time steps were 108 × 1 fs, 25 × 8 fs, 158 × 0.5 fs and 11 × 20 fs respectively.

In the second and third world, the monosaccharide structures hydroxyl, carbonyl, secondary amino methyl, and methylene groups were kept rigid. Spherical joints were placed between the rigid polysaccharide groups except for the groups that had only two atoms where cylindrical joints used. Each simulation was repeated three times.

On this model we performed four types of simulations: two types of MD - 1a) fully flexible MD and 1b) rigid body MD (RB-MD), and two types of HMC - 2a) fully flexible and 2b) multiple rigid body (RB) worlds, as described in [Table 2](#page-8-1). Rigid bodies covered monosaccharides or their lateral groups.

Trajectories were analyzed with VMD version 1.9.4 [[57\]](#page-10-40), NumPy 1.11.0 [[58\]](#page-10-41), and SciPy 0.17 [\[59](#page-10-42)]. Tridimensional structures were visualized with VMD and charts were generated using Matplotlib 1.5.1 [[60\]](#page-10-43) and Excel from Microsoft Office 365 [[61\]](#page-10-44).

The extent of the conformational space exploration was assessed for both MD and HMC simulations using all-vs-all RMSD analysis and cluster count analysis. For cluster count analysis, hierarchical agglomerative clustering was performed using Ward's method [\[62](#page-10-45)] on a series of RMSD and flattened at a threshold of 2 Å.

To directly compare the computational gain coming from rigid body simulations, time was measured in computer runtime seconds. All simulations were performed on an Intel(R) Core (TM) i9-9900K CPU (3.60GHz), 2 X GeForce RTX 2080Ti GPUs and 64 GB RAM with Ubuntu 18.04 OS. More specifically, E2N1 flexible MD simulations (747 degrees of freedom) take Robosample ~0.022 sec / step while RB-MD ones take ~0.008 sec / step (105 DOFs). CDHMC simulation benchmarks indicated ~0.007 sec / step for world 2 (74 DOFs), identical with RB-MD for world 3 and ~0.001 sec / step for world 4 (4 DOFs).

### 3. Results

#### 3.1. Mixed simulations accurately reconstruct the free energy surface of alanine dipeptide

Although alanine dipeptide - N-acetylalanine-N-methylamide - is small compared to other biomolecules, its potential energy surface (PES) is highly frustrated and presents a challenge to molecular sampling methods. The projection of its configuration space onto φ and ψ angles is able to retain most of the PES maxima and minima while smoothing the terrain. Therefore, representing the potential of mean force in the domain of φ and ψ angles is widely used to assess the quality and effectiveness of sampling methods [\[63](#page-10-46)].

To validate that the new software reproduces previous results, we ran a set of simulations with the same parameters as in Spiridon and

<span id="page-6-0"></span>![](_page_6_Figure_2.jpeg)

Fig. 7. The free energy landscape of alanine dipeptide, based on histograms of four types of HMC simulations with optimized parameters: (A) fully flexible, (B) mixed-RamaTD, (C) mixed-RamaCYL, (D) mixed-RamaBALL. Simulations were run for 120,000 Monte Carlo moves. Histograms were generated using a 15-degree spacing. Free energies are displayed in units of kJ/mol. Unsampled regions are shown in red. (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

Minh, 2017. Indeed, the results for the fully flexible and for the mixed-TD simulations closely resembled our previous results (Supplementary Table S1 – flexible and mixed-TD).

Based on the optimization, we used the following trajectory lengths

and time steps for all-bond dynamics: fully flexible moves used 108 steps of 1.87 fs (202.0 fs), torsional moves used 25 steps of 17.1 fs (427.5 fs), ball moves used 158 steps of 2.1 fs (331.8 fs), and cylindrical moves used 47 steps of 6.1 fs (286.7 fs). For the Ramachandran

<span id="page-6-1"></span>![](_page_6_Picture_7.jpeg)

Fig. 8. (A) Schematic representation of E2N1: N-AcGlc is in blue squares, mannose in green circles, galactose in yellow circles, sialic acid in purple diamonds; (B) E2N1 molecular model; (C) Superposition of first 400 ps of fully flexible MD simulation; (D) The C equivalent in the rigid body MD regimen. (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)

<span id="page-7-0"></span>![](_page_7_Figure_2.jpeg)

Fig. 9. RMSD timeseries for flexible and rigid body dynamics. Upper panel represents one of the flexible HMC simulation; the lower panel represents a CDHMC simulation over 4.8 ns and 52.1 ns respectively.

dynamics, torsional moves used 11 steps of 44.73 fs (492.0 fs), ball moves used 51 steps of 8.957 fs (456.8 fs), and cylindrical moves used 49 steps of 12.386 fs (606.9 fs). It appears that ball and cylinder dynamics can use a larger time step than flexible dynamics, but a much smaller time step than torsional dynamics. Our data also suggest that larger time steps are more feasible with the larger rigid bodies of Ramachandran dynamics.

All simulations were able to recover the free energy surface of alanine dipeptide ([Fig. 7](#page-6-0) and Supplementary Figs. S2 and S3). In our previous work on CDHMC, we showed that the alanine dipeptide free energy surface can be recovered by combining flexible HMC with torsional moves [\[25](#page-10-8)]. Here we show that CDHMC simulations based on the new cylindrical and spherical joints can also accurately recover the surface, irrespective of the size of the bodies. The free energy surfaces are not only consistent among the different worlds, but the current simulations at 625 K have basins in the same positions as previously published results at the same temperature [[25\]](#page-10-8) or at 800 K [\[24](#page-10-7)] ([Table 3](#page-8-2)).

Optimization of HMC increased the efficiency of sampling a rare event in all-bond simulations (Suppl. Table 1 vs Suppl. Table 2). For all types of simulation except all-bond Ball simulations, simulations with optimized parameters have faster transitions to and from the αL basin or crossing the phi = 0 barrier. The significant improvement demonstrates the importance of HMC parameter optimization. The all-bond Ball optimization procedure, although leading to slightly better results, had only a limited improvement (Suppl. Fig. S1). It appears that the Ball joint is less subject to diffusive motion than the others.

The mean first passage time between free energy basins provides insight into the relative efficiency of different types of alanine dipeptide simulations [\(Table 4](#page-9-7) and Suppl. Tables 1 and 2). Fully flexible simulations have the slowest transitions between C5, PPII, C7eq, and aL. In the all-bonds dynamics regime, Ball and Cylinder joints perform significantly better than the flexible Cartesian joints, while the torsional joint performs best (Suppl Table 2). All-bonds dynamics are all outperformed by different types of Ramachandran dynamics ([Table 4](#page-9-7)). Compared to fully flexible dynamics, the MPFT across the φ = 0 barrier is decreased by an order of magnitude. The MFPT between other basins within the two sides of the barrier is decreased by about a half. The efficiency of joints increases in the following order: Torsion, Cylinder and Ball.

#### 3.2. Using multiple worlds improve simulation efficiency

Glycans are medium- to large- sized polysaccharides attached to protein structures and are ubiquitously present along the evolutionary tree. The fact that all cells spend energy to build them suggest that they have critical functions. These may include protein physical protection, protein folding, specific recognition, and solubility [\[64](#page-10-47)]. However as they are flexible glycans are rarely present in solved structures [[4](#page-9-8),[9](#page-9-9)]. Hence flexibility makes them of special interest for conformational space analysis by molecular simulation. Moreover, the presence of multiple cyclic structures in glycans makes them a good target for rigid body simulation. Traditional simulations in Cartesian space are hindered by high-frequency motion within monosaccharides, which are mostly in the chair conformation.

To illustrate the effectiveness of using multiple worlds with different types of joints - both molecular dynamics (MD) and Hamiltonian Monte Carlo (HMC) simulations were performed using fully flexible vs other types of mixed movements were performed with Robosample on a glycan structure attached to a polypeptide [\(Fig. 8](#page-6-1)).

<span id="page-7-1"></span>![](_page_7_Figure_12.jpeg)

![](_page_7_Figure_13.jpeg)

Fig. 10. Cumulative cluster count analysis for the six MD simulation in the left panel and six HMC simulations in the right panel express in CPU time. Flexible simulations are represented with black lines while CDHMC simulations with red lines.

<span id="page-8-3"></span>![](_page_8_Figure_2.jpeg)

![](_page_8_Figure_3.jpeg)

Fig. 11. All vs all RMSD analysis for E2G1 simulations. The left panel compares flexible-MD with RB-MD while the right panel compares flexible-HMC with RB-HMC.

<span id="page-8-0"></span>Table 1 Simulations performed on alanine dipeptide.

| Dynamics                        | All-                                                        |                                                              | Rama                                                       |
|---------------------------------|-------------------------------------------------------------|--------------------------------------------------------------|------------------------------------------------------------|
| Optimization<br>Simulation type | NO (a)<br>fully-Flex<br>mixed-TD<br>mixed-Cyl<br>mixed-Ball | YES (b)<br>fully-Flex<br>mixed-TD<br>mixed-Cyl<br>mixed-Ball | YES (c)<br>mixed-RamaTD<br>mixed-RamaCyl<br>mixed-RamaBall |

<span id="page-8-1"></span>Table 2 Types of simulations of E2N1 glycan performed with Robosample.

|          | MD          | HMC          |
|----------|-------------|--------------|
| Flexible | flexible MD | flexible HMC |
| RB       | RB-MD       | CDHMC        |

<span id="page-8-2"></span>Table 3 Boundaries of the Four Major Conformational Regions in Alanine.

|      | ϕmin  | ϕmax | ψmax  | ψmax |
|------|-------|------|-------|------|
| C5   | −180° | −95° | 105°  | 180° |
| PPII | −96°  | −45° | 105°  | 180° |
| C7eq | −96°  | −45° | −25°  | 104° |
| αL   | 35°   | 85°  | −180° | 25°  |

#### 3.3. Robosample increases simulation efficiency

In the same amount of computer wall clock time, the mean value and fluctuations of the RMSD are significantly larger in representative mixed-world versus fully flexible HMC simulations ([Fig. 9](#page-7-0)).

The improved configuration space sampling of rigid-body dynamics is also corroborated by cluster count analysis [\(Fig. 10\)](#page-7-1). Compared to fully flexible simulations, HMC simulations that involve rigid bodies access a larger number of clusters in the same amount of computer run time. The contrast between rigid and fully flexible simulation is only evident with HMC, which involves multiple worlds with rigid bodies, opposed to MD, where there is only one rigid-body regime.

The full RMSD matrix of all MD and HMC simulations of E2N1 also demonstrates the improved sampling due to multiple rigid worlds ([Fig. 11](#page-8-3)). The multiple-world HMC simulations routinely access configurations in which the RMSD compared to the initial structure reaches up to 30 Å. Moreover, the structures accessed in the independent simulations are distinct from each other. In contrast, these large-RMSD structures are only evident in one of the MD simulations with one rigid world. All the fully flexible MD and HMC simulations maintain an RMSD with the initial structure around 15 Å or less, and the three independent simulations assess a similar part of configuration space.

#### 4. Discussion and future directions

Results presented herein indicate that Robosample is a good environment to perform constrained molecular simulations and achieve ergodicity through the rigorous framework of Gibbs sampling.

In addition to torsional and angle/torsion mobilities commonly used in constrained simulations [\[23](#page-10-6)], Robosample introduces other types of robotic mobility combinations by exploiting the mechanical joints implemented in Simbody. This allows examining arbitrary DOF couplings such as the weak bond/torsion couplings in cylindrical joints.

Using a simple HMC parameter optimization procedure and different joint types for alanine dipeptide, we obtained enhanced transitions to and from the αL basin in the free energy surface. With all-bond dynamics, all joint types outperformed flexible dynamics, with torsional dynamics providing the greatest acceleration. Ramachandran dynamics, in which only the φ and ψ angles are mobile, was found to be significantly more efficient than fully-flexible or all-bond dynamics with Ball dynamics resulting in the shortest MFPT for the rarest transition.

Besides alanine dipeptide, the benefits of rigid-body simulation are also evident in a large system, the E2N1 glycan attached to its (−12 to +6) peptide. Results on this system confirm the ability of alternating constraints to explore larger volumes of the configuration space than those reached in fully flexible simulations. The effectiveness of using multiple constrained worlds compared to fully flexible MD is clearly shown by a cluster count analysis [\(Fig. 10\)](#page-7-1) and RMSD matrix [\(Fig. 11](#page-8-3)).

In addition to simulations performed to estimate ensemble averages, Robosample may be useful for structure prediction. In homology modelling, the starting structure is usually refined based on a multiscale approach that starts from sidechains and gradually goes through loops until hopefully attaining the global minimum. Rigid body simulations can provide another approach to this classic multi-scaling procedure. By switching different types of degrees of freedom through the implemented joints, one can first optimize the side chains using cylindrical moves and then perform loop modelling using spherical moves. In this way, a large volume of configuration space may be sampled in the search for the deepest minimum.

Although rigid body simulation has many advantages, it should be mentioned that they have a disadvantage in the inability to accurately simulate time-dependent quantities. MD simulations are often used to study time-dependent quantities such as correlated motion, diffusion coefficients, or binding kinetics. Once a coarsening is used, time-dependent functions start losing accuracy due to entropic effects encoded in the lost degrees of freedom. For example, if monosaccharides are

<span id="page-9-7"></span>Table 4 Mean first passage times of Ramachandran dynamics simulations expressed in molecular dynamics steps, divided by 200 (to be comparable with our last paper which measured the MFPT in cycles of 200 MD steps). Within each matrix, row i and column j indicate the MFPT for moving from state i to state j. The quantities in parentheses are standard deviations over three independent simulations.

|                | Fully flexible |               |               |                |
|----------------|----------------|---------------|---------------|----------------|
| from \ to      | C5             | PP2           | C7eq          | αL             |
| C5             | 3.2 ± 0.2      | 9.1 ± 0.5     | 14.3 ± 1.7    | 5380.0 ± 285.0 |
| PP2            | 9.6 ± 0.9      | 5.7 ± 0.4     | 10.9 ± 1.8    | 5380.0 ± 285.0 |
| C7eq           | 16.5 ± 1.1     | 12.6 ± 0.6    | 2.8 ± 0.2     | 5370.0 ± 285.0 |
| αL             | 502.0 ± 195.0  | 499.0 ± 194.0 | 491.0 ± 194.0 | 14.1 ± 4.4     |
| Mixed-RamaTD   |                |               |               |                |
| C5             | 1.8 ± 0.03     | 4.0 ± 0.1     | 2.8 ± 0.1     | 631.0 ± 92.2   |
| PP2            | 2.5 ± 0.1      | 3.2 ± 0.1     | 2.6 ± 0.1     | 631.0 ± 92.2   |
| C7eq           | 3.2 ± 0.2      | 4.7 ± 0.2     | 1.5 ± 0.02    | 630.0 ± 92.2   |
| αL             | 50.9 ± 1.1     | 52.4 ± 1.2    | 49.4 ± 1.1    | 8.3 ± 1.3      |
| Mixed-RamaCyl  |                |               |               |                |
| C5             | 2.3 ± 0.03     | 5.2 ± 0.1     | 3.5 ± 0.02    | 615.0 ± 33.6   |
| PP2            | 3.2 ± 0.1      | 4.1 ± 0.1     | 3.3 ± 0.005   | 615.0 ± 33.4   |
| C7eq           | 4.1 ± 0.1      | 5.9 ± 0.1     | 1.9 ± 0.01    | 614.0 ± 33.6   |
| αL             | 46.4 ± 4.5     | 48.4 ± 4.6    | 44.7 ± 4.7    | 12.0 ± 1.4     |
| Mixed-RamaBall |                |               |               |                |
| from \ to      | C5             | PP2           | C7eq          | αL             |
| C5             | 2.4 ± 0.01     | 5.2 ± 0.01    | 3.1 ± 0.1     | 518.0 ± 18.6   |
| PP2            | 3.1 ± 0.01     | 4.1 ± 0.04    | 3.0 ± 0.1     | 518.0 ± 18.5   |
| C7eq           | 3.7 ± 0.01     | 5.6 ± 0.02    | 1.9 ± 0.02    | 518.0 ± 18.4   |
| αL             | 37.0 ± 1.4     | 39.1 ± 1.5    | 35.8 ± 1.6    | 12.8 ± 1.0     |
|                |                |               |               |                |

treated as rigid bodies, then chair – envelope - boat transitions will not be observed.

Further improvements in sampling may be obtained by optimizing the combination of various worlds. HMC moves based on flexible dynamics may be followed by multiple HMC moves (opposed to only one) based on constrained dynamics. Spherical and cylindrical moves might be especially complementary, as they both drive torsions while alternatively allowing coupling to bond lengths and angles.

Sampling in Robosample may be further improved by combination with other enhanced sampling techniques such as replica exchange [[65](#page-10-48)[,66](#page-10-49)]. Another possible direction for the improvement of Robosample is the implementation of additional joints that are already available in Simbody, such as the Slider.

# Declaration of Competing Interest

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

# Acknowledgements

This work was supported by UEFISCDI grant PN-III-P1 - 1.1-TE-2016-1852 and by the Romanian Academy Program 2 of the Institute of Biochemistry. The authors would like to thank Eliza Martin for involvement in software development, optimization and very insightful ideas and to Victor Gabriel Ungureanu for software development. AJP would like to thank Jeremy Smith for his instrumental contribution in setting the foundation of the Department of Bioinformatics and Structural Biochemistry in the Institute of Biochemistry of the Romanian Academy back in 1999.

# Appendix A. Supplementary data

Supplementary data to this article can be found online at [https://](https://doi.org/10.1016/j.bbagen.2020.129616) [doi.org/10.1016/j.bbagen.2020.129616](https://doi.org/10.1016/j.bbagen.2020.129616).

#### References

- <span id="page-9-0"></span>[1] K.M. Ruff, R.V. Pappu, A.S. Holehouse, Conformational preferences and phase behavior of intrinsically disordered low complexity sequences: insights from multiscale simulations, Curr. Opin. Struct. Biol. 56 (2019) 1–10, [https://doi.org/10.](https://doi.org/10.1016/j.sbi.2018.10.003) [1016/j.sbi.2018.10.003.](https://doi.org/10.1016/j.sbi.2018.10.003)
- [2] J. Habchi, P. Tompa, S. Longhi, V.N. Uversky, Introducing protein intrinsic disorder, Chem. Rev. 114 (2014) 6561–6588, [https://doi.org/10.1021/cr400514h.](https://doi.org/10.1021/cr400514h)
- <span id="page-9-1"></span>[3] M.R. Wormald, A.J. Petrescu, Y.-L. Pao, A. Glithero, T. Elliott, R.A. Dwek, Conformational studies of oligosaccharides and Glycopeptides: complementarity of NMR, X-ray crystallography, and molecular modelling, Chem. Rev. 102 (2002) 371–386, [https://doi.org/10.1021/cr990368i.](https://doi.org/10.1021/cr990368i)
- <span id="page-9-8"></span>[4] A.J. Petrescu, S.M. Petrescu, R.A. Dwek, M.R. Wormald, A statistical analysis of Nand O-glycan linkage conformations from crystallographic data, Glycobiology. 9 (1999) 343–352, [https://doi.org/10.1093/glycob/9.4.343.](https://doi.org/10.1093/glycob/9.4.343)
- <span id="page-9-2"></span>[5] Z. Peng, J. Yan, X. Fan, M.J. Mizianty, B. Xue, K. Wang, G. Hu, V.N. Uversky, L. Kurgan, Exceptionally abundant exceptions: comprehensive characterization of intrinsic disorder in all domains of life, Cell. Mol. Life Sci. 72 (2014) 137–151, [https://doi.org/10.1007/s00018-014-1661-9.](https://doi.org/10.1007/s00018-014-1661-9)
- <span id="page-9-3"></span>[6] R. Van Der Lee, M. Buljan, B. Lang, R.J. Weatheritt, G.W. Daughdrill, A.K. Dunker, M. Fuxreiter, J. Gough, J. Gsponer, D.T. Jones, P.M. Kim, R.W. Kriwacki, C.J. Oldfield, R.V. Pappu, P. Tompa, V.N. Uversky, P.E. Wright, M.M. Babu, Classification of intrinsically disordered regions and proteins, Chem. Rev. 114 (2014) 6589–6631, [https://doi.org/10.1021/cr400525m.](https://doi.org/10.1021/cr400525m)
- <span id="page-9-4"></span>[7] M.B. Marin, S. Ghenea, L.N. Spiridon, G.N. Chiritoiu, A.J. Petrescu, S.M. Petrescu, Tyrosinase degradation is prevented when EDEM1 lacks the intrinsically disordered region, PLoS One 7 (2012), [https://doi.org/10.1371/journal.pone.0042998.](https://doi.org/10.1371/journal.pone.0042998)
- <span id="page-9-5"></span>[8] A.J. Petrescu, T.D. Butters, G. Reinkensmeier, S. Petrescu, F.M. Platt, R.A. Dwek, M.R. Wormald, The solution NMR structure of glucosylated N-glycans involved in the early stages of glycoprotein biosynthesis and folding, EMBO J. 16 (1997) 4302–4310, [https://doi.org/10.1093/emboj/16.14.4302.](https://doi.org/10.1093/emboj/16.14.4302)
- <span id="page-9-9"></span>[9] A.J. Petrescu, A.L. Milac, S.M. Petrescu, R.A. Dwek, M.R. Wormald, Statistical analysis of the protein environment of N-glycosylation sites: implications for occupancy, structure, and folding, Glycobiology. 14 (2004) 103–114, [https://doi.org/](https://doi.org/10.1093/glycob/cwh008) [10.1093/glycob/cwh008.](https://doi.org/10.1093/glycob/cwh008)
- [10] C. Paduraru, L. Spiridon, W. Yuan, G. Bricard, X. Valencia, S.A. Porcelli, P.A. Illarionov, G.S. Besra, S.M. Petrescu, A.J. Petrescu, P. Cresswell, An N-linked glycan modulates the interaction between the CD1d heavy chain and β2-microglobulin, J. Biol. Chem. 281 (2006) 40369–40378, [https://doi.org/10.1074/jbc.](https://doi.org/10.1074/jbc.M608518200) [M608518200.](https://doi.org/10.1074/jbc.M608518200)
- <span id="page-9-6"></span>[11] J.J. Ziarek, D. Baptista, G. Wagner, Recent developments in solution nuclear magnetic resonance (NMR)-based molecular biology, J. Mol. Med. 96 (2018), [https://doi.org/10.1007/s00109-017-1560-2.](https://doi.org/10.1007/s00109-017-1560-2)
- [12] J. Trewhella, Small-angle scattering and 3D structure interpretation, Curr. Opin. Struct. Biol. 40 (2016) 1–7, [https://doi.org/10.1016/j.sbi.2016.05.003.](https://doi.org/10.1016/j.sbi.2016.05.003)
- [13] A.-J. Petrescu, P. Calmettes, D. Durand, V. Receveur, J.C. Smith, Change in backbone torsion angle distribution on protein folding, Protein Sci. 9 (2000) 1129–1136, [https://doi.org/10.1110/ps.9.6.1129.](https://doi.org/10.1110/ps.9.6.1129)
- [14] A.J. Petrescu, V. Receveur, P. Calmettes, D. Durand, M. Desmadril, B. Roux, Smallangle neutron scattering by a strongly denatured protein: analysis using random

- polymer theory, Biophys J. 72 (1997) 335–342, [https://doi.org/10.1016/S0006-](https://doi.org/10.1016/S0006-3495(97)78672-7) [3495\(97\)78672-7.](https://doi.org/10.1016/S0006-3495(97)78672-7)
- [15] W.F. van Gunsteren, M. Karplus, W.F. van Gunsteren, Effect of constraints on the dynamics of macromolecules, Macromolecules. 15 (1982) 1528–1544, [https://doi.](https://doi.org/10.1021/ma00234a015) [org/10.1021/ma00234a015.](https://doi.org/10.1021/ma00234a015)
- [16] A.-J. Peterscu, V. Receveur, P. Calmettes, D. Durand, J.C. Smith, Excluded volume in the configurational distribution of a strongly-denatured protein, Protein Sci. 7 (1998) 1396–1403, [https://doi.org/10.1002/pro.5560070616.](https://doi.org/10.1002/pro.5560070616)
- <span id="page-10-0"></span>[17] U.R. Shrestha, P. Juneja, Q. Zhang, V. Gurumoorthy, J.M. Borreguero, V. Urban, X. Cheng, S.V. Pingali, J.C. Smith, H.M. O'Neill, L. Petridis, Generation of the configurational ensemble of an intrinsically disordered protein from unbiased molecular dynamics simulation, Proc. Natl. Acad. Sci. 116 (2019) 20446–20452, [https://doi.org/10.1073/pnas.1907251116.](https://doi.org/10.1073/pnas.1907251116)
- <span id="page-10-1"></span>[18] J.P. Ryckaert, G. Ciccotti, H.J.C. Berendsen, Numerical integration of the cartesian equations of motion of a system with constraints: molecular dynamics of n-alkanes, J. Comput. Phys. 23 (1977) 327–341, [https://doi.org/10.1016/0021-9991\(77\)](https://doi.org/10.1016/0021-9991(77)90098-5) [90098-5.](https://doi.org/10.1016/0021-9991(77)90098-5)
- <span id="page-10-2"></span>[19] H.C. Andersen, Rattle: a "velocity" version of the shake algorithm for molecular dynamics calculations, J. Comput. Phys. 52 (1983) 24–34, [https://doi.org/10.](https://doi.org/10.1016/0021-9991(83)90014-1) [1016/0021-9991\(83\)90014-1.](https://doi.org/10.1016/0021-9991(83)90014-1)
- <span id="page-10-3"></span>[20] L.M. Rice, A.T. BrüNger, Torsion angle dynamics: reduced variable conformational sampling enhances crystallographic structure refinement, Proteins Struct. Funct. Bioinforma. 19 (1994) 277–290, [https://doi.org/10.1002/prot.340190403.](https://doi.org/10.1002/prot.340190403)
- <span id="page-10-4"></span>[21] E.G. Stein, L.M. Rice, A.T. Brünger, Torsion-angle molecular dynamics as a new efficient tool for NMR Structure calculation, J. Magn. Reson. 124 (1997) 154–164, [https://doi.org/10.1006/jmre.1996.1027.](https://doi.org/10.1006/jmre.1996.1027)
- <span id="page-10-5"></span>[22] G.S. Balaraman, I.-H. Park, A. Jain, N. Vaidehi, Folding of small proteins Using constrained molecular dynamics, J. Phys. Chem. B 115 (2011) 7588–7596, [https://](https://doi.org/10.1021/jp200414z) [doi.org/10.1021/jp200414z.](https://doi.org/10.1021/jp200414z)
- <span id="page-10-6"></span>[23] N. Vaidehi, A. Jain, Internal coordinate molecular dynamics: a foundation for multiscale dynamics, J. Phys. Chem. B 119 (2015) 1233–1242, [https://doi.org/10.](https://doi.org/10.1021/jp509136y) [1021/jp509136y.](https://doi.org/10.1021/jp509136y)
- <span id="page-10-7"></span>[24] S. Kandel, R. Salomon-Ferrer, A.B. Larsen, A. Jain, N. Vaidehi, Overcoming potential energy distortions in constrained internal coordinate molecular dynamics simulations, J. Chem. Phys. 144 (2016) 1–14, [https://doi.org/10.1063/1.4939532.](https://doi.org/10.1063/1.4939532)
- <span id="page-10-8"></span>[25] L. Spiridon, D.D.L. Minh, Hamiltonian Monte Carlo with constrained molecular dynamics as Gibbs sampling, J. Chem. Theory Comput. 13 (2017) 4649–4659, [https://doi.org/10.1021/acs.jctc.7b00570.](https://doi.org/10.1021/acs.jctc.7b00570)
- <span id="page-10-9"></span>[26] S. Geman, D. Geman, Stochastic relaxation, Gibbs distributions, and the Bayesian restoration of images, IEEE Trans. Pattern Anal. Mach. Intell. PAMI-6 (1984) 721–741, [https://doi.org/10.1109/TPAMI.1984.4767596.](https://doi.org/10.1109/TPAMI.1984.4767596)
- <span id="page-10-10"></span>[27] S. Duane, A.D. Kennedy, B.J. Pendleton, D. Roweth, Hybrid Monte Carlo, Phys. Lett. B. 195 (1987) 216–222, [https://doi.org/10.1016/0370-2693\(87\)91197-X.](https://doi.org/10.1016/0370-2693(87)91197-X)
- <span id="page-10-11"></span>[28] S.C. Flores, M.A. Sherman, C.M. Bruns, P. Eastman, R.B. Altman, Fast flexible modeling of RNA structure using internal coordinates, IEEE/ACM Trans. Comput. Biol. Bioinforma. 8 (2011) 1247–1257, [https://doi.org/10.1109/TCBB.2010.104.](https://doi.org/10.1109/TCBB.2010.104)
- <span id="page-10-12"></span>[29] D.D.L. Minh, Alchemical grid dock (AlGDock): binding free energy calculations between flexible ligands and rigid receptors, J. Comput. Chem. (2019), [https://doi.](https://doi.org/10.1002/jcc.26036) [org/10.1002/jcc.26036.](https://doi.org/10.1002/jcc.26036)
- <span id="page-10-13"></span>[30] M. Fixman, Classical statistical mechanics of constraints: a theorem and application to polymers, Proc. Natl. Acad. Sci. U. S. A. 71 (1974) 3050–3053, [https://doi.org/](https://doi.org/10.1073/pnas.71.8.3050) [10.1073/pnas.71.8.3050.](https://doi.org/10.1073/pnas.71.8.3050)
- <span id="page-10-14"></span>[31] M. Fixman, Simulation of polymer dynamics. I. General theory, J. Chem. Phys. 69 (1978) 1527–1537, [https://doi.org/10.1063/1.436725.](https://doi.org/10.1063/1.436725)
- <span id="page-10-15"></span>[32] [W.K. Hastings, Monte Carlo Sampling Methods Using Markov Chains and their](http://refhub.elsevier.com/S0304-4165(20)30128-8/rf0160) [Applications, \(1970\).](http://refhub.elsevier.com/S0304-4165(20)30128-8/rf0160)
- <span id="page-10-16"></span>[33] A. Jain, Compensating mass matrix potential for constrained molecular dynamics, J. Comput. Phys. 136 (1997) 289–297, [https://doi.org/10.1006/jcph.1997.5731.](https://doi.org/10.1006/jcph.1997.5731)
- <span id="page-10-17"></span>[34] A. Jain, S. Kandel, J. Wagner, A. Larsen, N. Vaidehi, Fixman compensating potential for general branched molecules, J. Chem. Phys. 139 (2013), [https://doi.org/10.](https://doi.org/10.1063/1.4851315) [1063/1.4851315.](https://doi.org/10.1063/1.4851315)
- <span id="page-10-18"></span>[35] G. Rodriguez, K. Kreutz, A. Jain, A spatial operator algebra for manipulator modeling and control, Proceedings, 1989 Int. Conf. Robot. Autom., IEEE Comput. Soc. Press, 2020, pp. 1374–1379, , [https://doi.org/10.1109/ROBOT.1989.100171.](https://doi.org/10.1109/ROBOT.1989.100171)
- <span id="page-10-19"></span>[36] R.E. Kalman, A new approach to linear filtering and prediction problems, J. Fluids Eng. Trans. ASME. 82 (1960) 35–45, [https://doi.org/10.1115/1.3662552.](https://doi.org/10.1115/1.3662552)
- <span id="page-10-20"></span>[37] [A.E. Bryson, M. Frazier, Smoothing for linear and nonlinear dynamic systems, Proc.](http://refhub.elsevier.com/S0304-4165(20)30128-8/rf0185) [Optim. Syst. Synth. Conf. ASD-TDR-63-1 19, Aeronautical Systems Division, Wright](http://refhub.elsevier.com/S0304-4165(20)30128-8/rf0185) [Patterson AFB, Ohio, DTIC Document, 1963, pp. 353](http://refhub.elsevier.com/S0304-4165(20)30128-8/rf0185)–364.
- <span id="page-10-21"></span>[38] Rodriguez, Kalman filtering, smoothing and recursive robot arm foward and inverse dynamics, IEEE J. Robot. Autom. 3 (1987) 624–639, [https://doi.org/10.1109/JRA.](https://doi.org/10.1109/JRA.1987.1087147) [1987.1087147.](https://doi.org/10.1109/JRA.1987.1087147)
- <span id="page-10-22"></span>[39] G. Rodriguez, A. Jain, K. Kreutz-Delgado, A spatial operator algebra for manipulator Modeling and control, Int. J. Robot. Res. 10 (1991) 371–381, [https://doi.org/10.](https://doi.org/10.1177/027836499101000406) [1177/027836499101000406.](https://doi.org/10.1177/027836499101000406)
- <span id="page-10-23"></span>[40] A. Jain, Unifi[ed Formulation of Dynamics for Serial Rigid Multibody Systems,](http://refhub.elsevier.com/S0304-4165(20)30128-8/rf0200) [\(1991\).](http://refhub.elsevier.com/S0304-4165(20)30128-8/rf0200)
- <span id="page-10-24"></span>[41] M.A. Sherman, A. Seth, S.L. Delp, Simbody: multibody dynamics for biomedical research, Procedia IUTAM. 2 (2011) 241–261, [https://doi.org/10.1016/j.piutam.](https://doi.org/10.1016/j.piutam.2011.04.023)

- [2011.04.023.](https://doi.org/10.1016/j.piutam.2011.04.023)
- <span id="page-10-25"></span>[42] P. Eastman, J. Swails, J.D. Chodera, R.T. McGibbon, Y. Zhao, K.A. Beauchamp, L.- P. Wang, A.C. Simmonett, M.P. Harrigan, C.D. Stern, R.P. Wiewiora, B.R. Brooks, V.S. Pande, OpenMM 7: rapid development of high performance algorithms for molecular dynamics, PLoS Comput. Biol. 13 (2017) e1005659, , [https://doi.org/10.](https://doi.org/10.1371/journal.pcbi.1005659) [1371/journal.pcbi.1005659.](https://doi.org/10.1371/journal.pcbi.1005659)
- <span id="page-10-26"></span>[43] A. Jain, I.H. Park, N. Vaidehi, Equipartition principle for internal coordinate molecular dynamics, J. Chem. Theory Comput. 8 (2012) 2581–2587, [https://doi.org/](https://doi.org/10.1021/ct3002046) [10.1021/ct3002046.](https://doi.org/10.1021/ct3002046)
- <span id="page-10-27"></span>[44] C.S. Jensen, Blocking Gibbs Sampling for Inference in Large and Complex Bayesian Networks with Applications in Genetics, Aalborg University, 1997, [https://doi.org/](https://doi.org/10.1086/302524) [10.1086/302524.](https://doi.org/10.1086/302524)
- <span id="page-10-28"></span>[45] E. Cancès, F. Castella, P. Chartier, E. Faou, C. Le Bris, F. Legoll, G. Turinici, Longtime averaging for integrable Hamiltonian dynamics, Numer. Math. 100 (2005) 211–232, [https://doi.org/10.1007/s00211-005-0599-0.](https://doi.org/10.1007/s00211-005-0599-0)
- <span id="page-10-29"></span>[46] A. Beskos, N. Pillai, G. Roberts, J.M. Sanz-Serna, A. Stuart, Optimal tuning of the hybrid Monte Carlo algorithm, Bernoulli. 19 (2013) 1501–1534, [https://doi.org/](https://doi.org/10.3150/12-BEJ414) [10.3150/12-BEJ414.](https://doi.org/10.3150/12-BEJ414)
- <span id="page-10-30"></span>[47] M. Betancourt, Identifying the Optimal Integration Time in Hamiltonian Monte Carlo, (2016), p. 31 [http://arxiv.org/abs/1601.00225.](http://arxiv.org/abs/1601.00225)
- <span id="page-10-31"></span>[48] M. Betancourt, A Unified Treatment of Predictive Model Comparison, (2015), pp. 1–20 [http://arxiv.org/abs/1506.02273.](http://arxiv.org/abs/1506.02273)
- <span id="page-10-32"></span>[49] P. Echenique, I. Calvo, J.L. Alonso, Quantum mechanical calculation of the effects of stiff and rigid constraints in the conformational equilibrium of the alanine dipeptide, J. Comput. Chem. 27 (2006) 1733–1747, [https://doi.org/10.1002/jcc.](https://doi.org/10.1002/jcc.20467) [20467.](https://doi.org/10.1002/jcc.20467)
- <span id="page-10-33"></span>[50] A. Seth, M. Sherman, P. Eastman, S. Delp, Minimal formulation of joint motion for biomechanisms, Nonlinear Dyn. 62 (2010) 291–303, [https://doi.org/10.1007/](https://doi.org/10.1007/s11071-010-9717-3) [s11071-010-9717-3.](https://doi.org/10.1007/s11071-010-9717-3)
- <span id="page-10-34"></span>[51] J. Prentoe, R. Velázquez-Moctezuma, E.H. Augestad, A. Galli, R. Wang, M. Law, H. Alter, J. Bukh, Hypervariable region 1 and N-linked glycans of hepatitis C regulate virion neutralization by modulating envelope conformations, Proc. Natl. Acad. Sci. U. S. A. 116 (2019) 10039–10047, [https://doi.org/10.1073/pnas.1822002116.](https://doi.org/10.1073/pnas.1822002116)
- <span id="page-10-36"></span>[52] L. Deng, L. Ma, M.L. Virata-Theimer, L. Zhong, H. Yan, Z. Zhao, E. Struble, S. Feinstone, H. Alter, P. Zhang, Discrete conformations of epitope II on the hepatitis C virus E2 protein for antibody-mediated neutralization and nonneutralization, Proc. Natl. Acad. Sci. U. S. A. 111 (2014) 10690–10695, [https://doi.org/10.1073/](https://doi.org/10.1073/pnas.1411317111) [pnas.1411317111.](https://doi.org/10.1073/pnas.1411317111)
- <span id="page-10-35"></span>[53] B. Webb, A. Sali, Comparative Protein Structure Modeling Using MODELLER, Curr. Protoc. Bioinforma, John Wiley & Sons Inc, Hoboken, NJ, USA, 2016, pp. 5.6.1–5.6.37, , [https://doi.org/10.1002/cpbi.3.](https://doi.org/10.1002/cpbi.3)
- <span id="page-10-37"></span>[54] D.A. Case, R.M. Betz, D.S. Cerutti, T.E. Cheatham III, T.A. Darden, R.E. Duke, T.J. Giese, H. Gohlke, A.W. Goetz, N. Homeyer, S. Izadi, P. Janowski, J. Kaus, A. Lovalenko, T. Lee, S. LeGrand, P. Li, C. Lin, T. Luchko, R. Luo, B. Madej, D. Mermelstein, K. Merz, G. Monard, H. Nguyen, H. Nguyen, I. Omelyan, A. Onufriev, D. Roe, A. Roitberg, C. Sagui, C. Simmerling, W. Botello-Smith, J. Swails, R. Walker, J. Wang, R. Wolf, X. Wu, L. Xiao, P. Kollman, AMBER 2016, Univ. Calif, San Fr. CA, USA, 2016, pp. 1–923, [https://doi.org/10.13140/RG.2.2.](https://doi.org/10.13140/RG.2.2.27958.70729) [27958.70729.](https://doi.org/10.13140/RG.2.2.27958.70729)
- <span id="page-10-38"></span>[55] J.A. Maier, C. Martinez, K. Kasavajhala, L. Wickstrom, K.E. Hauser, C. Simmerling, ff14SB: Improving the Accuracy of Protein Side Chain and Backbone Parameters from ff99SB, J. Chem. Theory Comput. 11 (2015) 3696–3713, [https://doi.org/10.](https://doi.org/10.1021/acs.jctc.5b00255) [1021/acs.jctc.5b00255.](https://doi.org/10.1021/acs.jctc.5b00255)
- <span id="page-10-39"></span>[56] K.N. Kirschner, A.B. Yongye, S.M. Tschampel, J. González-Outeiriño, C.R. Daniels, B.L. Foley, R.J. Woods, GLYCAM06: A generalizable biomolecular force field. carbohydrates, J. Comput. Chem 29 (2008) 622–655, [https://doi.org/10.1002/jcc.](https://doi.org/10.1002/jcc.20820) [20820.](https://doi.org/10.1002/jcc.20820)
- <span id="page-10-40"></span>[57] W. Humphrey, A. Dalke, K. Schulten, VMD: Visual molecular dynamics, J. Mol. Graph. 14 (1996) 33–38, [https://doi.org/10.1016/0263-7855\(96\)00018-5.](https://doi.org/10.1016/0263-7855(96)00018-5)
- <span id="page-10-41"></span>[58] S. Van Der Walt, S.C. Colbert, G. Varoquaux, The NumPy array: a structure for efficient numerical computation, Comput. Sci. Eng. 13 (2011) 22–30, [https://doi.](https://doi.org/10.1109/MCSE.2011.37) [org/10.1109/MCSE.2011.37.](https://doi.org/10.1109/MCSE.2011.37)
- <span id="page-10-42"></span>[59] K.J. Millman, M. Aivazis, Python for scientists and engineers, Comput. Sci. Eng. 13 (2011) 9–12, [https://doi.org/10.1109/MCSE.2011.36.](https://doi.org/10.1109/MCSE.2011.36)
- <span id="page-10-43"></span>[60] J.D. Hunter, Matplotlib: a 2D graphics environment, Comput. Sci. Eng. 9 (2007) 99–104, [https://doi.org/10.1109/MCSE.2007.55.](https://doi.org/10.1109/MCSE.2007.55)
- <span id="page-10-44"></span>[61] [Microsoft Corporation, Microsoft Excel, \(2018\).](http://refhub.elsevier.com/S0304-4165(20)30128-8/rf0305)
- <span id="page-10-45"></span>[62] J.H. Ward, Hierarchical grouping to optimize an objective function, J. Am. Stat. Assoc. 58 (1963) 236–244, [https://doi.org/10.1080/01621459.1963.10500845.](https://doi.org/10.1080/01621459.1963.10500845)
- <span id="page-10-46"></span>[63] B. Montgomery Pettitt, M. Karplus, The potential of mean force surface for the alanine dipeptide in aqueous solution: a theoretical approach, Chem. Phys. Lett. 121 (1985) 194–201, [https://doi.org/10.1016/0009-2614\(85\)85509-3.](https://doi.org/10.1016/0009-2614(85)85509-3)
- <span id="page-10-47"></span>[64] A. Varki, Biological roles of glycans, Glycobiology. 27 (2017) 3–49, [https://doi.org/](https://doi.org/10.1093/glycob/cww086) [10.1093/glycob/cww086.](https://doi.org/10.1093/glycob/cww086)
- <span id="page-10-48"></span>[65] E. Marinari, G. Parisi, Simulated Tempering: A New Monte Carlo Scheme, (1992), [https://doi.org/10.1209/0295-5075/19/6/002.](https://doi.org/10.1209/0295-5075/19/6/002)
- <span id="page-10-49"></span>[66] R.H. Swendsen, J.S. Wang, Replica Monte Carlo simulation of spin-glasses, Phys. Rev. Lett. 57 (1986) 2607–2609, [https://doi.org/10.1103/PhysRevLett.57.2607.](https://doi.org/10.1103/PhysRevLett.57.2607)