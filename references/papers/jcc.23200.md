![](_page_0_Picture_2.jpeg)

# Advanced Techniques for Constrained Internal Coordinate Molecular Dynamics

Jeffrey R. Wagner,[a] Gouthaman S. Balaraman,[a] Michiel J. M. Niesen,[a] Adrien B. Larsen,[a] Abhinandan Jain,[b] and Nagarajan Vaidehi\*[a]

Internal coordinate molecular dynamics (ICMD) methods provide a more natural description of a protein by using bond, angle, and torsional coordinates instead of a Cartesian coordinate representation. Freezing high-frequency bonds and angles in the ICMD model gives rise to constrained ICMD (CICMD) models. There are several theoretical aspects that need to be developed to make the CICMD method robust and widely usable. In this article, we have designed a new framework for (1) initializing velocities for nonindependent CICMD coordinates, (2) efficient computation of center of mass velocity during CICMD simulations, (3) using advanced integrators such as Runge–Kutta, Lobatto, and adaptive CVODE for CICMD simulations, and (4) cancelling out the ''flying ice cube effect'' that sometimes arises in Nose–Hoover dynamics. The Generalized Newton–Euler Inverse Mass Operator (GNEIMO) method is an implementation of a CICMD method that we have developed to study protein dynamics. GNEIMO allows for a hierarchy of coarse-grained simulation models based on the ability to rigidly constrain any group of atoms. In this article, we perform tests on the Lobatto and Runge–Kutta integrators to determine optimal simulation parameters. We also implement an adaptive coarse-graining tool using the GNEIMO Python interface. This tool enables the secondary structure-guided ''freezing and thawing'' of degrees of freedom in the molecule on the fly during molecular dynamics simulations and is shown to fold four proteins to their native topologies. With these advancements, we envision the use of the GNEIMO method in protein structure prediction, structure refinement, and in studying domain motion. V*<sup>C</sup>* 2013 Wiley Periodicals, Inc.

DOI: 10.1002/jcc.23200

## Introduction

The desire to use larger integration timesteps has motivated the use of constraints in molecular dynamics (MD) simulations to eliminate high-frequency degrees of freedom.[1–9] The MD techniques that impose such bond length constraints include the SHAKE[10] and RATTLE[11] algorithms. Both of these algorithms solve the all-atom (referred to as Cartesian hereafter) equations of motion, followed by an iterative solution of the constraint equations to enforce the bond length constraints. These algorithms are available in widely used software packages such as CHARMM, AMBER, and NAMD. Internal coordinate molecular dynamics (ICMD) on the other hand use bond/ angle/torsional (BAT) coordinates which more naturally reflect the large- and small-scale motion degrees of freedom within the molecule. The well known ''torsional molecular dynamics'' technique is a special case of constrained ICMD (CICMD) models which freeze bond lengths and angles.[1,3,6,8,4,9,12] In CICMD models, the molecule is modeled as a collection of rigid clusters connected by hinges. Each cluster is a rigid collection of atoms, within which all bond lengths and bond angles are frozen.

An issue that has received considerable attention from researchers has been the increased complexity of the equations of motion for CICMD models. While the number of degrees of freedom is smaller, the mass matrix is dense and configuration dependent. The computational cost for some of the original approaches to solving the constrained equations of motion scale as the cube of the number of degrees of freedom.[1,2,5] We developed a spatial operator algebra- (SOA)[13,14] based Generalized Newton–Euler Inverse Mass Operator (GNEIMO) method for CICMD[3] that solves the same equations of motion exactly with OðN Þ computational cost, where N denotes the number of degrees of freedom. The SOA algorithm and its variants have also been used by other groups for CICMD simulations.[4,8,9,12,15]

In our recent work, we have demonstrated the application of CICMD models and trajectories for the folding of proteins,[16] the refinement of protein structure,[17] and the simulation of protein domain motion.[18] We have observed that CICMD models are able to fold small proteins faster and more reliably than using Cartesian models.[16] Further, we found that CICMD models are able to refine protein homology models to higher accuracy consistently and enrich the population of refined structure.[17] Finally, CICMD simulations of proteins known to undergo large-scale domain motion have managed

![](_page_0_Picture_18.jpeg)

<sup>[</sup>a] J. R. Wagner, G. S. Balaraman, M. J. M. Niesen, A. B. Larsen, N. Vaidehi Division of Immunology, Beckman Research Institute of the City of Hope, Duarte, California 91010 E-mail: nvaidehi@coh.org

<sup>[</sup>b] A. Jain

Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California 91109

Contract/grant sponsor: National Institute of Health; Contract/grant number: RO1GM082896; Contract/grant sponsor: NIH RO1GM082896 is the only grant that funded this project and not the SimBios grant.

V<sup>C</sup> 2013 Wiley Periodicals, Inc.

![](_page_1_Picture_2.jpeg)

to reproduce the expected conformational changes, while Cartesian simulations do not.[18] These applications exploit the qualitatively different nature of CICMD models to obtain performance improvements over traditional Cartesian models beyond just increased integration timestep size. One reason for the different behavior of CICMD models is that the number of degrees of freedom in CICMD models is approximately one order of magnitude smaller than that in traditional Cartesian MD models. We believe that the use of the more natural BAT coordinate models and the retention of the essential degrees of freedom play a significant role in the superior performance of CICMD models. We are continuing to work on extending and applying the ICMD and CICMD methodology to a broader range of MD applications.

The added complexity of the equations of motion resulting from the use of BAT coordinates and holonomic constraints is just one of several aspects, where the CICMD models differ from unconstrained Cartesian MD models. There remain several outstanding technical issues with CICMD that need to be addressed to make the method a robust and stable dynamics simulation tool for use with the wider variety of molecular systems. In this article, we have addressed some of these key issues by developing the theoretical framework and associated computational algorithms for implementation.

These are summarized below.

- The equipartition principle is used as the basis for the initialization of atom velocities in accordance with Boltzmann distributions in Cartesian MD. However, the traditional equipartition principle does not hold for CICMD models. This is apparent by observing that the velocities of atoms within the same cluster are completely correlated. In a companion work, we have derived a rigorous new equipartition principle for CICMD models that makes use of new modal velocity coordinates.[19] In the Methods section, we describe a low-cost computational procedure for carrying out velocity initialization based on the CICMD equipartition principle and the modal coordinates.
- While MD simulations are typically initialized with zero overall system linear and rotational momentum, there is inevitable drift in practice, and periodic nulling of linear and angular motion is needed. We have developed a composite body inertia-based technique for efficiently carrying out such momentum nulling for CICMD models, with the mathematical details included in the Methods section.
- We observed the emergence of a flying ice cube behavior during CICMD simulations that results in a bleeding of thermal energy from the configuration degrees of freedom. In the Methods section, we describe a mathematical explanation for this behavior and techniques for avoiding this phenomenon.
- The simplicity of Cartesian models has facilitated the use of energy conserving integration techniques such as Verlet integrators for Cartesian MD. Conversely, the presence of holonomic constraints in CICMD results in a nonseparable Hamiltonian (i.e., the kinetic energy depends on the configuration coordinates) and as a consequence, the need for alternative integration techniques. We have analyzed the performance of the Runge–Kutta 4 (RK4) and Lobatto integrators for long-term stable CICMD simulations in the Integrator Performance section.

• We have developed the GNEIMO simulation platform in which the user can freeze and thaw any degree(s) of freedom of the protein model to perform all-atom MD, all-torsion MD, or any other form of coarse-grained MD by placing rigid constraints on the appropriate degrees of freedom. Thus, the GNEIMO method provides an inherent coarse-graining feature that is valuable for enriching conformational sampling during MD simulations. We have demonstrated the use of the freeze and thaw coarse-graining method in studying protein dynamics.[18] In our previous work, however, the coarse-graining model was prescribed in the beginning of the simulation and remained fixed. In this article, we have used an adaptive onthe-fly technique that allows changes in the level of coarsegraining during the simulation. In the Methods section, we describe three levels of coarse-graining techniques for CICMD simulations that range from clustering strategies, to run-time manual freeze and thaw techniques, to automated dynamicclustering strategies. In the Results section, we demonstrate the use of the dynamic clustering algorithm for the folding of four proteins starting from their extended structures to their respective native states by treating the helical and b-strand regions of the proteins as clusters.

We have implemented all the techniques described here in the GNEIMO CICMD simulation package. GNEIMO applies userdefined holonomic constraints onto a Cartesian model of the system to produce a set of generalized coordinates for simulation. We include an architectural description of GNEIMO, including the addition of important MD capabilities such as constant temperature, (N,V,T) Nose–Hoover dynamics, solvent models, and replica exchange MD (REXMD) strategies. The supporting information includes an overview of the GNEIMO software design, as well as summary mathematical derivations of some of the algorithmic techniques used in GNEIMO.

## Method Development

The GNEIMO method for performing CICMD simulations includes several standard capabilities required for MD simulations such as:

- The constant energy (N,V,E) and Nose–Hoover (N,V,T) ensembles in the GNEIMO framework are described in a previous work.[6] For convenience, a summary of this method is included in the supporting information.
- An implementation of the generalized Born solvation method for implicit solvation.[20]
- Support for multiple molecules of any type, including explicit solvent.
  - Implementation of the Fixman correction potential.[21]
- A temperature-based REXMD method,[22] in which temperatures may be switched randomly or probabilistically using the Metropolis algorithm.[23]
  - Support for periodic boundary conditions.
- Architectural integration with the LAMMPS[24] tool. These include the generation of ICMD coordinate files to complement LAMMPS checkpoint files and allow for the restarting of simulations, utilization of Message Passing Interface (MPI) for parallelized force field calculations, and support for standard Cartesian simulation.

![](_page_2_Picture_1.jpeg)

- A GPU-accelerated OpenMM force field implementation, [25] as well as a general set of expected interface functions that could allow GNEIMO to interact with any atomistic force calculator.
- Soft constraints between atom pairs for applying user-defined restraints during simulations.

In the following sections, we describe new techniques and CICMD extensions that have been implemented within the GNEIMO method.

#### Initialization of velocities in modal coordinates

In this article, we focus on torsional MD for systems with tree topology and internal hinges with one degree of freedom. At the start of a simulation, initial velocities need to be assigned to the generalized velocity coordinates. Unlike the Cartesian dynamics case, the degrees of freedom in the constrained dynamics models are coupled. Hence, the Boltzmann distribution cannot be used to directly assign velocities for the constrained dynamics model.

One option is to ignore the cluster model, and to assign the atom velocities using the Cartesian approach and use an ad hoc approximation/projection step to obtain internal coordinate velocities that best approximate the assigned atom velocities. We describe below an alternate method that avoids such ad hoc steps, and instead identifies independent modal velocity degrees of freedom for the constrained dynamics model that can be used to assign velocities according to the Boltzmann distribution. [19]

For a desired temperature, T, the overall thermal energy in the system with  $\mathcal N$  degrees of freedom is defined as

$$\Re_e = \frac{1}{2}(\mathcal{N} - 6)\mathbf{k}\mathbf{T} \tag{1}$$

Here,  ${\bf k}$  is the Boltzmann constant and  $\Re_e$  is the kinetic energy. Using the expression for the mass matrix  ${\cal M}$  derived in the supporting information, the kinetic energy  $\Re_e$  in the system can be expressed as

$$\mathfrak{R}_{e} = \frac{1}{2}\dot{\theta}^{*}\mathcal{M}\dot{\theta} = \frac{1}{2}\dot{\theta}^{*}[I + H\phi\mathcal{K}]\mathcal{D}[I + H\phi\mathcal{K}]^{*}\dot{\theta}$$
(2)

Now define a new set of velocity coordinates,  $v \in \mathcal{R}^{\mathcal{N}}$  as

$$v \triangleq \mathcal{D}^{\frac{1}{2}}[I + H\phi\mathcal{K}]^*\dot{\theta} \tag{3}$$

The above equation represents a velocity coordinate transformation. This transformation is reversible and we can recover  $\dot{\theta}$  from v using the expression

$$\dot{\theta} = [I - H\psi \mathcal{K}]^* \mathcal{D}^{-\frac{1}{2}} v \tag{4}$$

Using this in eq. (2) results in the new kinetic energy expression:

$$\Re_e = \frac{1}{2} v^* v = \frac{1}{2} \sum_{k=1}^{\mathcal{N}} v^*(k) v(k)$$
 (5)

This kinetic energy expression is simply the sum of the kinetic energy contributions of the new v(k) velocity coordinates. We

refer to these independent velocity coordinates as the *modal* velocity coordinates for the constrained dynamics model. It has been shown that the equipartition principle for CICMD models holds for these modal velocity coordinates. <sup>[19]</sup> We assign random velocities according to the Boltzmann distribution in this canonical velocity coordinates, and then recover the  $\dot{\theta}$  value using eq. (4). The algorithm proceeds as follows:

- 1. For an initial desired temperature T, use eq. (1) to obtain a desired  $\Re_e$  target value for the system.
- **2.** Use a zero-mean, unit variance, normal distribution to randomly assign initial values to the  $\nu$  velocity coordinates (except for the six degrees of freedom for the base clusters).
- **3.** Compute  $\dot{\theta}$  from the v values using eq. (4). This can be done via a  $O(\mathcal{N})$  recursive base-to-tips scatter computational algorithm.
- **4.** Reset any nonzero center of mass (CM) velocity resulting from this velocity assignment using the method described below. This step will initialize the velocity of the base cluster's six degree of freedom hinge.
- **5.** Compute the overall kinetic energy in the system, and scale all the velocities so that kinetic energy in the system matches the desired temperature.

#### Resetting the CM velocity in spatial coordinates

In this section, we derive the expressions and algorithms for resetting the CM velocity. The notation here borrows heavily from that in Jain.<sup>[27]</sup>

Simulations are initialized so that the CM of the system has zero translational and rotational velocity. However, due to numerical errors, the system trajectory can accumulate nonzero linear and angular CM velocities. We describe here an efficient procedure developed to reset such nonzero CM velocities (both linear and angular). This method is used during initialization at the start of the simulation, as well as periodically during long simulation runs to remove any accumulated CM energy.

Computing the CM spatial velocity,  $\mathcal{V}_{CM}$  requires computing the overall spatial inertia of the system,  $M_S$  and the overall spatial momentum,  $\mathfrak{h}_S$ , of the system. Both of these quantities depend upon the atomic position coordinates, whereas the momentum also depends on the velocity coordinates. Once these quantities are available,  $\mathcal{V}_{CM}$  can be computed by solving the following linear matrix equation:

$$M_{\rm S} \mathcal{V}_{\rm CM} = \mathfrak{h}_{\rm S}$$
 (6)

Multiple chains are treated as a single system when computing the overall system spatial inertia and spatial momentum. Equation (6) allows us to compute the linear and angular velocity of a frame at the CM. The CM velocity can be nulled out by removing its contribution to the spatial velocity of the base clusters for each of the chains. We have tested this method to null the CM velocity for a single chain system and will be performing the tests for multiple chain systems in the future.

![](_page_2_Picture_31.jpeg)

![](_page_3_Picture_2.jpeg)

Computing system spatial inertia,  $M_s$ . The overall  $6\times 6$  spatial inertia matrix for the system, referenced to the base-cluster, is given by

$$\begin{aligned} M_{S} &= \sum_{k=1}^{n} \phi(n,k) M(k) \phi^{*}(n,k) = [\phi(n,1), \cdots \phi(n,n)] M \begin{bmatrix} \phi^{*}(n,1) \\ \vdots \\ \phi^{*}(n,n) \end{bmatrix} \\ &= E \phi M \phi^{*} E^{*} = E \Big[ \mathcal{R} + \tilde{\phi} \mathcal{R} + \mathcal{R} \tilde{\phi}^{*} \Big] E^{*} = E \mathcal{R} E^{*} = \mathcal{R}(n) \end{aligned}$$

We have used the following facts in the above derivation:

$$E \triangleq [0_6, \dots 0_6, l_6] \in \mathcal{R}^{6 \times 6n},$$

$$E \phi = [\phi(n, 1), \dots \phi(n, n)] \quad \text{and} \quad \tilde{\phi} \mathcal{R} E^* = 0$$
(7)

Here,  $\phi(j,k)$  denotes the 6  $\times$  6 rigid-body transformation matrix between the jth and kth cluster, M(k) is the spatial inertia of the kth cluster,  $\mathcal{R}(k)$  is the 6  $\times$  6 composite rigid body spatial inertia of the kth cluster and all its children, and  $\tilde{\phi}=\phi-I$ , where I denotes the identity matrix. E is also referred to as the base pick-off operator. The first moment vector of the  $\mathcal{R}(n)$  spatial inertia specifies the instantaneous location of the system's CM with respect to the base-cluster's reference frame.

Computing the system spatial momentum,  $\mathfrak{h}_{S}$ . The base-cluster frame referenced spatial momentum of the system is given by

$$\mathfrak{f}_{S} = \sum_{k=1}^{n} \phi(n, k) M(k) \mathcal{V}(k) = E \phi M \mathcal{V} = E \phi M \phi^{*} H^{*} \dot{\theta} 
= E \Big[ \mathcal{R} + \tilde{\phi} \mathcal{R} + \mathcal{R} \tilde{\phi}^{*} \Big] H^{*} \dot{\theta} = E \phi \mathcal{R} H^{*} \dot{\theta}$$

For the base cluster, that is, one whose nth hinge is a full six degree of freedom hinge,  $H^*(n) = I$  and  $\dot{\theta}(n) = \mathcal{V}(n)$ , and hence  $\mathfrak{h}_S$  can be rewritten as

$$\mathfrak{h}_{S} = \mathcal{R}(n)\mathcal{V}(n) + \sum_{k=1}^{n-1} \phi(n,k) \ \mathcal{R}(k)H^{*}(k)\dot{\theta}(k)$$

Computing  $\mathcal{V}_{CM}$  CM spatial velocity. The spatial momentum  $\mathfrak{h}_S$ , and the system level spatial inertia  $M_S$  and the CM velocity  $\mathcal{V}_{CM}$  are (referenced about the base cluster frame) related together by

$$\mathfrak{h}_{s} = M_{s} \mathcal{V}_{CM} = \mathcal{R}(n) \mathcal{V}_{CM} \tag{8}$$

Hence, for an isolated molecular system we have

$$\mathcal{V}_{CM} = \mathcal{R}^{-1}(n) \left[ \mathcal{R}(n) \mathcal{V}(n) + \sum_{k=1}^{n-1} \phi(n,k) \mathcal{R}(k) H^*(k) \dot{\theta}(k) \right]$$
$$= \mathcal{V}(n) + \mathcal{R}^{-1}(n) \sum_{k=1}^{n-1} \phi(n,k) \mathcal{R}(k) H^*(k) \dot{\theta}(k) \tag{9}$$

From eq. 9, it follows that adding  $\delta_V$  to the  $\mathcal{V}(n)$  base-cluster spatial velocity will result in an additional  $\mathcal{R}(n)\delta_V$  of spatial angular momentum to the system. For the spatial momentum to be zero, we must have

$$0 = \mathfrak{h}_{S} + \mathcal{R}(n)\delta_{V} \Rightarrow \delta_{V} = -\mathcal{R}^{-1}(n)\mathfrak{h}_{S} = -\mathcal{V}_{CM}$$

Thus, an additional spatial velocity of  $\delta_V$  for the base-cluster adds  $\mathcal{R}(n)\delta_V$  of base-cluster frame referenced spatial momentum for an isolated molecule. Hence, an additional spatial velocity of  $-\mathcal{V}_{\mathit{CM}}$  at the base-cluster is applied to nullify and reset the system's spatial momentum.

### Flying ice cube effect in Nosé-Hoover dynamics

In MD simulations, an artifact of numerical integration known as the *flying ice cube* effect is a phenomenon, where the energy from high-frequency modes is drained into zero- (translational and rotational) and low-frequency modes.<sup>[28]</sup> As a result of this, the molecule gradually loses internal kinetic energy until the point at which the molecule drifts in space as a rigid body in a frozen conformation. The flying ice cube effect is known to arise in velocity rescaling thermostats such as the Berendsen thermostat<sup>[28]</sup> due to repeated scaling of the kinetic energy of the molecule to maintain the temperature of the simulation.

Here, we report a similar effect that arises in Nosé–Hoover thermostats. For simplicity, we present the derivations in Cartesian coordinates. The velocity equation for the Nosé–Hoover thermostat multiplied by the mass of the *i*th particle is given as

$$\frac{dm_i\mathbf{v}_i}{dt} = \mathbf{F}_i - \zeta m_i\mathbf{v}_i \tag{10}$$

where  $\zeta$  is friction coefficient associated with the bath,  $F_i$ ,  $m_i$  and  $v_i$  are the internal force, mass, and velocities associated with the *i*th particle, respectively. Summing eq. 10 over *i* yields

$$\frac{d\mathcal{M}\mathcal{V}_{CM}}{dt} = -\zeta \mathcal{M}\mathcal{V}_{CM} \tag{11}$$

where  $\mathcal{M}=\sum_i m_i$  is the total mass of the system and  $\mathcal{V}_{CM}$  is the CM velocity of the system. Note that we make use of the fact that the sum of internal force acting on a system is zero, that is,  $\sum_i \mathbf{F}_i = 0$ . Solving eq. (11) yields

$$\mathcal{V}_{CM}(t) = \mathcal{V}_{CM}(t=0) \exp\left[-\int_{0}^{t} \zeta(t')dt'\right]$$

$$= \mathcal{V}_{CM}(t=0) \exp\left[-\ln(s_{t})\right], \tag{12}$$

where  $ln\ (s_t)$  is the bath variable in Nosé formulation that is proportional to the bath potential energy. Ideally, if  $ln\ (s_t)$  the CM velocity at the start of simulation is zero and if remains positive, the CM velocity of the system will remain damped at all times. But the expression for the CM velocity in eq. 12 points to the fact that negative  $ln\ (s_t)$  values can lead to growth in CM velocity leading to a flying ice cube effect in Nosé–Hoover thermostats.

The CM kinetic energy can be written using eq. (12) as

$$KE_{CM} = c \exp(-2 \ln s_t), \tag{13}$$

where *c* is a constant proportional to initial CM kinetic energy.

The CM kinetic energy from eq. (13) is compared for cluster model and all-atom model with and without GB/SA solvation

![](_page_4_Picture_1.jpeg)

![](_page_4_Figure_2.jpeg)

Figure 1. Plot of center of mass (CM) kinetic energy as a function of time for a) cluster model (Cl) using integration timestep 20 fs and b) all atom (AA) Cartesian simulations using integration timestep 1 fs. The kinetic energy of CM from simulations with GB/SA solvation (red line) and vacuum simulation (blue-dashed line) for both cluster and all-atom models have been shown here. The formula c exp(2 In st) for a fitted value of c is shown as squares and triangles for cases with and without GB/SA solvation, respectively. All runs shown here use RK4 integrator. [Color figure can be viewed in the online issue, which is available at wileyonlinelibrary.com.]

in Figures 1a and 1b. The constant c in eq. 13 is determined as a fitting parameter from the simulations. From Figure 1, we observe that the rate of growth of CM kinetic energy is lower in the presence of GB/SA solvation compared with vacuum simulations. We also observe that the rate of growth of CM kinetic energy is far greater in the cluster model than in allatom Cartesian simulations. A potential cause for the larger error could be the larger integration step size (dt ¼ 20 fs) used with cluster model simulations compared with the 1 or 2 fs step size all-atom Cartesian runs. In Figure 2, the value of KECM at the end of 90 ps is shown for different integration step sizes of the cluster model with and without GB/SA solvation. We see that value of KECM for 20 fs timestep has grown by a factor of 10 compared with 1 fs timestep. This shows that the flying ice cube effect is a problem that can occur in constrained dynamical simulations with higher integration timesteps. To overcome this issue, we reset the CM velocity at periodic intervals using the method discussed above.

#### Coarse-graining methods

GNEIMO allows for coarse graining of the dynamic model by freezing any desired degree(s) of freedom. The choice of degrees of freedom to be frozen depends on the nature of the simulation being performed. We have implemented the following three methods to define various scenarios of freezing or thawing chosen degrees of freedom.

Automated clustering model. GNEIMO requires a definition of the degrees of freedom to be fixed in the form of a basic cluster model file. This file is used to construct the protein model and define the generalized coordinates to be used in the simulation. It may contain a completely unconstrained system for all-atom dynamics, frozen bonds and angles for torsional dynamics, or any other level of constraints for specific applications. As a default option, we have developed an automated clustering scheme in which all terminal atoms are added into a rigid cluster with their nonterminal neighbor atom, leading to dynamics using all the torsional angles in the protein. Figure 3 shows the automated default clustering scheme in GNEIMO. The peptide bond dihedral angles are not rigidly constrained in this scheme, and side chains have all torsions free (with the exception of closed rings).

Manual freeze and thaw method. The ''freeze and thaw'' method allows the user to freeze or thaw degrees of freedom during the run beyond those defined in the basic cluster model file. The user may make these constraints take effect at any point during the simulation using the Python interface for GNEIMO. For example, this model may be used to treat a helical regions as rigid clusters connected to the loops with flexible torsions.

Dynamic coarse graining. In addition to the manual ''freeze and thaw'' user interface, we have developed an automated toolkit referred to as ''Dynamic Clustering'' within the GNEIMO framework. It can be used to automatically apply hierarchical clustering schemes on proteins during the MD simulations. The generic all-torsion cluster model of a protein is not always adequate for studying large-scale conformational changes in proteins, and manual control is often impractical. In simulations of the long-timescale processes which are governed by low-frequency modes of a protein, one needs to allow for changes in the clustering model to observe meaningful motion. There are several strategies that can be used when applying constraints during the course of the simulations. Poursina et al. have discussed an adaptive clustering strategy for RNA simulations.[29] Here, we have developed a dynamic

![](_page_4_Figure_10.jpeg)

Figure 2. Plot of KECM at the end of 90 ps for different integration timesteps with GB/SA solvation (red full lines) and in vacuum (dashed blue curve). Nose–Hoover thermostat at a bath temperature of 300 K and bath relaxation constant 500 fs was used with RK4 integrator. [Color figure can be viewed in the online issue, which is available at wileyonlinelibrary.com.]

![](_page_4_Picture_12.jpeg)

![](_page_5_Picture_2.jpeg)

![](_page_5_Figure_3.jpeg)

Figure 3. Colored representation of the automated GNEIMO clustering scheme. Each group of same-colored atoms represents a rigid ''cluster'' which is connected to its neighbors via torsional hinges.

clustering scheme wherein the secondary structure elements of the protein (as detected by STRIDE[30]) may be clustered into larger rigid bodies as they are formed, or released to a finer level of dynamical detail during the simulation depending on a number of criteria. This ''freeze and thaw'' clustering scheme is applied on the fly during the dynamics simulations.

As the default clustering model in GNEIMO leaves all nonring torsions free, we are able to lock the torsional hinges between these clusters to rigidify desired sections of a molecule. Our locking approach is taken from methods described in a previous work.[31] At a user-defined frequency, a STRIDE[30] scan is run on the current state of the protein simulation. This scan generates a list of residues involved in secondary structure. The Dynamic Clustering software then maps the residue numbers of these detected motifs to the specific clusters that define their backbone torsion angles and locks them in the simulation. This software is tied to the GNEIMO–REXMD implementation and the REXMD temperature may be used as a threshold, so that secondary structure will only be locked at high temperature. The user may choose to impose an upper threshold for angular velocity during locking as well, to ensure that no torsion which is undergoing significant motion is locked. Further, the build up of stress forces at frozen hinges can be monitored and used to unlock them.

## Results

Our goal is to develop constrained dynamics techniques in GNEIMO that would allow simulations ranging from constraint-free all-atom to large numbers of constraints that give rise to rigid body clusters for the study of domain motion in proteins. We envision using these tools for long-timescale MD simulations to study protein dynamics,[18] protein structure prediction and refinement,[17] folding of simple proteins,[16] and protein–protein complexes. Towards this effort, we have validated the GNEIMO method with dynamic clustering for ab initio prediction of protein structures. We have also performed long-time dynamics simulations starting from crystal structures for different choices of integrators and integration timesteps.

## Protein structure prediction using the dynamic clustering algorithm

In this section, we discuss the utility of the dynamic clustering scheme for the folding of four proteins starting from an extended structure. We first predict the secondary structure region using the ''PSIPRED'' tool[32] and build the predicted helical regions of the sequence. We then start the GNEIMO–REXMD simulations using this sequence in a linear conformation (e.g., the extended structure shown in the top-right corner of Fig. 4(a)) with the helical regions treated as clusters. We used 12 replicas with temperatures ranging from 300 to 1050 K. Temperatures switch randomly within this range at 7.5 ps intervals. Helices and b-strands that are detected by STRIDE[30] during temperature exchange events have their backbone torsions frozen if the REXMD temperature chosen is greater than 400 K. The temperature range was chosen to cover high temperature for improving topological sampling and low temperature to encourage the growth of low energy native like structures. The protein folds to the molten globule state in this process.

The adaptive timestep CVODE[33] integrator using the Adams–Moulton method was used in simulations for dynamic clustering. The adaptive timestep allows for rapid, stable sampling of topological space by large-clustered bodies, while taking smaller steps to allow for more precise dynamics as the protein collapses toward a possible fold. To produce meaningful and comparable results, we use a test set of proteins similar to a previous study on structure prediction.[34] We tested the dynamic clustering scheme on the experimentally resolved portions of four proteins with PDB IDs: 1BDD (res 11-56), 1EON (res 7-31), 1PRB (res 11-53), and 1UBQ (res 1-35). The experimental crystal structures of these proteins are shown in Figures 5d–5g. The starting structure for each dynamic clustering simulation contained only predicted secondary structure[32] and was in the extended conformation otherwise.

Figure 4 shows the population density histogram for the four proteins, taken from the GNEIMO–REXMD trajectory that led to the best folded structure. The quantity on the x-axis is the root mean square deviation in Cartesian coordinates of the

![](_page_6_Picture_1.jpeg)

![](_page_6_Figure_2.jpeg)

Figure 4. a) Backbone CRMSD histogram of a ''Dynamic Clustering'' replica-exchange simulation of 1BDD, beginning from an extended conformation containing only predicted secondary structure elements. Helices that are treated as rigid bodies are shown as broad ribbons. b–d) Backbone CRMSD histograms of dynamic clustering simulations of proteins from predicted helical structure. [Color figure can be viewed in the online issue, which is available at wileyonlinelibrary.com.]

backbone atoms to the corresponding crystal structures (henceforth denoted as CRMSD). It is seen that maximum population of the conformations falls between 5 and 7 A˚ for 1BDD, between 6 and 10 A˚ for 1UBQ, between 7 and 8 A˚ for 1EON, and between 8 and 10 A˚ for 1PRB. The closest structure to the crystal structure in each case has a backbone CRMSD of 4.007 A˚ (1BDD), 4.198 A˚ (1EON), 3.726 A˚ (1PRB), and 4.325 A˚ (1UBQ).

Twelve randomly-switched replicas for each protein were simulated. Each replica of 1PRB and 1UBQ was run for 3 ns (total sim time ¼ 36 ns). Since velocity reinitialization is required at each replica exchange when the freezing or thawing of clusters occurs, the protein has a chance of losing some secondary structure before equilibration is reached. It was found that most of the starting secondary structure was lost early in these simulations, leading to poor refinement after 1 ns. For this reason, each replica of 1BDD and 1E0N was run for just 0.3ns (total sim time ¼ 3.6 ns). Similar folding results were observed for these proteins. Since less overall processor time is required for the shorter simulations, more replicas can be run using the same resources.

Figure 4(a) shows representative structures of one of the 1BDD replicas that folded. The protein begins at an extended structure containing only predicted helices, then begins to sample small numbers of interhelical contacts at backbone CRMSDs in a broad range from 12 to 16 A˚. At 8–11 A˚, the protein samples incorrect packings of all three helices. Below 7 A˚, the helices are sampling the proper native topology. This demonstrates that the dynamic clustering approach samples the correct three dimensional topology within 40 ns of total simulation time. It should

![](_page_7_Picture_2.jpeg)

![](_page_7_Picture_3.jpeg)

Figure 5. a–c) The crystal structures of the three proteins of various sizes used for equilibrium dynamics simulations and d–g) the crystal structures of the four proteins (captioned with the experimentally-resolved residue subrange used in simulations) used for ab initio structure prediction. [Color figure can be viewed in the online issue, which is available at wileyonlinelibrary.com.]

be noted that while GNEIMO–REXMD can be used for structure prediction, the dynamics are not reflective of the pathway of folding since REXMD enhanced sampling techniques have been used and the effect of the dynamic clustering method on the folding pathway has yet to be studied. However, dynamic clustering can be used in other applications such as studying conformational transitions in proteins, a task for which the standard GNEIMO clustering scheme has been shown to outperform unconstrained Cartesian simulations.[18] While using dynamic clustering or treating full helices as rigid bodies, we suggest the use of the CVODE[33] integrator for stable dynamics.

### MD simulations of crystal structures of proteins

Here, we discuss the results of GNEIMO torsional dynamics simulations of three proteins of various sizes with high-resolution crystal structures. The three proteins are Crambin (PDB ID: 1CRN and resolution 1.50 A˚), Defensin (PDB ID: 1DFN and resolution 1.9 A˚), and Bovine Pancreatic Trypsin Inhibitor (BPTI, PDB ID: 4PTI and resolution 1.50 A˚), shown in Figures 5a–5c. GNEIMO (N,V,T) all torsion dynamics simulation for 5 ns at room temperature (310 K) was performed for all the three proteins starting from their respective crystal structures after 500 ps of simulated equilibration. A Hoover thermostat bath relaxation constant of 250 fs was used. We tested various timestep sizes (1–30 fs) and two integrators (RK4 and Lobatto) to evaluate

the integrator performance for room temperature all torsion dynamics. These simulations were performed to analyze the performance of GNEIMO in temperature equilibration and maintaining structural stability.

Integrator performance. Here, we discuss the accuracy of the fixed-step Lobatto and RK4 integrators for the simulations. The Lobatto integrator is the implicit Lobatto IIIa-b partitioned RK method which is an adaptation of the explicit Stormer–Verlet symplectic method for Cartesian MD.[15] The fixed step-size RK4 method is an instance of a standard high-order explicit method. Their accuracies have been assessed by monitoring the standard deviation in simulation temperature for the (N,V,T) canonical ensemble and deviations of the structure from the starting crystal structures as detailed below. Figure 6 shows the standard deviation in temperature for 5 ns simulations using various timestep sizes for the Lobatto and RK4 integrators. We have also shown the results for all-atom Cartesian simulations in the same figure. While multiple-chain thermostats are commonly used for improving the ergodicity of simulations,[35] all simulations in this article used the single-chain Nose–Hoover thermostat with a relaxation time of 250 fs.

We must keep in mind that the Lobatto integrator is second order, whereas the RK4 is fourth order. Though the Lobatto integrator must compute thermostat and Coriolis forces twice per timestep, it only requires one expensive position-dependent

![](_page_8_Picture_1.jpeg)

![](_page_8_Figure_2.jpeg)

Figure 6. Standard deviation of temperature vs. timestep size over a 5 ns simulation for various molecule-integrator combinations. [Color figure can be viewed in the online issue, which is available at wileyonlinelibrary.com.]

molecular force field calculation per timestep. Thus, the RK4 integrator makes four atomistic force field computations each timestep compared to the Lobatto integrator's one. To provide a meaningful comparison, we plot simulation metrics against ''normalized'' timesteps (fs per force computation), as this metric corresponds more closely to wall time. In the tests, the RK4 integrator was stable using simulation timesteps of up to 16 fs, compared to the Lobatto's 10 fs. When converted to normalized timesteps, however, this RK4 simulation takes a normalized timestep of 4 fs compared to the Lobatto's 10 fs.

In Figure 6, we see that integration error (as measured by variation in temperature) is higher for the GNEIMO torsional dynamics model than for the Cartesian model, and that the Cartesian model is not generally stable for timesteps larger than 2 fs. It is seen that the average standard deviation in temperature for the RK4 integrator is relatively flat for normalized timesteps under 4 fs. The Lobatto integrator error rises directly with timestep size, and the simulations begin to fail at normalized timesteps larger than 9 or 10 fs.

When using the Nose–Hoover thermostat, we observed a slow, long-term increase in the magnitude of the bath potential energy. We believe that this is due to the dissipative nature of the nonconserving integrators was used for the CICMD simulations. Future work will investigate the use of multiple Nose– Hoover chains,[35] energy conserving integrators for CICMD simulations, and the use of Nose–Poincare[36] techniques for constant temperature (N,V,T) simulations.

Structural properties. The structural properties of the proteins from the simulation trajectories were analyzed relative to their respective crystal structures. We expect that stable dynamics of these proteins should lead to sampling of the energetically allowed regions of torsional conformation space. We have calculated the average CRMSD of the backbone atoms of the protein from the crystal structures for each of the trajectories. Figure 7 shows the average CRMSD from simulations run with various timestep for the three proteins: crambin, BPTI, and defensin. Simulations which encountered a crash before the full 5 ns were reached have been excluded from the plot. It is observed that the mean CRMSD is less than 2.5 A˚ for most torsional dynamics simulations and that, for the Lobatto and RK4 integrators, the

average simulation CRMSD is not correlated with timestep size. The Cartesian simulations of crambin and BPTI moved away from the crystal structure, whereas the Cartesian simulation of 1DFN remained close to the folded state.

RMSF. Another test of simulation quality is to compare the fluctuation in various regions of the structure to that found experimentally in the crystal structures of the proteins. We have calculated a quantity called root mean square fluctuation (RMSF) for each residue in the simulation. The RMSF is calculated by taking the deviation of each residue's position from its average position during the simulation. We compare this to the RMSF that would be expected from the crystal structure using the relation B ¼ <sup>8</sup>p<sup>2</sup> <sup>3</sup> RMSF2, where B is the crystallographic B-factor for each residue.

Figure 8 shows the difference between the simulation RMSF and the RMSF derived from the crystal structure's B-factor for each residue, timestep size, and for both the RK4 and Lobatto integrators. Each horizontal slice of these plots represents a single simulation, showing the protein sequence from N- to C-terminus with color indicating the difference between the simulated and experimentally-derived RMSF. Horizontal slices that are consistently close to zero represent integrator-timestep combinations that replicate the experimentally expected dynamics. We note again that there are few trends that correlate consistently with timestep size. This result implies that the choice of simulation timestep and integrator have little bearing on the results of GNEIMO simulations, supplied that the timestep does not lead to a crash. Since all simulations appear to have similar outcome, we see that Lobatto is the most efficient choice for performing accurate simulations.

## Conclusions

There are several theoretical and computational issues unique to CICMD models that need to be addressed for their proper use as an MD technique. Some of these issues are addressed in this article, including a rigorous method to initialize velocities in generalized coordinates, an efficient algorithm for nulling the CM velocity, and a proposed explanation and solution

![](_page_8_Figure_14.jpeg)

Figure 7. Average CRMSD in coordinates vs. timestep size over a 5 ns simulation for various molecule-integrator combinations. [Color figure can be viewed in the online issue, which is available at wileyonlinelibrary.com.]

![](_page_8_Picture_16.jpeg)

![](_page_9_Picture_2.jpeg)

![](_page_9_Figure_3.jpeg)

Figure 8. Simulation RMSF—B-factor derived RMSF for various proteins, integrators, and timesteps. RMSFs are calculated from the entire trajectory of each 5 ns simulation. A value of zero indicates that the RMSF of the residue on the x-axis observed in a simulation using the timestep on the y-axis matches the RMSF derived from the crystallographic B-factor. All distances in A˚.

for the flying ice cube effect in CICMD methods. We have developed a rigorous theoretical and computational framework and implemented these features in the GNEIMO CICMD method. The GNEIMO method supports various force fields (CHARMM and AMBER), thermostats (Nose–Hoover, Berendsen, and temperature rescaling thermostats), and solvation (GB/SA and distance dependent dielectric). Long-timescale torsional dynamics simulations of proteins using large timesteps lead to the flying ice cube effect with the Nose–Hoover thermostat unless CM kinetic energies are periodically nulled.

We have implemented and tested advanced integrators such as RK4, Lobatto, and adaptive timestep CVODE.[33] Simulations of three crystal structures were performed with various timestep sizes for 5 ns using the Lobatto and RK4 integrators. Both integrators were found to be stable and replicated crystal flexibility for normalized (fs per force field computation) timesteps as large as 10 fs for Lobatto and 4 fs for RK4. Within this range, the measured stability metrics from both types of integrator appeared to be largely independent of timestep size used. We conclude that the first-order Lobatto integrator allows for the largest normalized timestep and is thus the more efficient integration method for GNEIMO simulations.

We have also developed a toolkit to apply various clustering strategies which provides the user with several levels of coarse

![](_page_10_Picture_1.jpeg)

graining for the dynamic model of the protein. We have implemented the automated clustering of secondary structure elements as they are formed during the folding of proteins. This dynamic clustering tool is an adaptive scheme for coarse graining. Using the dynamic clustering tool, we have folded four different proteins starting from their extended structure to molten globule-like native structures within 4–5 A˚ of the crystal. In the GNEIMO software, we have created a framework for CICMD simulations that can be used in applications such as protein structure refinement, the study of domain motions in proteins, and ab initio structure prediction.

Some of the advancements required for robust CICMD simulations have been developed and implemented in the GNEIMO code but have yet to be tested numerically. These are: (1) dynamics of multiple chains including explicit solvent molecules and (2) use of the Fixman correction potential[21,37] to eliminate the systematic bias that stems from the holonomic constraints placed on the models in the calculation of thermodynamic properties from CICMD simulations.[5,38] We have derived a computational framework for calculating the Fixman compensating potential using the SOA techniques.[39] The numerical testing of the effect of this correction potential and comparison of its performance to force fields that have been tailored for CICMD simulations[40,4] are planned for future work.

## Acknowledgments

The authors thank Dr. Amy Swain (NIGMS) and Dr. Paul Brahznik (NIGMS) for their support and encouragement. Part of the research described in this article was performed at the Jet Propulsion Laboratory (JPL), California Institute of Technology, under contract with the National Aeronautics and Space Administration. The authors thank Simbios for providing us with the GB/SA solvation module, Mark Friedrichs for his help with validating our GB/SA module, and Michael Sherman and Christopher Bruns for helping us with running Simbody. The Simbios software was made freely available on http://simtk.org/home/openmm by the Simbios NIH National Center for Biomedical Computing.

Keywords: constrained dynamics torsional dynamics ab-initio folding molecular dynamics integration techniques GNEI-MO

How to cite this article: J. R. Wagner, G. S. Balaraman, M. J. M. Niesen, A. B. Larsen, A. Jain, N. Vaidehi, J. Comput. Chem. 2013, 34, 904–914. DOI: 10.1002/jcc.23200

Additional Supporting Information may be found in the online version of this article.

- [1] A. Mazur, A. R. Abagyan, J. Biomol. Struct. Dyn. 1989, 6, 815.
- [2] A. K. Mazur, V. E. Dorofeev, R. A. Abagyan, J. Comput. Phys. 1991, 92, 261.
- [3] A. Jain, N. Vaidehi, G. Rodriguez, J. Comput. Phys. 1993, 106, 258.
- [4] J. Chen, W. Im, C. L. Brooks, III, J. Comput. Chem. 2005, 26, 1565.
- [5] K. D. Gibson, H. A. Scheraga, J. Comput. Chem. 1990, 11, 468.
- [6] N. Vaidehi, A. Jain, W. Goddard, J. Phys. Chem 1996, 100, 10508.
- [7] N. Vaidehi, W. A. Goddard, J. Phys. Chem. A 2000, 104, 2375.
- [8] C. D. Schwieters, G. M. Clore, J. Magn. Reson. 2001, 152, 288. [9] S. C. Flores, M. A. Sherman, C. M. Bruns, P. Eastman, R. B. Altman, IEEE/ ACM Trans. Comput. Biol. Bioinform. 2011, 8, 1247.
- [10] J. Ryckaert, G. Ciccotti, H. J. C. Berendsen, J. Comput. Phys. 1977, 23, 327.
- [11] H. C. Andersen, J. Comput. Phys. 1983, 52, 24.
- [12] D. J. Tobias, C. L. Brooks, III, J. Comput. Phys. 1988, 89, 5115.
- [13] A. Jain, J. Guid. Control Dynam. 1991, 14, 531.
- [14] G. Rodriguez, K. Kreutz-Delgado, A. Jain, Int. J. Robot. Res. 1991, 10, 371.
- [15] H. M. Chun, C. E. Padilla, D. N. Chin, M. Watanabe, V. I. Karlov, H. E. Alper, K. Soosaar, K. B. Blair, O. M. Becker, L. S. D. Caves, R. Nagle, D. N. Haney, and Barry L. Farmer, J. Comput. Chem. 2000, 21, 159.
- [16] G. Balaraman, I. Park, A. Jain, N. Vaidehi, J. Phys. Chem. B 2011, 115, 7588.
- [17] I. Park, V. Gangupomu, J. Wagner, A. Jain, N. Vaidehi, J. Phys. Chem. B 2012, 116, 2365.
- [18] V. Gangupomu, I. -H. Park, J. Wagner, A. Jain, N. Vaidehi, submitted, 2012.
- [19] A. Jain, I. -H. Park, N. Vaidehi, J. Chem. Theory Comput. 2012, 8, 2581.
- [20] W. C. Still, A. Tempczyk, R. C. Hawley, T. Hendrickson, J. Am. Chem. Soc. 1990, 112, 6127.
- [21] M. Fixman, Proc. Natl. Acad. Sci. USA 1974, 71, 3050.
- [22] Y. Sugita, Y. Okamoto, Chem. Phys. Lett. 1999, 314, 141.
- [23] N. Metropolis, A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller, E. Teller, J. Comput. Phys. 1953, 21, 1087.
- [24] S. Plimpton, J. Comput. Phys. 1995, 117, 1.
- [25] M. S. Friedrichs, P. Eastman, V. Vaidyanathan, M. Houston, S. LeGrand, A. L. Beberg, D. L. Ensign, C. M. Bruns, V. S. Pande, J. Comput. Chem. 2009, 30, 864.
- [26] A. Mathiowetz, A. Jain, N. Karasawa, W. Goddard, Proteins: Struct. Funct. Bioinform. 1994, 20, 227.
- [27] A. Jain, Robot and Multibody Dynamics: Analysis and Algorithms, 1st ed.; Springer, 2010. p. 35.
- [28] S. Harvey, R. Tan, T. Cheatham, J. Comput. Chem. 1998, 19, 726.
- [29] M. Poursina, K. Bhalerao, S. Flores, K. Anderson, A. Laederach, Methods Enzymol. 2011, 487, 73.
- [30] D. Frishman, P. Argos, Proteins 1995, 23, 566.
- [31] A. Jain, G. Rodriguez, J. Guid. Control Dynam. 1993, 16, 830.
- [32] D. Jones, J. Mol. Biol. 1999, 292, 195.
- [33] S. D. Cohen, A. C. Hindmarsh, Comput. Phys. 1996, 10, 138.
- [34] S. B. Ozkan, G. A. Wu, J. D. Chodera, K. A. Dill, Proc. Natl. Acad. Sci. USA 2007, 104, 11987.
- [35] G. J. Martyna, M. L. Klein, M. Tuckerman, J. Chem. Phys. 1992, 97, 2635.
- [36] S. Bond, B. Leimkuhler, B. Laird, J. Comput. Phys. 1999, 151, 114.
- [37] M. Fixman, J. Comput. Phys. 1978, 69, 1527.
- [38] N. Go, H. Scheraga, Macromolecules 1976, 9, 535.
- [39] A. Jain, J. Comput. Phys. 1997, 136, 289.
- [40] V. Katritch, M. Totrov, R. Abagyan, J. Comput. Chem. 2003, 24, 254.

Received: 29 August 2012 Revised: 26 October 2012 Accepted: 4 November 2012 Published online on 23 January 2013

![](_page_10_Picture_49.jpeg)