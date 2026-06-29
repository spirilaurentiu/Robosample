# **Application of Torsion Angle Molecular Dynamics for Efficient Sampling of Protein Conformations**

#### **JIANHAN CHEN, WONPIL IM, CHARLES L. BROOKS, III**

*Department of Molecular Biology and Center for Theoretical Biological Physics, The Scripps Research Institute, 10550 North Torrey Pines Road, La Jolla, California 92037*

> *Received 5 July 2005; Accepted 8 July 2005 DOI 10.1002/jcc.20293 Published online in Wiley InterScience (www.interscience.wiley.com).*

**Abstract:** We investigate the application of torsion angle molecular dynamics (TAMD) to augment conformational sampling of peptides and proteins. Interesting conformational changes in proteins mainly involve torsional degrees of freedom. Carrying out molecular dynamics in torsion space does not only explicitly sample the most relevant degrees of freedom, but also allows larger integration time steps with elimination of the bond and angle degrees of freedom. However, the covalent geometry needs to be fixed during internal coordinate dynamics, which can introduce severe distortions to the underlying potential surface in the extensively parameterized modern Cartesian-based protein force fields. A "projection" approach (Katritch et al. J Comput Chem 2003, 24, 254 –265) is extended to construct an accurate internal coordinate force field (ICFF) from a source Cartesian force field. Torsion crossterm corrections constructed from local molecular fragments, together with softened van der Waals and electrostatic interactions, are used to recover the potential surface and incorporate implicit bond and angle flexibility. MD simulations of dipeptide models demonstrate that full flexibility in both the backbone / and side chain <sup>1</sup> angles are virtually restored. The efficacy of TAMD in enhancing conformational sampling is then further examined by folding simulations of small peptides and refinement experiments of protein NMR structures. The results show that an increase of several fold in conformational sampling efficiency can be reliably achieved. The current study also reveals some complicated intrinsic properties of internal coordinate dynamics, beyond energy conservation, that can limit the maximum size of the integration time step and thus the achievable gain in sampling efficiency.

© 2005 Wiley Periodicals, Inc. J Comput Chem 26: 1565–1578, 2005

**Key words:** force field; generalized Born; implicit solvent; internal coordinate; NMR refinement; replica exchange; structure prediction

## **Introduction**

Sampling is one of the major bottlenecks in molecular dynamics (MD) studies of many problems in chemistry and biology. In particular, the time scale of traditional MD simulations in Cartesian coordinates is severely limited by the femtosecond (fs) integration time steps required by high-frequency degrees of freedom (bonds and angles). It is clear that these hard degrees of freedom are only mildly excited in dynamical processes at room temperature, and have little impact on slow processes such as conformational transitions in macromolecules.1–3 A natural and popular approach to increase the time step in conventional (Cartesian) MD, and thus extend the accessible simulation time scale, is to eliminate the hardest degrees of freedom. For example, the SHAKE algorithm4 has become a standard approach for constrained MD with fixed bond lengths and/or angle, where the holonomic constraints are handled by iterative corrections of violations at each integration step. However, the SHAKE algorithm cannot handle networks of holonomic constraints efficiently.5 As such, SHAKE is typically applied only to fix lengths of bonds that involve hydrogen atoms, which generally allows an increase of time step from 1 fs to 2 fs. Various alternative approaches have also been developed to extend

*Correspondence to:* C. L. Brooks, III; e-mail: brooks@scripps.edu Contract/grant sponsor: La Jolla Interfaces in Science (to C.J.H.)

Contract/grant sponsor: National Institutes of Health; contract/grant numbers: RR12255 and GM48807

Contract/grant sponsor: National Science Foundation; contract/grant numbers: PHY0216576 and PHY0225630

This article includes Supplementary Material available from the authors upon request or via the Internet at http://www.interscience.wiley.com/ jpages/0192-8651/suppmat

![](_page_1_Figure_2.jpeg)

**Figure 1.** Backbone  $\phi/\psi$  potential energy surfaces for the alanine dipeptide (Ace-Ala-NMe) in vacuum using the CHARMM PARAM22 force field with (a) fixed ideal covalent geometry and (b) fully relaxed covalent geometry (by energy minimization in Cartesian space). The angles are in degrees and energy values in kcal/mol

the effective time steps such as multiple time scale methods<sup>6,7</sup> and related Langevin stabilization of molecular dynamics.<sup>8</sup>

A particularly attractive alternative is to solve the equations of motion directly in torsion space. This approach completely eliminates bond and angle degrees of freedom and, in principle, allows much larger time steps. The equations of motion for such an internal variable system can be expressed as follows:<sup>2</sup>

$$\Omega(\theta)\ddot{\theta} + C(\theta, \dot{\theta}) = T(\theta), \tag{1}$$

where  $\theta$  denotes an *n*-dimensional vector of internal coordinates (torsion angles):  $\Omega(\theta)$  is an  $n \times n$  mass matrix:  $C(\theta, \dot{\theta})$  is an *n*-dimensional vector of Coriolis forces; and  $T(\theta)$  is an *n*-dimensional vector of generalized force. Note that the number of torsional degrees of freedom, n, is often much smaller than the 3N – 6 Cartesian internal degrees of freedom (N being the total atom number). Computation of the  $\Omega$ , C, and T matrices is nontrivial but possible.<sup>2,3</sup> However, because the mass matrix is no longer diagonal, the computational cost of solving eq. (1) directly scales as  $O(n^3)$ . The sharp increase of computational cost for large systems had severely limited the applicability of the internal variable approach, until the development of several efficient algorithms with computational efforts that scale linearly with respect to the system size. 10-14 In particular, Jain et al. 13 described a fast recursive algorithm based on Newton–Euler inverse mass operator (NEIMO) for evaluating the  $\Omega$ , C, and T matrices and solving the equations of motion, which is suitable for torsion angle molecular dynamics (TAMD) simulations of macromolecules using traditional Cartesian-based force fields. A comprehensive review of internal coordinate simulation methods can be found elsewhere. 15 TAMD has been successfully applied to study the equilibrium properties of very large systems, 16 to model nucleic acids, 17 and to refine NMR and X-ray structures. 18-21 Simulated annealing in torsion space has been especially successful in NMR structure calculation and refinement, and now becomes a standard technique that demonstrates significant enhancement in conformational sampling efficiency compared to previous protocols based on Cartesian MD (CMD). Note that simplified force fields are typically employed in these calculations, for example, ignoring electrostatic

interactions and employing a purely repulsive soft sphere interaction. Another interesting recent application is coupling TAMD with global optimization techniques for molecular structure prediction,<sup>22</sup> where both a simplified force field as well as realistic ECEPP/3 torsion-based force field<sup>23</sup> were employed. It is reasonable to expect that TAMD can be used to substantially enhance conformational sampling even with more realistic modern all-atom force fields. However, its initial application to general purpose MD<sup>16,17,24–27</sup> has only achieved limited success, and little further development in this direction has, to the best of our knowledge, been seen in the literature. The main reason is probably linked to several consequences associated with the covalent rigid body approximation, discussed below.

Significant improvements in physics-based force field development and parameterization have been seen during the last twenty years, 28,29 especially in the widely used general purpose force fields including CHARMM, 30,31 Amber, 32 OPLS, 33 and GRO-MOS.34 Particular advances have been recently made in the development of efficient generalized Born (GB) implicit solvent models. 35-37 These empirical force fields have been very successful in various application to biological systems; active developments are still ongoing, and continued advances are expected.<sup>29</sup> Therefore, it is highly desirable to exploit the potential of TAMD in enhancing conformational sampling of proteins in the context of these Cartesian-based force fields, which have been parameterized with flexible covalent geometry. Unfortunately, several wellknown difficulties exist for such an application. The primary obstacle arises as a consequence of fixing the covalent geometry during internal coordinate dynamics. 38-40 First, van der Waals repulsions between 1-4 atoms (atoms separated by three bonds) are not relaxed in the (fixed) covalent geometry, leading to higher rotational barriers. Second, severe clashes can occur between atoms separated by more than three bonds, leading to prohibitively high barriers on the internal coordinate potential energy surface. For example, Figure 1 illustrates the differences in the backbone  $\phi/\psi$  potential energy surfaces of the alanine dipeptide with fixed and fully relaxed covalent geometries. Note that, when the covalent geometry is fixed, not only are the energy minima distorted but the energy barriers in many regions are also dramatically elevated

due to long-range atom clashes. The barrier increase around  $(\phi, \psi) = (0^{\circ}, 0^{\circ})$  is attributed to a clash between the N-terminal carbonyl oxygen and C-terminal amide hydrogen (separated by six bonds). The barrier increase around  $\phi = 120^{\circ}$  is due to the clash between the N-terminal carbonyl oxygen and side chain methyl group, separated by four or five bonds. Note that the surface in Figure 1a was computed with rigid peptide planes. Corresponding surface with flexible peptide planes show less increase in barrier heights but the overall distortions are qualitatively similar (data not shown). Finally, identifying the optimal fixed covalent geometry is often problematic. Such a geometry not only needs to comply with the force field parameters<sup>39</sup> but also depends on the specific problems to be addressed.

Recently, Katritch et al. 38 described a "projection" approach to construct an internal coordinate force field (ICFF) in close approximation to the source Cartesian force field. Using six-fold Fourier series torsion energy correction terms (to correct 1–4 repulsion) and a soft polynomial repulsion function (to relieve 1-5 and 1-6 clashes), some degree of implicit flexibility is restored in the rigid body approximation. The fixed covalent geometry is determined by Cartesian energy minimization using a local part of the Cartesian force field. Using these techniques, Katritch et al. were able to reproduce the torsion potential energy profile around the energy minima to a satisfactory accuracy and successfully applied the resulting ICFF to evaluate the conformation energies of a set of organic molecules. However, to study large conformational changes such as protein folding using MD (instead of Monte Carlo) simulations, one does not only need to reproduce the potential surface near energy minima but also to recover the correct energy maxima (barriers). Furthermore, as illustrated in Figure 1, increase of energy barriers on the  $\phi/\psi$  potential surface with fixed covalent geometry is mainly due to clashes between atoms that are separated by four bonds or more, and such increase depends on more than two or more torsion angles. Therefore, the original projection approach cannot be directly applied to construct an ICFF that is suitable for TAMD simulations of large conformational changes for proteins. Furthermore, although the fixed covalent geometry determined by local Cartesian energy minimization might be appropriate for accurate evaluation of conformational energies of small molecules near their optimal conformational states, it is not likely to be suitable for exploring large conformation changes where the fixed covalent geometry needs to be representative for a broad range of conformational states. In the current study, by extending the basic ideas of the projection approach, we describe an improved method of constructing an ICFF with implicit flexibility that allows effective use of TAMD to enhance conformational sampling of proteins with modern Cartesian-based force fields. The new approach is designed specifically for modeling proteins, but the basic ideas can be extended to modeling other macromolecules. The quality of resulted ICFF is carefully examined by comparing the underlying backbone  $\phi/\psi$ and side chain  $\chi_1$  potential surfaces as well as the dynamics properties of dipeptide models. Then the efficacy as well as limitations of TAMD in augmenting conformational sampling will be examined through folding simulations of peptides and refinement of protein NMR structures.

#### Methods

The recursive NEIMO algorithm that we employ here to solve the equations of motion in torsion space has been described in detail elsewhere <sup>13</sup> and will be first summarized here to illustrate the essential features. Following, an effective approach of constructing an ICFF with implicit flexibility from a Cartesian force field will be described.

## Recursive Solution of Equations of Motion

The efficient recursive NEIMO algorithm of Jain et al. 13 has been implemented in a CHARMM TAMD module. Forces arising from the Cartesian force field can be used directly to propagate the time development of the system in internal coordinates, and evaluation of the gradients of the Cartesian potential function with respect to internal variables is not necessary. The algorithm starts by representing the molecule with a branched tree structure. The tree consists of rigid bodies of atoms with invariable relative positions, clusters, connected by hinges where the permissible relative motion between adjoined clusters can be partially constrained. The internal coordinates are associated with the hinges, and describe the relative position and motion of clusters. In the current implementation, we are mainly interested in hinges that allow only a single torsional degree of freedom. One particular cluster in the tree is designated the base cluster, whose choice has only slight impact on the computational efficiency. Adjoining clusters are referred to as parent/child, with the one on the path to the base being the parent and the other being the child. The tip clusters have no children. For each cluster there is one origin atom, where the hinge connection to its parent cluster originates. Additionally, nontip clusters have one or more branching atoms where hinges to their child clusters originate. Efficient expressions for evaluating the relevant matrices and solving the equations of motion for such a tree-topology molecule can be derived using spatial operator algebra. 12,13 Structures forming closed loops such as disulfide bonds can also be handled efficiently by extending the algorithm described below.41

To illustrate these principles, consider the simplest case of a serial chain of clusters connected by hinges that only allow one rotational degree of freedom, as shown in Figure 2. The conformation of such a serial chain is uniquely defined by the position and orientation of the base cluster, designated as the (n + 1)th cluster, and n torsion angles through recursive relationships. The velocity, acceleration and hinge force of clusters satisfy the Newton–Euler recursive relations. <sup>13</sup>

$$V_{k} = \phi_{k+1,k}^{T} V_{k+1} + H_{k}^{T} \dot{\theta}_{k}, \tag{2}$$

$$\alpha_k = \phi_{k+1,k}^T \alpha_{k+1} + H_k^T \ddot{\theta}_k + a_k, \tag{3}$$

$$F_k = \phi_{k,k-1} F_{k-1} + M_k \alpha_k + b_k - F_k^{(c)}, \tag{4}$$

$$T_k = H_k F_k, (5)$$

where the spatial velocity  $V_k$  is a six-dimensional vector consisting of the angular and linear velocity vectors,  $V_k = \text{Col}[\omega_k, v_k]$ ;  $\phi_{x,y}$ 

![](_page_3_Figure_2.jpeg)

**Figure 2.** A serial chain of clusters connected by hinges. Note that the kth hinge connects the (k + 1)th and kth clusters.  $q_{k,0}$  denotes the position of the origin atom of the kth cluster and  $q_{k,b}$  that of the branching atom.

is a  $6 \times 6$  spatial transformation matrix that transforms quantities between two frames with origins x and y; hinge matrix  $H_k^T =$  $\operatorname{Col}[\hat{h}_k, 0, 0, 0]$  with  $\hat{h}_k$  being the unit vector in the direction of the kth hinge; the spatial acceleration  $\alpha_k$  is the time derivative of the spatial velocity;  $F_k^{(c)}$  is a six-dimensional effective Cartesian spatial force computed by combining the atomic spatial force of all the atoms of the kth cluster. The atomic spatial force consists of the torque (with respect to the cluster origin) and force as given by the potential;  $F_k$  is the hinge spatial force of interactions between the (k+1)th and kth clusters;  $M_k$  is the spatial inertia with respect to the cluster origin; the spatial gyroscopic force,  $a_k$ , and Coriolis acceleration term,  $b_k$ , are functions of spatial velocities, inertia and hinge velocity;  $T_k$  is the projection of the hinge force along the allowable degree(s) of freedom of the kth hinge, whose value vanishes when all forces from the force field are included in  $F_{\nu}^{(c)}$ . The superscript T denotes matrix transpose.

For a serial chain with n + 1 clusters, n such recursive relations as in eqs. (2)–(5) exist. Spatial operators can be used to express the recursive relations concisely.<sup>12</sup> For example, eq. (2) can be summarized as.

$$V = \Phi^T H^T \dot{\theta}, \tag{6}$$

in which the lower triangular  $6n \times 6n$  spatial operator  $\Phi$  is constructed from  $\phi_{k+1,k}$ ; V and  $\dot{\theta}$  are  $6n \times 1$  and  $n \times 1$  vectors constructed from  $nv_k$  vectors and  $\theta_k$ ; and operator  $H^T$  is a  $6n \times n$  block diagonal matrix,  $H^T = \text{diag}\{H_1^T, \ldots, H_n^T\}$ . Equations (3)–(5) can be rewritten in a similar fashion. Finally, by substituting the spatial operator forms of eqs. (2)–(4) into that of eq. (5), we arrive at the equations of motion in internal coordinates of eq. (1) with

$$\Omega(\theta) \equiv H \Phi M \Phi^T H^T, \tag{7}$$

$$C(\theta, \dot{\theta}) \equiv H\Phi(M\Phi^{T}a + b - F^{(c)}). \tag{8}$$

The mass matrix  $\Omega(\theta)$  has a dimensionality of  $n \times n$  and its factorization in eq. (7) is called the Newton-Euler operator fac-

torization. <sup>14</sup>  $C(\theta, \dot{\theta})$  is a  $n \times 1$  vector that contains contributions of Coriolis, centrifugal, gyroscopic, and Cartesian forces.

Solving the equations of motion for the acceleration vector requires the inversion of the mass matrix. An efficient recursive algorithm for its inverse is possible through an alternative operator factorization of the mass matrix, named the *innovations operator factorization*. <sup>13</sup>

$$\Omega(\theta) = [1 + H\Phi K]D[1 + H\Phi K]^{T}, \tag{9}$$

where the factors are square and invertible such that

$$\Omega^{-1} = [1 - H\Phi K]^T D^{-1} [1 - H\Phi K]^T.$$
 (10)

Note that matrix D is positive definite and diagonal when all hinges allow only a single degree of freedom. Equation (10) represents a recursive algorithm for inverting the mass matrix with O(n) computational cost. In addition, explicit computation of the mass matrix is not required as all the auxiliary matrices can be evaluated directly. The extension of the above formalism to a true tree with branches is straightforward. During the recursion from tips to base, at each cluster, quantities from each of its children are summed up before proceeding; during a base to tips recursion, the recursions continue separately along each of the branches of each cluster. The detailed computational steps involved in solving the equations of motion can be found in the original reference.<sup>13</sup>

The base cluster of a tree is connected to the inert frame (or another tree) by a virtual hinge that allows all six translational and rotational degrees of freedom. Such a hinge requires special treatment in updating the base orientation during the dynamics. The hinge matrix for six degrees of freedom hinges is simply a  $6 \times 6$  unit matrix. One can use the differentials of quasi-coordinates to express the associated  $6 \times 1$  generalized velocity vector,  $\dot{\theta} = V_k = \text{Col}[\omega_k, v_k]$ . Here, the quasi-coordinates refer to the integrals of the angular velocity, which do not have physical meanings and only their time derivatives have physical meanings. Generalized coordinates such as the Euler angles or a quaternion representation can be used to describe the orientation of the base cluster. The quaternions of a rigid body satisfy the quaternion equations of motion, Which can be solved iteratively using Fincham's implicit quaternion algorithm.

## Implementation in CHARMM

The TAMD module of CHARMM contains three major facilities: tree-topology setup, energy minimization, and molecular dynamics. The required tree representation of the system can be set up automatically for proteins with CHARMM standard residues, or semiautomatically for nonstandard residues. The internal coordinates are defined automatically to be consistent with the tree representation. More detailed description of the interface will be included in the forthcoming CHARMM releases. The molecular dynamics is implemented using a modified leapfrog formulation of the Verlet algorithm, <sup>16,43</sup> which can be considered as a simplification of the more rigorous implicit leap-frog integrator. <sup>15</sup> As demonstrated by the energy conservation properties, such simplification seems to be sufficient (see Results and Discussion). At the

beginning of the dynamics, velocities are first assigned from a Gaussian distribution in the Cartesian space and then "mapped" onto the internal coordinates. 21 A simple Berendsen thermostat 45 in combination with velocity reassignment is used for temperature control. Although more sophisticated thermostats such as Nosé-Hoover chains 46,47 can be implemented within the recursive TAMD algorithm, <sup>24</sup> it does not appear that they will significantly impact on the dynamic properties of the system, particularly in the context of studying large conformational changes in combination with advanced sampling techniques such as the replica exchange method (REX).48 Validation of the complex dynamics algorithm was based on examining the energy conservation properties of constant energy simulations as well as dynamical properties of peptides of the constant temperature simulations and REX MD simulations (see Results and Discussion section). Note that increase of the computational cost with TAMD is minimal, and the energy and force evaluations are still the limiting computational steps. In fact, none of the bond and angle related internal energy terms need to be evaluated in TAMD, which slightly reduces the overall computational cost of a single dynamics step propagation. Energy minimization does not involve the mass matrix. However, the expressions that are derived above can be used to evaluate the gradients with respect to the internal coordinates efficiently from the Cartesian forces.<sup>21</sup> Once the gradients in internal coordinates are available, standard minimization algorithms such as steepest descent and Powell's Methods<sup>49</sup> can be used.

## Towards an ICFF for Proteins with Implicit Flexibility

As discussed previously, several obstacles exist in using TAMD directly with Cartesian-based force fields, which are key to the efficacy of TAMD in enhancing the conformational sampling. One of the primary focuses of this work is to effectively address these problems. In this section, we describe an extension of the projection idea of Katritch and coworkers<sup>38</sup> to construct an ICFF that is in accurate agreement with the original Cartesian force field and suitable for MD simulation of large conformational changes of proteins.

# Local Fragments and Dihedral Crossterm Corrections

Both minima and maxima of the potential energy surface in the projected ICFF need to be sufficiently reproduced to study large conformational changes using dynamics simulations techniques. As illustrated in Figure 1, clashes between atoms separated by four or more bonds can lead to a dramatic barrier increase on the potential surface and such increase depends on more than one torsion angle. To correct these, we first use softened van der Waals and electrostatic interactions, where the interactions are not altered until the energies exceed certain thresholds.30 The interactions have linear forms when the soft core potential is active. Such softening introduces minimal disturbance to the interactions normally seen in proteins and peptides (e.g., the backbone  $\phi/\psi$  potential surfaces are essentially identical in Cartesian coordinates with and without softening), while effectively reducing the barrier increase to manageable range (e.g., see Fig. 4a). Then, the residual increase of barriers is corrected by dihedral crossterms using the CMAP facility in CHARMM.<sup>50–52</sup> These crossterm corrections

![](_page_4_Picture_7.jpeg)

**Figure 3.** Local fragments used to construct dihedral crossterm correction maps for proteins. Atoms circled by the dashed line are included in calculation of the backbone  $\phi/\psi$  correction map. Additional  $\phi/\chi_1$  and  $\psi/\chi_1$  maps are computed for all possible combinations of X and Y atom types.

are computed from local torsion crossterm fragments, described in details below.

The fragment used for computing the backbone  $\phi/\psi$  crossterm correction consists of —CO—NH—C<sub>o</sub>H—NH—CO—, shown in Figure 3. Two-dimensional grid-based correction maps are constructed by computing the difference between two potential energy surfaces with fully relaxed and fixed covalent geometries. The relaxation of covalent geometry is achieved by extensive restrained energy minimization in Cartesian space. A single  $\phi/\psi$ correction map is shared by all residues except proline. Note the use of such dihedral crossterms corrections requires fixed peptide planes. The consequences of rigid peptide planes are not completely clear. However, simulations of dipeptide models, simple peptides, and small proteins appear to indicate that the protein dynamics is not significantly altered (see Results and Discussion section). The effect of  $\phi/\psi$  correction is illustrated using valine dipeptide, shown in Figure 4b. Distortion to the potential surface around  $(\phi, \psi) = (0^{\circ}, 0^{\circ})$  is effectively reduced.

Clashes between backbone and side chain atoms also lead to barrier increase on the potential surface, and these are corrected by two additional crossterm corrections,  $\phi/\chi_1$  and  $\psi/\chi_1$  maps. The standard residues are divided into several groups based on the  $C_\beta$  connectivity (see Fig. 3). For example, Val and Ile share the same set of  $\phi/\chi_1$  and  $\psi/\chi_1$  correction maps. Common  $\phi/\chi_1$  and  $\psi/\chi_1$  correction maps. Common  $\phi/\chi_1$  and  $\psi/\chi_1$  correction maps are then constructed for each group using appropriate local fragments, which are —CO—NH—C $_\alpha$ —C $_\beta$ H(X—)(Y—) and C $_\beta$ H(X—)(Y—)—C $_\alpha$ —NH—CO—, respectively. All interactions are included in the correction map construction. Atomic charges are neutralized when bonded atoms are deleted from the local fragment. For example, charges on  $C_\gamma$  atoms of Valine residues need to be set to 0.0 (from -0.27e) as

![](_page_5_Figure_2.jpeg)

**Figure 4.** Backbone  $\phi/\psi$  potential surfaces for the valine dipeptide in vacuum with (a–c) fixed, and (d) fully relaxed covalent geometry. Softened van der Waals and electrostatic interactions were used in (a–c). When all three crossterm corrections  $(\phi/\psi, \phi/\chi_1, \text{ and } \psi/\chi_1)$  are applied, shown in (c), the potential energy surface is quite accurately recovered over the entire space.

three bonded hydrogens are deleted from the fragment. Similarly, charges on  $O_{\gamma}$  of Serine and Threonine residues are set to be  $-0.23\mathrm{e}$  from the original values of  $-0.66\mathrm{e}$  after deleting the bonded hydrogens with charges of 0.43e. The effects  $\phi/\chi_1$  and  $\psi/\chi_1$  corrections on the backbone  $\phi/\psi$  surface as well as sidechain  $\chi_1$  torsion profiles are demonstrated in Figure 4c and Figure 5.

Combination of the three correction maps described above can effectively reproduce the backbone  $\phi/\psi$  and side-chain  $\chi_1$  torsional potential surfaces of the original Cartesian force fields. Note that the corrections are applied on residue level and no additional tuning is necessary for individual sequence. The need to read in larger parameter files that contain the precalculated correction maps can lead to a noticeable increase (typically a few seconds) in the overhead of starting a new CHARMM job. Computational cost associated with softening the van der Waals and electrostatic interactions and evaluating additional correction terms is comparable to those of evaluating bond and angle terms and is thus minimal compared to the speed limiting steps of evaluating nonbonded interactions (especially with GB implicit solvent).

The torsion barrier increase due to 1-4 van der Waals repulsion can affect the flexibility of side chains (except  $\chi_1$ ), which might be corrected by the original projection approach.<sup>38</sup> However, such effects appear to be minor for protein side chains and the side-chain conformational equilibrium and dynamics flexibility appear to be minimally impacted. This is demonstrated by 100 nanosecond (ns) room temperature MD simulations of dipeptide models.

Side-chain conformation properties of two representative residues, Lys and Trp, are shown in the supporting materials. Furthermore, small changes of side-chain flexibility is not expected to have a

![](_page_5_Figure_8.jpeg)

**Figure 5.** Comparison of side-chain  $\chi_1$  torsion energy profiles for the alanine dipeptide at two representative backbone configurations. The solid lines were obtained with full energy minimization in Cartesian space. The dashed and dotted lines were obtained in torsion space with and without the crossterm corrections.

large impact on proteins, especially in protein folding and unfolding studies. Therefore, at the current stage, we have mainly focused on obtaining correct backbone and side-chain  $\chi_1$  properties, and no further correction has been applied for other side-chain torsions

#### Fixed Covalent Geometry

The optimal fixed covalent geometry might depend on the particular application and its identification is not always obvious. The fixed covalent geometry optimized through local energy minimization in Cartesian space<sup>38</sup> is not representative for molecular structures undergoing large conformational changes. In particular, specific backbone covalent geometries are required for forming compact secondary structures like an  $\alpha$ -helix. Fixing the covalent geometry at different configurations can severely hinder the ability to form such compact structures. As such, we currently optimize the fixed backbone geometry to allow  $\alpha$ -helix formation. Preliminary studies show that this geometry appears to have minimal impact on the ability to form  $\beta$ -strands and the stability of native  $\beta$ -structures (data not shown). However, more investigation is necessary to fully assess the consequences of fixing the covalent geometry at any configuration in the context of protein folding studies.

#### Simulation Protocols

All simulations were carried out using CHARMM30 with PARAM22 all-atom force field.<sup>31</sup> Room-temperature dynamics simulations of dipeptide models were carried out to examine the backbone and side-chain conformation equilibrium and dynamics properties in Cartesian and internal coordinates. The initial structures were first built in CHARMM and then subject to 100-ns constant temperature (NVT) simulations. Nosé-Hoover thermostat<sup>46,47</sup> was used in CMD, and Berendsen thermostat<sup>45</sup> in TAMD. Constant energy (NVE) simulations were also carried out to study energy conservation properties (namely, fluctuation and drift) using several systems with varying sizes and topologies: pentapeptide Met-enkephalin (Ace-Tyr-Gly-Gly-Phe-Met-NMe) in a collapsed coil conformation; Polyvaline (Val)<sub>10</sub> in canonical  $\alpha$ -helical conformation; 56-residue B1 domain of protein G (GB1) in its native  $\alpha/\beta$  fold<sup>53</sup> and 159-residue Dihydrofolate reductase  $(DHFR)^{54}$  in native  $\alpha/\beta$  fold with several long loops. Initial structures (experimental structures or models built in CHARMM) were first equilibrated at 300 K before NVE production simulations. All these simulations were carried out in vacuum with a constant dielectric constant of 1.0. As the native conformations of GB1 and DHFR are not stable in vacuum and weak harmonic restraints were applied to backbone heavy atoms (0.1 kcal/mol/Å<sup>2</sup>) to prevent dramatic conformational changes. No cutoff was used for the nonbonded interactions unless otherwise specified. SHAKE<sup>4</sup> was applied to fix the length of bonds involving hydrogens in all CMD simulations.

The peptide folding and NMR protein structure refinement REX MD simulations were enabled by the MMTSB Tool Set (available from http://mmtsb.scripps.edu).<sup>55,56</sup> The GBSW implicit solvent model<sup>57</sup> was used in the refinement simulation and a related implicit membrane model<sup>58</sup> was used in the folding simu-

lation of the WALP transmembrane peptides.<sup>59,60</sup> In replica exchange simulations, multiple copies (replicas) of the system are simulated at different temperatures independently and simultaneously. The temperatures are usually distributed exponentially within a specified range, and there is always one single replica simulated at each temperature. Replicas attempt to exchange simulation temperatures according to a Metropolis type algorithm after a number of steps of simulation. In the course of an REX simulation, replicas can travel up and down the temperature space automatically in a self-regularized fashion, which in turn, induces a nontrivial walk in temperature space. REX can greatly reduce the probability of being trapped in states of local energy minima and sample a larger conformation space.

#### **Results and Discussion**

In this section, the quality of the ICFF constructed using the extended projection approach is first examined by direct comparison of the backbone  $\phi/\psi$  and side-chain  $\chi_1$  torsion energy profiles with those in the original Cartesian force field. Then backbone and side-chain conformational equilibrium and dynamics properties of representative dipeptide models are examined to check whether the implicit flexibility is sufficiently and accurately recovered in the ICFF. Finally, folding of simple peptides and refinement of NMR structures of a protein domain are carried out to study the efficacy of TAMD in enhancing conformational sampling. We will also discuss potential limiting factors and remaining problems in the application of TAMD to augment conformational sampling.

#### Quality of $\phi/\psi$ and $\chi_I$ Potential Surfaces in ICFF

Using dipeptide models, we have examined the recovery of Cartesian backbone  $\phi/\psi$  and side-chain  $\chi_1$  torsion energy profiles for all residues (except proline) in the ICFF. As illustrated using the valine dipeptide, shown in Figure 4, the approach described in the previous section can accurately reconstruct the whole potential surface including both energy minima and maxima. The dramatic increase in energy barriers is first significantly reduced by the use of softened van der Waals and electrostatic interactions (see Fig. 4a). The remaining distortions of the energy surface are then effectively removed by the crossterm corrections, demonstrated by Figure 4b-c. Note that the  $\phi/\psi$  crossterm correction is only responsible for reducing barriers around  $\phi \sim 0^{\circ}$ , while the  $\phi/\chi_1$  and  $\psi/\chi_1$  terms are responsible for regions around  $\phi = 120^{\circ}$ . The root-mean-square difference (RMSD) between the final ICFF map and original map is 0.67 kcal/mol for the low energy regions (within 5 kcal/mol from the global minimum) and 1.6 kcal/mol for the whole surface. Similar results were obtained for all the standard residue types (except proline).

Flexibility of the side chain along  $\chi_1$  is also important in large conformational reorganizations. The  $\chi_1$  torsion profile explicitly depends on the backbone  $\phi$ ,  $\psi$  conformation, rendering straightforward corrections<sup>38</sup> completely inapplicable. However, as demonstrated in Figure 5, the ICFF  $\chi_1$  torsion energy profiles are essentially identical to the original ones for all possible backbone conformations. These results indicate that implicit flexibility can

![](_page_7_Figure_2.jpeg)

**Figure 6.** Backbone / distribution for the valine dipeptide during 100 ns CMD and TAMD simulations at 300 and 1000 K. All simulations were carried out in vacuum using the CHARMM PARAM22 force field or projected ICFF. No cutoff was used for nonbonded interactions.

be successfully restored in the rigid covalent geometry approximation by the extended projection approach.

## *Conformational Equilibrium and Dynamics Properties of Dipeptide Models*

Room-temperature MD simulations of 100 ns have been carried out for representative dipeptide models in both torsion and Cartesian space to examine the backbone and side-chain conformational equilibrium and dynamics. The hydrogen mass was increased to be 6.0 in all TAMD simulations (see next section). Both / and <sup>1</sup> distribution and transition properties were examined. The results show that conformational equilibrium as well as flexibility are very well restored in the ICFF for all the dipeptide models with the exception of proline, where flexibility along and <sup>1</sup> is completely lost. For example, Figures 6 and 7 compare the dynamic properties of the valine dipeptide. They demonstrate that the dipeptide in the ICFF is at least as flexible as in the original Cartesian force field. In fact, more transitions between states of energy minima were observed in TAMD simulations, especially at 300 K. This is probably an indication of enhanced conformational sampling due to exclusive sampling of the torsional degrees of freedom in TAMD. However, definitive conclusions are difficult to be made as the underlying energy surfaces in CMD and TAMD, while being highly similar, still differ slightly. Although no correction is applied to other side-chain torsions, dynamics simulations of dipeptide models indicate that the conformational and dynamics properties are only minimally impacted (see supporting materials). In fact, more transitions between rotamer states were consistently observed in TAMD simulations, despite the slightly elevated rotational barrier due to lack of correction for 1– 4 van der Waals repulsion.

#### *Energy Conservation Properties*

Once an ICFF in close approximation to the source Cartesian force field is available, one can start to investigate the efficacy of TAMD in conformational sampling enhancement. One of the main advantages of carrying out molecular dynamics in torsion space is the possibility to use larger integration time steps. The first property that needs to be examined is the accuracy of the dynamics, which can be measured by the fluctuation in the total energy during a constant energy (NVE) simulation.16,43 The energy fluctuation can be computed as

$$\delta E = \frac{\langle E^2 \rangle - \langle E \rangle^2}{\langle E_k \rangle},\tag{11}$$

where *E* and *Ek* are the total energy and total kinetic energy. Figure 8 shows the energy conservation as a function of integration time step in CMD and TAMD of four systems ranging from the pentapeptide Met-enkephalin to the 159-residue enzyme DHFR. Several observations can be made from these results. First, as shown in Figure 8a, TAMD conserves the total energy much better than CMD. The energy fluctuations in TAMD are smaller by about two to three orders of magnitude in all cases. With SHAKE applied

![](_page_7_Figure_12.jpeg)

**Figure 7.** Side-chain <sup>1</sup> angle transition and distribution properties for the valine dipeptide during 100 ns CMD (black traces) and TAMD (red traces) simulations at 300 and 1000 K. Note that only results from the first few nanoseconds are shown in (a) and (b). In (c), the solid traces were computed from simulations at 300 K and dashed traces from those at 1000 K.

![](_page_8_Figure_2.jpeg)

**Figure 8.** Energy fluctuations during CMD and TAMD simulations as functions of time steps for Met-Enk, (Val)10, GB1, and DHFR. The fluctuations were computed from over 1.0 ps NVE simulations.

to fix the length of hydrogen attached bonds, CMD simulations conserve energy reasonably well with time steps up to 3 fs and diverge quickly with time steps over 5 fs. TAMD conserves energy well with time steps up to 10 fs for small peptides and up to 5 fs for larger proteins, and does not diverge even with time steps as large as 20 fs. Second, the energy fluctuation of CMD and TAMD simulations show very different dependence on the size and degree of compactness of the system (see below). Third, simply increasing the mass of all hydrogen atoms in TAMD can significantly improve the integration accuracy for small to medium systems, making it feasible to use up to 10 fs time steps even for compact proteins like GB1. However, such a trick is not as effective for large, compact systems like DHFR. In Figure 9, we also examine the long-term energy drift of CMD and TAMD simulations. The drift is computed by linear regression fitting of the energies during 100 ps NVE simulations. It shows that long-term energy drift in TAMD is also very small and lower than in CMD with a time step of 2 fs or larger. Even with a time step of 5 fs, the energy drift of TAMD is below 0.1 kcal/ps for compact proteins like GB1.

The different dependence of CMD and TAMD (with and without increasing the mass of hydrogen atoms) on time step size is a consequence and reflection of the complicated hierarchy of fast motions present in biological systems.61 Constraining the hydrogen-attached bonds removes the fastest motions in this hierarchy and allows time steps up to 3 fs. Constraints on all bond and angle degrees of freedom further extends allowable time steps up to about 5 fs, beyond which rotations of groups with small inertia (e.g., hydroxyl and methyl groups) and collisions between nonbonded heavy atoms become the dominating factors that limit the integration time step. In CMD simulations with SHAKE applied to hydrogen-attached bonds, the accuracy limiting factor is mainly bond and angle degrees of freedom, which is independent of the system size and topology. As such, the energy fluctuation of CMD simulations only becomes dependent on the size and topology of the system when the time steps approach 5 fs, where factors such as nonbonded heavy atom collisions become significant. On the contrary, the accuracy of TAMD simulations depends strongly on the system size and topology for all time steps, as all bond and angle degrees of freedom are constrained and the remaining accuracy limiting factors are all system specific. It is interesting to note that DHFR does not suffer much from the increased size compared to GB1 when the hydrogen mass is not increased. The reason is probably because rotations of light tip groups becomes a dominating factor when the systems are similarly compact. One can simply increase the mass of all hydrogen atoms to suppress the impact of rotations of the tip groups and increase the dynamics accuracy, demonstrated in Figure 8b. In this case, the accuracy is mainly limited by nonbonded heavy atom collisions. This is the reason why smaller time steps are required for large systems like DHFR compared to medium-size systems like GB1 even though the degree of compactness is comparable.

#### *Folding of Simple Peptides*

A key objective in the current work is to use TAMD to extend our sampling capability and thus study biologically important processes that involve large conformational changes such as the problems of protein folding and *ab initio* structure prediction. Arguably, TAMD is most useful when combined with efficient implicit solvent models that have been rapidly advanced in the last few years.37 Although reaching the above-stated goal probably requires further optimization of the implicit solvent models, particularly in rebalancing the hydrogen bonding interactions and solvent–solute dispersion interactions,62 the efficacy of TAMD in augmenting conformational sampling, as well as possible limiting factors can be examined by folding simulations of small peptides and refinement experiments of protein NMR structures (see next section). Several peptides were used, including poly-alanine, (Ala)10, poly-valine, (Val)10, and the WALP transmembrane peptides.59,60 The results show that larger time steps can be effectively used in TAMD and increase the conformation space sampled for given computational time. In particular, for simple peptides like poly-alanine, as shown in Figures 10 and 11, time steps up to 20 fs can be reliably used, leading to enhancement of at least 10-fold compared to typical CMD. Note that the total number of degrees of freedom in internal coordinates is typically only about one-tenth

![](_page_8_Figure_9.jpeg)

**Figure 9.** Energy drifts during 100 ps NVE simulations of GB1 in Cartesian and internal coordinates. The same setup was used as in Figure 8. Note that CMD simulation with 5-fs time step diverged shortly after 3 ps and the drift shown was computed using only the first 3 ps of the trajectory.

![](_page_9_Figure_2.jpeg)

**Figure 10.** RMSD values of the lowest temperature ensembles as a function of the simulation time during several REX folding simulations of the (Ala) $_{10}$  peptide. Eight replicas were used over a temperature range of 300 to 800 K. All simulations were carried out in vacuum with no cutoff for nonbonded interactions. Time steps of 2 to 20 fs were used in the TAMD simulations. Note that the total computational cost is inversely proportional to the time step size. Control repeat REX simulations (TAMD with time steps of 5 and 20 fs) show that the average "folding" time of (Ala) $_{10}$  in vacuum with the above REX setup is  $29 \pm 5$  ps.

of that in Cartesian space. As such, the computational cost of REX simulations can be further reduced as smaller numbers of replicas are required for the same temperature range (optimal number of replicas is proportional to  $\sqrt{N_{\rm dof}}^{63}$ ). Furthermore, there is also some potential increase in conformational sampling efficiency, because only the most relevant degrees of freedoms are explicitly sampled in TAMD.

Increase in the integration time steps is more limited for more complex peptides, especially those containing residues with bulky side chains. For example, as shown in Figure 11 for (Val)<sub>10</sub>, even though a more rapid initial decrease in RMSD values is consis-

![](_page_9_Figure_6.jpeg)

**Figure 11.** RMSD values of the lowest temperature ensembles as a function of the simulation time during several REX folding simulations of the (Val)<sub>10</sub> peptide. The same simulation setup was used as in Figure 10. Control repeat REX simulations (CMD and TAMD with a 2-fs time step) show that the average "folding" time of (Val)<sub>10</sub> in vacuum with current REX setup is about 250–500 ps if backbone RMSD value of 1.0 Å is used as a criteria for folded helical states.

![](_page_9_Figure_8.jpeg)

**Figure 12.** REX folding simulations of WALP16 in GBSW implicit membrane. Sixteen replicas were used. The temperature range was 300 to 600 K for the CMD simulations and 300 to 1000 K for the TAMD simulations. Nonbonded interactions were smoothly switched off at 20 Å. The thickness of the implicit membrane was set to be 28.0 Å. The initial extended peptide chain was placed in perpendicular to the membrane plane.

tently observed in the TAMD simulations, the peptides appear to spend more time in "nearly folded" intermediate states (with backbone RMSD values around 2.0 Å) before reaching the fully folded helical states. It is evident that such traps in nearly folded compact states become more severe when larger time steps are used. More severe limitations were observed in the REX folding simulations of the WALP16 trans-membrane peptide (sequence: Ace-GWW(LA)<sub>5</sub>WWA-NMe) using a GSBW implicit membrane model,<sup>58</sup> as shown in Figure 12. It appears that TAMD REX simulations with a 5-fs time step have problems in sampling the fully folded helical states. Previous results (see above) indicate that the accuracy of the dynamics is likely to be sufficient for small peptides even with time steps up to 10 fs. To verify that the accuracy of the dynamics is not the reason that prevents sampling of the native states, energy conservation properties were reexamined for WALP16 in the implicit membrane over a whole range of representative conformations during the folding process. The impact of nonbonded interaction truncation was also examined. The results are shown in Figure 13. Clearly, the accuracy of the dynamics at 300 K does not deteriorate appreciably when the time step is increased from 1 to 5 fs. Even with time steps as large as 10 fs, a decrease in the accuracy is only significant for the most compact state (Fig. 13If). At 800 K, the accuracy of the dynamics decreases more rapidly but still appears to be reasonable with time steps up to 5 fs. Therefore, these folding simulations seem to suggest that to maintain efficient sampling of the conformational

![](_page_10_Figure_2.jpeg)

**Figure 13.** Energy conservation properties as a function of integration time step for different conformational states of the WALP16 peptide. Membrane GBSW was used and the hydrogen mass was increased to 6.0 in all calculations. The conformational states were first extracted from the REX folding simulations and then equilibrated at the appropriate temperatures. The energy fluctuation was computed in the same way as in Figure 8.

space the time steps might be more limited than previous estimations based purely on energy conservation properties.<sup>2,16</sup>

#### Refinement of NMR Structures of a Protein Domain

As a realistic test of the efficacy of TAMD in enhancing conformational sampling, refinement of NMR structures of a heat-shock protein (HSP) redox-switch domain was carried out. The protein domain has a novel fold, with two helices at right angles to each other, a two-stranded  $\beta$ -hairpin and a third helix at the C terminus.<sup>64</sup> In this test, we used an NMR restraint set obtained during very early stages of the structure determination.<sup>65</sup> Conventional structure calculations using this initial sparse NOE set were unable to identify a unique topology due to a lack of unambiguous long range NOEs (see Table 1). An accurate overall topology was

eventually obtained through laborious hand identification of more long-range NOEs. However, it was demonstrated previously that high-quality native-like models could be obtained through REX refinement with GB implicit solvent (REX/GB).<sup>65</sup> In this test, the same ensemble of initial structures were used (run I of Table 2 in ref. 66). Even though the optimal REX setup (e.g., number of replicas and the temperature range) might be different with TAMD, the same REX setup was used in this comparison, that is, 16 replicas were used from 300 to 600 K. The MD length of each REX step was 1.0 ps. A total of 1000 REX steps were carried out with the last 200 steps as production cycles. As shown in Table 1, the refinement results demonstrate that larger time steps can be effectively used to achieve similar quality of structure convergence. Note that average structures from the TAMD ensembles

**Table 1.** Comparison of Structural and NOE Statistics of Various Structure Ensembles in the REX/GB Refinement Tests.

| Ensembles | RMSD <sup>a</sup> | NOE <sup>b</sup> |
|-----------|-------------------|------------------|
| Initial   | $8.8 \pm 5.8$     | 2.1/0.021        |
| CMD/2 fs  | $2.2 \pm 2.7$     | 4.4/0.020        |
| TAMD/2 fs | $3.0 \pm 2.4$     | 5.0/0.024        |
| TAMD/3 fs | $3.1 \pm 0.9$     | 4.5/0.021        |
| TAMD/4 fs | $3.3 \pm 2.7$     | 5.7/0.026        |
| TAMD/5 fs | $3.1 \pm 4.3$     | 5.1/0.025        |

<sup>&</sup>lt;sup>a</sup>Backbone RMSD of the ensemble average from the final NMR structure  $^{64}$  ± backbone RMS fluctuation around the average (in Å). Only structured regions (residues 7–60) were included in the RMSD calculations.

have slightly larger backbone RMSD values from the final NMR structure (e.g., greater than 3.0 Å) compared to previous ensembles obtained using CMD (between 2.0 to 3.0 Å). However, the difference is mainly in the compactness of the structures rather than the tertiary fold. Also note that smaller structural variations in the TAMD ensemble with a 3-fs time step are due to the fact that the lowest temperature ensemble was occupied mainly by a single replica during the production period.<sup>66</sup> In addition, when time steps of 5 fs or larger were used, some limitations were observed both in sampling the native states and in ranking different conformers (i.e., most native-like models do not have as high occupancy in the lowest temperature windows). Such limitations are consistent with what was observed in the folding simulations of small peptides. Similarly, energy conservation does not seem to be the limiting factor. Instead, the overall flexibility of the molecules with the rigid body approximation might still not be as good as in Cartesian space. Highly collective motions in torsion space and long-range collisions have larger impacts on the overall flexibility with larger time steps, which could limit the ability of making small local conformational adjustment in compact states. More studies are necessary to fully understand the limitations of TAMD, especially in the optimal ways of using TAMD to augment the conformational sampling.

### Conclusion

An efficient recursive NEIMO algorithm<sup>13</sup> for molecular dynamics in torsion space has been implemented in a CHARMM TAMD module to address one of the major limitations in traditional MD simulation, that accessible simulation time scale of MD in Cartesian coordinates is severely limited by the femtosecond time step required by bond and angle degrees of freedom. Carrying out dynamics simulation in torsion space eliminates these high-frequency degrees of freedom and allows larger time steps. In addition, only the most relevant degrees of freedom are explicitly sampled in TAMD, which can also potentially expedite the conformational search. However, with Cartesian-based force fields,

the rigid covalent geometry approximation introduces significant distortions on the energy surface, rendering direct application of TAMD untenable. We have extended a "projection" approach<sup>38</sup> to accurately reconstruct local potential surfaces and restore local covalent flexibility implicitly. The new approach first employs softened van der Waals and electrostatic interactions to reduce the dramatic barrier increases on the potential surface. Torsion crossterm corrections, constructed from local molecular fragments, are then applied using the CMAP facility in CHARMM.50-52 Backbone  $\phi/\psi$  and side chain  $\chi_1$  torsion profiles and dynamics properties of dipeptide models demonstrate that the resulting ICFF is an accurate projection of the original Cartesian force field. Application of TAMD to folding simulations of small peptides and refinement of NMR structures have been very encouraging, showing that it is possible to use much larger time steps effectively, up to 20 fs for simple peptides and more limited for compact and large systems. Additional benefits come from about a 10-fold reduction of the number of degrees of freedom in internal coordinates and exclusive sampling of conformationally relevant degrees of free-

Current studies have also suggested that more investigation is necessary for better understanding of the intrinsic dynamic properties in torsion space. In particular, the impacts and limitations of larger time steps need to be further studied not just in the context of integration accuracy but also in terms of conformational sampling efficiency. Current results indicate that the time steps might be more limited than previous estimations based purely on the energy conservation. 2,16 The high collectivity of motions in torsion space and collisions between long-range atoms can reduce the overall flexibility of the system, especially when larger time steps are employed. In addition, the choice of fixed covalent geometry is an important issue that needs further investigation. The current choice of backbone geometry that allows for stable  $\alpha$ -helices appears to be a reasonable compromise, but it is not clear that such a choice is optimal. Last, no systematic studies have been carried out to understand intrinsic properties and limitations of dynamics in reduced variable models. Such knowledge is essential for the effective use of TAMD. Important thermodynamic properties such as entropy and free energy need to be compared with and without the rigid body approximation to identify possible thermodynamic consequences of freezing bond and angle degrees

One of the most powerful ways of applying TAMD is to combine it with efficient implicit solvent models and advanced sampling techniques like the REX method. A particular useful application of such combination is to predict tertiary structures of proteins from their primary sequences. For this, one will need not only sampling capability but also accurate force fields. Recent years have seen many successful examples of accurate *ab initio* prediction of structures for several miniproteins, such as protein A,  $^{67,68}$  Trp-cage,  $^{69}$   $\alpha\beta\beta$ -motifs  $^{70,71}$  and the fd Coat protein. The existing GB implicit solvent models need to be further optimized, particularly in rebalancing the hydrogen bonding interactions and solvent—solute dispersion interactions. Folding simulations of these miniproteins can be used as test cases for identifying the optimal sampling protocol and optimizing the force field. TAMD

<sup>&</sup>lt;sup>b</sup>Average number of NOE restraints violated by more than 0.2 Å/RMSD of NOE restraints for all structures in the ensemble (in Å). Note that there is no NOE violated by over 0.5 Å in any ensemble.

will also be able to augment enhanced sampling in refining NMR structures with more accurate force fields. It has been increasingly recognized that refinement with GB implicit solvent can lead to significant improvement in protein NMR structures,73 especially when the experimental data is limited such as in early stages of NMR structure determination and for large proteins.65,66 Using TAMD instead of CMD in these refinement calculations is straightforward and can lead to several fold reduction of computational cost.

# **Acknowledgments**

The authors acknowledge helpful discussions with Abhinandan Jain.

## **References**

- 1. Pear, M.; Weiner, J. J Chem Phys 1979, 71, 212.
- 2. Mazur, A. K.; Abagyan, R. A. J Biomol Struct Dyn 1989, 6, 815.
- 3. Gibson, K. D.; Scheraga, H. A. J Comput Chem 1990, 11, 468.
- 4. Ryckaert, J. P.; Ciccotti, G.; Berendsen, H. J. C. J Comput Phys 1977, 23, 327.
- 5. Barth, E.; Kuczera, K.; Leimkuhler, B.; Skeel, R. D. J Comput Chem 1995, 16, 1192.
- 6. Tuckerman, M.; Berne, B. J.; Martyna, G. J. J Chem Phys 1992, 97, 1990.
- 7. Watanabe, M.; Karplus, M. J Chem Phys 1993, 99, 8063.
- 8. Izaguirre, J. A.; Catarello, D. P.; Wozniak, J. M.; Skeel, R. D. J Chem Phys 2001, 114, 2090.
- 9. Mazur, A. K.; Dorofeev, V. E.; Abagyan, R. A. J Comput Phys 1991, 92, 261.
- 10. Bae, D. S.; Haug, E. J. Mech Struct Mach 1987, 15, 359.
- 11. Bae, D. S.; Haug, E. J. Mech Struct Mach 1988, 15, 481.
- 12. Jain, A. J Guidance Control Dyn 1991, 14, 531.
- 13. Jain, A.; Vaidehi, N.; Rodriguez, G. J Comput Phys 1993, 106, 258.
- 14. Rodriguez, G.; Jain, A.; Kreutz–Delgado, K. Int J Rob Res 1991, 81, 371.
- 15. Mazur, A. K. In Computational Biochemistry and Biophysics; Beker, O. M., Ed.; Marcel Dekker: New York, 2001.
- 16. Mathiowetz, A. M.; Jain, A.; Karasawa, N.; Goddard, W. A. Proteins 1994, 20, 227.
- 17. Mazure, A. K. J Am Chem Soc 1998, 120, 10928.
- 18. Rice, L. M.; Bru¨nger, A. T. Proteins 1994, 19, 277.
- 19. Stein, E. G.; Rice, L. M.; Bru¨nger, A. T. J Magn Reson 1997, 124, 154.
- 20. Gu¨ntert, P.; Mumenthaler, C.; Wu¨thrich, K. J Mol Biol 1997, 273, 283.
- 21. Schwieters, C. D.; Clore, G. M. J Magn Reson 2001, 152, 288.
- 22. Klepeis, J. L.; Floudas, C. A. Comput Chem Eng 2000, 24, 1761.
- 23. Nemethy, G.; Gibson, K. D.; Palmer, K. A.; Yoon, C. N.; Paterlini, G.; Zagari, A.; Rumsey, S.; Scheraga, H. A. J Phys Chem 1992, 96, 6472.
- 24. Vaidehi, N.; Jain, A.; Goddard, W. A., III. J Phys Chem 1996, 100, 10508.
- 25. Mazur, A. K. J Comput Chem 1997, 18, 1354.
- 26. Bertsch, R. A.; Vaidehi, N.; Chan, S. I.; Goddard, W. A. Proteins 1998, 33, 343.
- 27. Vaidehi, N.; Floriano, W. B.; Trabanino, R.; Hall, S. E.; Freddolino, P.; Choi, E. J.; Zamanakos, G.; Goddard, W. A. Proc Natl Acad Sci USA 2002, 99, 12622.
- 28. Ponder, J. W.; Case, D. A. Adv Protein Chem 2003, 66, 27.

- 29. MacKerell, A. D., Jr. J Comput Chem 2004, 25, 1584.
- 30. Brooks, B. R.; Bruccoleri, R. E.; Olafson, B. D.; States, D. J.; Swaminathan, S.; Karplus, M. J Comput Chem 1983, 4, 187.
- 31. MacKerell, A. D., Jr., Bashford, D., Bellott, M., Dunbrack, R. L., Evanseck, J. D., Field, M. J., Fischer, S., Gao, J., Guo, H., Ha, S., Joseph-McCarthy, D., Kuchnir, L., Kuczera, K., Lau, F. T. K., Mattos, C., Michnick, S., Ngo, T., Nguyen, D. T., Prodhom, B., Reiher, W. E., III, Roux, B., Schlenkrich, M., Smith, J. C., Stote, R., Straub, J., Watanabe, M., Wiorkiewicz-Kuczera, J., Yin, D., Karplus, M. J Phys Chem B 1998, 102, 3586.
- 32. Cornell, W.; Cieplak, P.; Bayly, C.; Gould, I.; Merz, K., Jr.; Ferguson, D.; Spellmeyer, D.; Fox, T.; Caldwell, J.; Kollman, P. J Am Chem Soc 1995, 117, 5179.
- 33. Jorgensen, W. L.; Tirado–Rives, J. J Am Chem Soc 1988, 110, 1657.
- 34. van Gunsteren, W. F. GROMOS. Groningen Molecular Simulation Program Package; University of Groningen: Groningen, 1987.
- 35. Still, W. C.; Tempczyk, A.; Hawley, R. C.; Hendrickson, T. J Am Chem Soc 1990, 112, 6127.
- 36. Bashford, D.; Case, D. A. Annu Rev Phys Chem 2000, 51, 129.
- 37. Feig, M.; Brooks, C. L., III. Curr Opin Struct Biol 2004, 14, 217.
- 38. Katritch, V.; Totrov, M.; Abagyan, R. J Am Chem Soc 2003, 24, 254.
- 39. Dunfield, L. G.; Burgess, A. W.; Scheraga, H. A. J Phys Chem 1978, 82, 2609.
- 40. Hinsen, K.; Kneller, G. R. Phys Rev E 1995, 52, 6868.
- 41. Rodriguez, G.; Jain, A.; Kreutz–Delgado, K. J Astronaut Sci 1992, 40, 27.
- 42. Meirovitch, L. Methods of Analytical Dynamics; McGraw–Hill: New York, 1970.
- 43. Allen, M. P.; Tildesley, D. J. Computer Simulation of Liquids; Clarendon Press: Oxford, 1987.
- 44. Fincham, D. Mol Simulat 1992, 8, 165.
- 45. Berendsen, H. J. C.; Postma, J. P. M.; van Gunsteren, W. F.; Dinola, A.; Haak, J. R. J Chem Phys 1984, 81, 3684.
- 46. Nose´, S. J Chem Phys 1984, 81, 511.
- 47. Hoover, W. G. Phys Rev A 1985, 31, 1695.
- 48. Sugita, Y.; Okamoto, Y. Chem Phys Lett 1999, 314, 141.
- 49. Press, W. H.; Teukolsky, S. A.; Vetterling, W. T.; Flannery, B. P. Numerical Recipies in C: The Art of Scientific Computing; Cambridge University Press: Cambridge, UK, 1992, 2nd ed.
- 50. Feig, M.; MacKerell, A. D., Jr.; Brooks, C. L., III. J Phys Chem 2003, 107, 2831.
- 51. MacKerell, A. D., Jr.; Feig, M.; Brooks, C. L., III. J Am Chem Soc 2004, 126, 698.
- 52. MacKerell, A. D., Jr.; Feig, M.; Brooks, C. L., III. J Comp Chem 2004, 25, 1400.
- 53. Gronenborn, A. M.; Filpula, D. R.; Essig, N. Z.; Achari, A.; Whitlow, M.; Wingfield, P. T.; Clore, G. M. Science 1991, 253, 657.
- 54. Sawaya, M. R.; Kraut, J. Biochemistry 1997, 36, 586.
- 55. Feig, M.; Karanicolas, J.; Brooks, C. L., III. 2001 MMTSB Tool Set, MMTSB NIH Research Resource, The Scripps Research Institute, 2001.
- 56. Feig, M.; Karanicolas, J.; Brooks, C. L., III. J Comp Graph Modl 2004, 22, 337.
- 57. Im, W.; Lee, M. S.; Brooks, C. L., III. J Comput Chem 2003, 24, 1691.
- 58. Im, W.; Feig, M.; Brooks, C. L., III. Biophys J 2003, 85, 2900.
- 59. de Planque, M. R. R.; Killian, J. A. Mol Membr Biol 2003, 20, 271.
- 60. Im, W.; Brooks, C. L., III. Proc Natl Acad Sci USA 2005, 102, 6771.
- 61. Mazur, A. K. J Phys Chem 1998, 102, 473.
- 62. Im, W.; Chen, J.; Brooks, C. L., III. In Peptide H-Bonds and Peptide Solvation, Baldwin, R. E., Baker, D. J., Ed.; Elsevier, 2005.

- 63. Hukushima, K.; Nemoto, K. J Phys Soc Jpn 1996, 65, 1604.
- 64. Won, H.-S.; Low, L. Y.; De Guzman, R.; Martinez–Yamout, M.; Jakob, U.; Dyson, H. J. J Mol Biol 2004, 341, 893.
- 65. Chen, J.; Won, H.-S.; Im, W.; Dyson, H. J.; Brooks, C. L., III. J Biomol NMR 2004, 31, 243.
- 66. Chen, J.; Im, W.; Brooks, C. L., III. J Am Chem Soc 2004, 126, 16038.
- 67. Jang, S. M.; Kim, E.; Shin, S.; Pak, Y. J Am Chem Soc 2003, 125, 14841.
- 68. Vila, J. A.; Ripoll, D. R.; Scheraga, H. A. Proc Natl Acad Sci USA 2003, 100, 14812.
- 69. Pitera, J. W.; Swope, W. Proc Natl Acad Sci USA 2003, 100, 7587.
- 70. Abagyan, R. A.; Totrov, M. J Chem Phys 1999, 151, 402.
- 71. Jang, S.; Shin, S.; Pak, Y. J Am Chem Soc 2002, 124, 4976.
- 72. Im, W.; Brooks, C. L., III. J Mol Biol 2004, 337, 513.
- 73. Xia, B.; Tsui, V.; Case, D. A.; Dyson, H. J.; Wright, P. E. J Biomol NMR 2002, 22, 317.