# **Software News and Updates**

# ICFF: A New Method to Incorporate Implicit Flexibility into an Internal Coordinate Force Field

## VSEVOLOD KATRITCH,1,\* MAXIM TOTROV,2 RUBEN ABAGYAN1,2

<sup>1</sup>Department of Molecular Biology, The Scripps Research Institute, 10550 North Torrey Pines, TPC-28, La Jolla, California 92037 <sup>2</sup>Molsoft, 3366 Torrey Pines Court, La Jolla, California 92037

Received 10 October 2001; Accepted 12 February 2002

**Abstract:** We introduce a new method to accurately "project" a Cartesian force field onto an internal coordinate molecular model with fixed-bond geometry. The algorithm automatically generates the Internal Coordinate Force Field (ICFF), which is a close approximation of the "source" Cartesian force field. The ICFF method reduces the number of free variables in a model by at least 10-fold and facilitates the fast convergence of geometry optimizations, an advantage that is critical for many applications such as the docking of flexible ligands or conformational modeling of macromolecules. Although covalent geometry is fixed in an ICFF model, implicit flexibility is incorporated into the force field parameters in the following two ways. First, we formulate an empirical torsion energy term in ICFF as a sixfold Fourier series and develop a procedure to calculate the Fourier coefficients from the conformational energy profiles of the fully flexible Cartesian model. The ICFF torsion parameters thus represent not only torsion component of the source force field, but also bond bending, bond stretching, and "1-4" van der Waals interactions. Second, we use a soft polynomial repulsion function for "1-5" and "1-6" interactions to mimic the flexibility of bonds, connecting these atoms. Also, we suggest a way to use a local part of the Cartesian force field to automatically generate fixed covalent geometries, compatible with the ICFF energy function. Here, we present an implementation of the ICFF algorithm, which employs the MMFF94s Cartesian force field as a "source." Extensive benchmarking of ICFF with a representative set of organic molecules demonstrates that the implicit flexibility model accurately reproduces MMFF94s equilibrium conformational energy differences (RMSD  $\sim 0.64$  kcal) and, most importantly, detailed torsion energy profiles (RMSD  $\sim 0.37$  kcal). This accuracy is characteristic of the method, because all the ICFF parameters (except one scaling factor in the "1-5,1-6" repulsion term) are derived directly from the source Cartesian force field and do not depend on any particular molecular set. In contrast, the rigid geometry model with the MMFF94s energy function yields highly biased estimations in this test with the RMSD exceeding 1.2 kcal for the equilibrium energy comparisons and ~3.4 kcal for the torsion energy profiles.

© 2002 Wiley Periodicals, Inc. J Comput Chem 24: 254-265, 2003

Introduction Cartesian Force Fields

Force field-based molecular mechanics (MM) and molecular dynamics (MD) modeling provide a number of practical tools for structural biology and drug design. These applications include the prediction of protein structures<sup>1–4</sup> and conformational changes,<sup>5,6</sup> studies of protein–protein and protein–nucleic acid interactions,<sup>7–9</sup> ligand docking,<sup>10,11</sup> and virtual screening,<sup>7,12</sup> and many more. The predictive power of these computational tools largely relies on the accuracy of the employed force field, i.e., the ability to discriminate correctly between low-energy conformations and high-energy conformations. Equally important for biologically relevant predictions is the ability of a molecular model to reach convergence of the global energy optimizations in MM or conformational equilibrium in MD, which requires ample sampling of the conformational space on a reasonable time scale.<sup>13</sup>

Recent years have seen dramatic improvements in the methodology of molecular force field development and the scope of molecular structures covered by parameterization. Hampy of the recent force fields have been parameterized using extensive quantum mechanics (QM) calculations as well as a variety of experimental results, such as X-ray geometries, molecular spectra, liquid state thermodynamic data, etc. Among the widely parameterized force field functions one can mention several modifications of CHARMM, 15–17 MM3 and MM4, 18–20 and Amber, 21–23 as well as

Correspondence to: R. Abagyan; e-mail: abagyan@scripps.edu

<sup>\*</sup>Present address: Plexus Vaccine Inc., 11770 Bernardo Plaza Ct., San Diego, CA 92128

the more recently developed CFF95,<sup>24</sup> OPLS-AA,<sup>25–27</sup> and MMFF94 force fields.<sup>28–35</sup> These, and other potential functions, are employed in various molecular mechanics and molecular dynamics applications yielding accuracy, comparable in some cases to QM calculations. For example, a general-purpose MMFF94 force field reproduces conformational energy differences with an RMSD of about 0.3–0.5 kcal for a set of 146 diverse organic compounds included in the parameterization of this force field.<sup>32,33</sup> For the molecules not included in their training sets, the accuracy of the best available force fields can still be estimated at about 1–1.5-kcal level.

The accuracy and transferability of the force fields are constantly improving through the development of novel algorithms for large scale self-consistent parameterization of potential functions, <sup>24,28,36</sup> as well as by extending individual training sets to new classes of molecules. <sup>20,37–39</sup>

## Modeling in Internal Coordinates

One of the main challenges in the molecular modeling applications is the giant size of conformational space and the complexity of the energy landscapes, which make extensive sampling and convergent global energy optimizations for most biologically interesting molecular systems prohibitive. An obvious way to reduce the conformational space is to constrain high-frequency modes of a molecular model. 6,40-42 Thus, switching from a Cartesian to an internal coordinate representation with fixed bond lengths and bond angles one can cut the number of free variables by about sixfold. Also, it is often natural to enforce either planar or tetrahedral covalent geometry for most atoms, as well as to impose planarity on some double-bonded structures and aromatic rings. Commonly used in protein modeling, 43-46 such a rigid covalent geometry approximation accounts for about a 10-fold decrease in the number of free variables. Further dramatic reductions of conformational space in internal coordinates can be achieved naturally by "activating" only a specific subset of torsion angles of the molecular system. Thus, one can unfix only torsion angles of protein side chains and loop backbones for homology modeling<sup>47,48</sup> or binding site residues and ligands for docking.<sup>49,50</sup> The advantage of a torsion angle representation is evident not only in the smaller dimensionality of the sampling space and faster energy evaluations at each step, but also in more reliable local minimizations, which have much larger radii of convergence than Cartesian space local minimizations.<sup>44</sup>

Combined with effective energy optimization procedures, for example, Monte Carlo Minimization<sup>51</sup> (MCM) or Optimal Bias MCM,<sup>52,53</sup> the internal coordinate approach yields reproducible global convergence in a reasonable time scale for such computationally demanding problems as *ab initio* folding of a 23-residue beta-beta-alpha peptide<sup>53</sup> (a few hours) or the flexible docking of organic compounds<sup>12,49</sup> (a few minutes).

# Force Field Implementation in Internal Coordinates

The development of accurate and transferable force fields compatible with the rigid geometry approximation is currently the major bottleneck for the torsion coordinate model.

Historically, one of the first widely available force fields for peptides, ECEPP, was derived in torsion angle space. 54,55 The ECEPP parameters were adjusted empirically to reproduce X-ray geometries and experimental conformational energies of peptides in terms of the rigid geometry model. A more recent version of this potential function, ECEPP/3,56 is effectively used in various aspects of protein modeling, including peptide folding, 43,53 modeling by homology, 4,47 protein-protein, 57 and protein-peptide docking. 58-60 The success of ECEPP-based applications demonstrates the adequateness of rigid geometry approximations, at least for peptide modeling. Although accurate in protein modeling, the parameterization of the ECEPP force field is very hard to extend to other biopolymers and organic compounds. In the absence of an accurate general-purpose torsion force field one has to use the torsion component of the Cartesian force fields for modeling of nonpeptide molecules in torsion coordinate space. 50,60

Unfortunately, several well-known issues make the application of Cartesian force fields to the torsion angle modeling questionable. 61 First, van der Waals repulsion between atoms, separated by three covalent bonds ("1-4" contacts), which is partially relaxed in the Cartesian model due to bond and valence angle flexibility, would result in inappropriate high-energy clashes if the covalent geometry is fixed. One would have to adjust the torsion energy function to compensate for these discrepancies.<sup>54</sup> Second, large energy anomalies in the rigid model can originate from steric clashes between atoms separated by four or five covalent bonds ("1-5" and "1-6" clashes). Unlike "1-4" contacts, these interactions depend on more than one torsion angle, and the corresponding energy deviations cannot be compensated by the torsion energy function.<sup>62</sup> The third problem is how to define an adequate covalent geometry for the molecule. Conformational energy in a rigid model is hypersensitive to the bonded geometry because small variations in bond length and bond angles may cause severe van der Waals clashes. To reduce this effect, the covalent geometry for each molecule has to be derived individually to comply with the force field parameters.62

Here we address these issues by introducing a new concept of the accurate *projection* of a "source" Cartesian force field into the torsion angle conformational energy function, ICFF.

Our projection algorithm implicitly treats "1-4" van der Waals interactions as a part the ICFF composite torsion energy term, together with the bond bending, bond stretching, and out-of-plane energy terms. We also consider "1-5,1-6" van der Waals contacts as a part of bonded energy potential and treat "1-5,1-6" repulsion with a soft empirical potential to mimic bond flexibility in clashed conformations. The energy function with softer "1-5,1-6" van der Waals repulsion and no explicit "1-4" term is much less sensitive to the fixed bond geometry definition. Finally, we show that a fast Cartesian optimization using only local terms of the "source" force field provides consistent covalent geometry for ICFF molecular models. Parameterization of the internal coordinate energy function in our algorithm is derived completely from a source Cartesian force field, with only one adjustable scaling factor. This suggests high transferability of the method to all molecular structures covered by the Cartesian force field. There is always a choice between a static file with precomputed torsion parameters for all combinations of four atom types and a run-time calculation of the parameters. In this article we argue for a run-time generation, because the

![](_page_2_Figure_2.jpeg)

**Figure 1.** Parameterization of the ICFF torsion energy function for the torsion angle  $\Omega$  (C—N—C—H) in the alanine dipeptide analog. Labels are shown for seven atoms of the "torsion fragment," considered in Cartesian energy optimizations. MMFF94 atom types are given in brackets.

ICFF parameterization algorithm is fast enough to be performed "on the fly" for each torsion type in a molecular object and, most importantly, overcomes the limitation of exceedingly large combinatorics of possible torsion types, eliminating the need to store the torsion force field parameters.

The current implementation of the ICFF method is based on the MMFF94s Cartesian force field.<sup>28</sup> The accuracy of the resulting energy function is compared to the accuracy of the source force field on a large conformational set of diverse organic compounds.<sup>33</sup>

# **ICFF** Formulation

# Torsion Energy Term

The definition of the composite torsion energy term, which accounts implicitly for bond bending, bond stretching and "1-4" van der Waals interactions, is a central part of our force field projection algorithm.

We define the ICFF empirical torsion term as a sixfold Fourier series:

$$E_{\text{tor}} = C_0 + \sum_{k=0}^{k=6} (A_k \cos k\theta + B_k \sin k\theta), \tag{1}$$

where  $\theta$  is a torsion angle for rotation about a specific bond. Fourier coefficients for the bond are calculated in the following manner:

1. The two atoms flanking the bond and all their immediate neighbors (4 to 8 atoms overall) are copied into a separate object called a *torsion fragment* with original atom types and bond topology preserved (Fig. 1).

- A grid search from 0° to 360°, using a 12° increment, of the torsion angles specifying the rotation of the bond is performed.
  This results in 30 conformations of the torsion fragment.
- 3. For each value of the torsion angle *constrained*, Cartesian *minimization* of the fragment object is performed with the *local* part of the "source" force field potential. The local energy function here includes all bonded terms and "1-4" van der Waals interactions, but leaves out other van der Waals interactions and all electrostatics.
- The best energy values found by the Cartesian minimization for each increment of the torsion angle, constitute an *energy profile* for this torsion fragment.
- 5. This torsion energy profile is subsequently approximated by the Fourier series function (1), yielding torsion parameters {A<sub>k</sub>, B<sub>k</sub>}. Though a direct Fourier transform usually produces satisfactory coefficients for eq. (1), we found a general least-squares approach, where the merit function is weighted by the inverse relative energy value, more appropriate for this approximation to ensure better accuracy in the low-energy part of the torsion profile (see Results and Discussion).

Note, that the torsion parameters produced by this algorithm do not depend on the initial covalent geometry of the model, because the "torsion fragment" geometries are reoptimized at step 3 for each value of the torsion angle. Therefore, the initial model can be built with any generic algorithm, provided that topology of atom connections is correct. We use the above protocol to derive torsion parameters for each ICFF *bond type* in the molecule, a bond type being defined by atom types and connection topology of the atoms included in the torsion fragment.

## Empirical Term for "1-5" and "1-6" Interactions

The van der Waals repulsion between atoms connected by four or five bonds ("1-5,1-6" interactions) is treated in ICFF by a special empirical term:

$$E_{\text{vw}(1-5,1-6)} = \varepsilon_{IJ} \left[ C \left( \left( \frac{R_{IJ} - R_{IJ}^*}{R_{IJ}^0 - R_{IJ}^*} \right)^2 - 1 \right) + (1 - C) \left( \left( \frac{R_{IJ} - R_{IJ}^*}{R_{IJ}^0 - R_{IJ}^*} \right)^3 - 1 \right) \right], \text{ if } R_{IJ} < R_{IJ}^*, \quad (2)$$

where  $\varepsilon_{IJ}$  is the minimal energy in the van der Waals surface for the corresponding atom types I and J;  $R_{IJ}^*$  is the atom-atom distance at which this minimum occurs;  $R_{IJ}^0$  is the distance at which vdW surface crosses zero. The potential energy surface is defined here by the van der Waals function of the source Cartesian force field. A scaling factor C (0  $\leq C \leq$  1) controls the balance between harmonic and cubic contributions, and is the only adjustable parameter in the "1-5,1-6" repulsion term.

Note that the energy function (2) is formulated to satisfy the following criteria: (1)  $E_{\text{vw}(1-5,1-6)}$  is equal to  $\varepsilon_{IJ}$  and has zero derivative when  $R_{IJ} = R_{IJ}^*$ ; (2)  $E_{\text{vw}(1-5,1-6)}$  crosses zero at the same point,  $R_{IJ} = R_{IJ}^0$ , as the original van der Waals function.

Thus defined, the empirical repulsion term forms a continuous curve with the attractive part of the van der Waals potential,  $R_{IJ} > R_{IJ}^*$ , which is kept unchanged. The repulsion potential behaves as

a soft cubic polynomial to mimic the flexibility of bonds connecting 1-5 or 1-6 atom pairs. The van der Waals repulsion between hydrogen bond atoms are treated as a special case by the original "hard" van der Waals term to compensate strong attraction between donor protons and acceptor atoms. Note, that the soft  $E_{\text{vw}(1-5,1-6)}$  term has a finite value at  $R_{IJ}=0$ ; therefore, it requires an appropriately capped electrostatic term to prevent highly charged atoms form occasionally collapsing into each other in MC or MD simulations.

#### Covalent Geometry

Fixed covalent geometry for an ICFF model is obtained by Cartesian optimization of the molecule with the local part of the "source" potential function, as defined in the Torsion Energy Term Section. In the absence of nonbonded interactions the geometry optimization of the whole molecule is essentially divided into many virtually independent small optimization problems. To generate a unique initial guess for this geometry optimization we set the torsion angles to the corresponding best energy values found in the calculation of torsion energy profiles. This choice of starting configuration and simplified energy function guarantees fast convergence of the Cartesian optimization, quickly producing a unique covalent geometry for ICFF models. Moreover, such geometries used with ICFF energy functions also result in a better average accuracy in our tests (see Results and Discussion) than more computationally expensive covalent geometries, including those optimized with the full Cartesian force field or MP2 quantum mechanics.

The "local energy" geometry optimization procedure can also be applied to the construction of the ICFF "residue library" for peptides, nucleic acids, or other biopolymers by storing and reusing optimal covalent geometries for the biopolymer residues.

## MMFF94s as a Source Cartesian Force Field

The calculation of ICFF parameters is implemented here on the basis of the MMFF94 Cartesian force field. <sup>28–34</sup> Although the ICFF algorithm formulated above is generic and can be applied to other Cartesian force fields in a similar way, we had several reasons to test the approach with MMFF94. First, it is a "general-purpose" force field parameterized in self-consistent manner for a wide class of organic molecules. <sup>28</sup> Second, MMFF94 functional form, parameters, universal atom type assignment protocol, and molecular validation suite are described in great detail, <sup>28–34</sup> allowing an exact implementation of the force field. Third, an extensive benchmark of MMFF94 with a diverse set of small organic molecules was published recently. <sup>33</sup> Detailed results of these calculations are available in electronic format, and can be used as a reference point for our analysis.

The full MMFF94 potential energy function  $E_{\rm MMFF}$  can be written as

$$E_{\text{MMFF}} = \sum EB_{ij} + \sum EA_{ijk} + \sum EBA_{ijk} + \sum EOOP_{ijk;l} + \sum ET_{iikl} + \sum EvdW_{ij} + \sum EQ_{ij},$$

with bond stretching, bond bending, stretch-bend, out-of plane bending, torsion, van der Waals, and electrostatic terms, respectively. Definitions of these terms are given in eqs. (1)–(14) of ref. 28. Here, we use a version of MMFF94 parameter set, MMFF94s, developed specifically for molecular mechanics applications.<sup>34</sup> Both MMFF94 and MMFF94s parameters are available for public use, and are implemented in a number of molecular modeling packages.<sup>28,34</sup>

## **Computational Methods**

We used the ICM molecular modeling package from Molsoft (www.molsoft.com) for the development the of ICFF torsion coordinate force field and the analysis of ICFF accuracy.

The MMFF94s Cartesian force field energy function, atom type assignment procedure, and parameterization are implemented in the ICM code according to the published guidelines. MMFF94s implementation was verified with a validation suite available at the Computational Chemistry List Web site. 63

All the Cartesian optimizations in this work are performed with the quasi-Newton method in ICM, using the MMFF94s energy function with analytically calculated first derivatives. To constrain a torsion angle  $\theta$  in the torsion energy profile calculations, an effective potential  $(E_r = C_r(\theta - \theta^0)^2, C_r = 10000 \text{ kcal})$  was applied, so that deviation of the angle from the specified value  $\theta^\circ$  usually did not exceed 0.1 degree.

The ICFF module in ICM consists of four major parts. The first procedure assigns an ICFF type for a torsion fragment of the molecule, taking into account the MMFF94 types of the torsion fragment atoms and their connection topology. If torsion parameters already exist for this ICFF type either in the current ICM session or in a saved parameter file, these precomputed parameters are used. Otherwise, the second procedure calculates torsion parameters  $\{A_k, B_k\}$  for the new ICFF torsion type, as described in the ICFF Formulation. The third routine generates covalent geometry by Cartesian optimization of the molecule with the local part of the MMFF94s energy function. An ICM procedure subsequently builds a molecular tree for internal coordinates modeling and fixes bond lengths, bond angles, and phase angles, leaving only torsion coordinates free.<sup>64</sup> Finally, the fourth function calculates the ICFF torsion, "1-5,1-6" interaction, van der Waals, and electrostatic terms with the corresponding first derivatives in internal coordinates. This ICFF energy function can be used in the ICM local minimization or Monte Carlo global energy minimization procedure.65

## Results and Discussion

The accuracy of the Cartesian force field projection onto the torsion coordinate space was assessed by comparing ICFF conformational energies with MMFF94s, both with the original Cartesian MMFF94s<sup>28–34</sup> and an artificial version of MMFF94s with fixed geometry. The latter comparison is used here to illustrate general problems, common for any "flexible" Cartesian energy function when it is applied to a rigid geometry model. We performed two complementary accuracy benchmarks—first using equilibrium conformational energy differences only, and second involving the calculation of detailed torsion energy profiles for the molecules of

**Table 1.** Summary of Equilibrium Energy Comparisons of the Torsion Coordinate Models vs. MMFF94s (see Supplementary Material Table SM-I for Individual ICFF Results), and Force Fields vs. MP4SDQ/TZP.

|                                       | Rigid models vs. Cartesian (MMFF94s <sup>a</sup> ) |                          | Force field models vs. quantum mechanics (MP4SDQ/TZP <sup>b</sup> ) |                      |                  |
|---------------------------------------|----------------------------------------------------|--------------------------|---------------------------------------------------------------------|----------------------|------------------|
|                                       | ICFF                                               | MMFF94rigid <sup>c</sup> | ICFF                                                                | MMFF94s <sup>b</sup> | MM3 <sup>b</sup> |
| Number of                             |                                                    |                          |                                                                     |                      |                  |
| comparisons                           | 120                                                | 118                      | 120                                                                 | 121 <sup>d</sup>     | 108              |
| Number of conformers,                 |                                                    |                          |                                                                     |                      |                  |
| having no local                       |                                                    |                          |                                                                     |                      |                  |
| minimum                               | 1                                                  | 3                        | 1                                                                   | 0                    | 13               |
| Geometry                              |                                                    |                          |                                                                     |                      |                  |
| Mean RMSD, Å                          | 0.04                                               | 0.06                     | 0.04                                                                | 0.04                 | _                |
| Max RMSD, Å                           | 0.27                                               | 0.74                     | 0.48                                                                | 0.31                 | _                |
| Number of favoring                    |                                                    |                          |                                                                     |                      |                  |
| wrong conformer                       | 7                                                  | 8                        | 6                                                                   | 6                    | 18               |
| Max ener. deviation                   | 4.01                                               | 8.12                     | 4.2                                                                 | 1.56                 | 3.47             |
| 1–2 kcal/mol                          | 12                                                 | 14                       | 14                                                                  | 1                    | 22               |
| 2–3 kcal/mol                          | 0                                                  | 11                       | 2                                                                   | 0                    | 5                |
| 3–5 kcal/mol                          | 1                                                  | 1                        | 1                                                                   | 0                    | 0                |
| Energy RMSD, kcal<br>RMSD, kcal, when | 0.64                                               | 1.21                     | 0.74                                                                | 0.31                 | 1.03             |
| $E_{\rm MMFF94s} < 4 \text{ kcal}$    | 0.46                                               | 0.92                     | 0.57                                                                | 0.30                 | 0.91             |

<sup>&</sup>lt;sup>a</sup>MMFF94s data compiled from Table SM-I of ref. 33.

the set. Both approaches use a large diverse set of organic molecules 33,66 compiled by Thomas Halgren to compare the performance of MMFF94 and other major force fields. The only conformational pairs we excluded from the validation set are those with different ring configurations, because ring movements in ICFF have to be treated by a different algorithm. 41,67

# **Equilibrium Energy Comparisons**

The ability of ICFF to reproduce optimal geometries and energy differences between the equilibrium conformers of small organic molecules is summarized in Table 1, while detailed ICFF results are presented in Supplementary Material Table SM-I.

Fixed covalent geometries for each molecule were defined by the local MMFF94s Cartesian optimization as described in the ICFF Formulation. To obtain different equilibrium conformers, we set torsion angles of these molecular models to the corresponding values in the published MP2/6-31G\* geometries<sup>33,66</sup> and then optimized the torsion angles with the ICFF energy function. Note, that it would be inappropriate to use MP2/6-31G\* conformations directly as starting conformations for ICFF optimization, because they have slightly different covalent geometry. ICFF torsion energy optimizations were performed first with some angular constraints to ensure that the molecule would remain in a particular local energy minimum, and then optimized again without these constraints.

Conformations optimized with the ICFF torsion energy function accurately reproduce equilibrium MMFF94s geometries in the conformational set, as measured by the average distance RMSD of about 0.04 Å. Some geometry deviations (maximum RMSD = 0.27 Å) occur in molecules with tri-coordinate nitrogens, where ICFF enforces planar configuration for these nitrogen groups, although this approximation does not significantly affect the accuracy of the ICFF energy measurement.

In Table 1, columns 1 and 2 show that the ICFF potential reproduces the relative equilibrium energies much better than the MMFF94s potential in the rigid geometry model, with RMS errors of  $\sim$ 0.64 kcal versus  $\sim$ 1.21 kcal, respectively. The conformational minimum was not found by ICFF in only one case when the flat nitrogen approximation hindered a shallow conformational minimum for the alanine dipeptide (compare Fig. 4d below). In the other case the error in energy calculations exceeded 2 kcal for the alanine dipeptide analog. As expected, ICFF is more accurate for conformations with low to moderate conformational stress. Thus, the energy RMSD is only 0.46 kcal when eight data points with high MMFF94s energy ( $\Delta E_{\rm MMFF94s} > 4$  kcal) are excluded. Also, one can see that the value of  $RMSD_{(ICFF\ vs.\ QM)} = 0.74\ kcal$  (Table 1, column 3) is close to the linear combination of  $RMSD_{(ICFF \ vs. \ MMFF94s)} = 0.64 \ kcal \ and \ RMSD_{(MMFF94s \ vs. \ QM)} =$ 0.31 kcal (columns 1 and 4), suggesting independent contributions of ICFF approximations into the overall error.

Let us consider the quality of the ICFF projection in the context of Cartesian force field's performance for this molecular set. Note here that the extremely low MMFF94s RMSD = 0.31 kcal in the third column is specific for this molecular set, used for the deri-

<sup>&</sup>lt;sup>b</sup>MP4SDQ/TZP, MMFF94s and MM3 results from Table SM-I of ref. 33.

<sup>&</sup>lt;sup>c</sup>Fixed geometry model with MMFF94s energy function.

<sup>&</sup>lt;sup>d</sup>All conformational pairs without ring conformational changes from ref. 33.

vation of MMFF94.33 Instead, an appropriate unbiased measure of Cartesian force fields accuracy for this molecular set is the second best RMSD (1.03 kcal) reported for the MM3 force field33 (column 4). Given a Cartesian force field RMSD 1.03 kcal, and assuming that the ICFF projection yields an independent RMSD contribution of 0.64 kcal, one can estimate that ICFF projection of a Cartesian force field into the torsion coordinate space adds only about 18% (( 1.032 0.642 )/1.03 - 1.18) to the overall energy RMSD.

# *Torsion Profile Comparisons*

Equilibrium conformation analysis, however, gives little information about the shape of the energy function, which is particularly important in the most statistically populated lower energy part of the conformational profile. To study this aspect of the ICFF performance we carried out a detailed analysis of the torsion term alone, and then the full ICFF conformational energy function, by comparing ICFF and MMFF94s torsion profiles. All 110 dihedral angles responsible for conformational changes in the above molecular set (Table SM-I) were chosen for the energy profiling.

## *Accuracy of the ICFF Torsion Term*

We defined a "local energy" as a part of the full conformational energy without electrostatics, van der Waals (beyond "1-4" atom pairs) and any other interactions that involve atoms, separated by four or more covalent bonds. Thus, the ICFF local energy here is represented by the torsion term alone, while the MMFF94s local energy is comprised of torsion, bond bending, bond stretching, out-of-plane and "1-4" van der Waals terms.

Figure 2 displays examples of typical local energy profiles calculated in terms of the ICFF rigid (dashed), MMFF94s rigid (dotted) and MMFF94s Cartesian (solid lines) models. The first two plots in Figure 2 represent a most abundant group (75 out of 110 molecules) of high accuracy torsion profiles with the energy RMSD between ICFF and MMFF94s Cartesian models less then 0.2 kcal. Although agreement between ICFF and the "source" force field is almost perfect in this group, the rigid model with the Cartesian potential often fails even in these relatively simple cases. Thus, Figure 2a shows that the flexibility of NOCOCON bonds in the glycine dipeptide analog allows a local minimum in the *C*5 conformation when accounted for explicitly in the Cartesian model, or implicitly in ICFF. The rigid model with the MMFF94s energy function does not reproduce this minimum at all. Similarly, the second example (Fig. 2b) shows a perfect agreement between ICFF and MMFF94s models, but strong qualitative deviations for the rigid MMFF94s model. Good overall agreement between ICFF and MMFF94s is largely conserved for the cases with higher RMSD, like in Figure 2c and d, and even for the two cases with the worst RMSD in the set (Fig. 2e and f) (see Supplementary material Fig. SM1 for all the 110 profiles).

Note that here we calculate RMSD for those values of the torsion angle where the MMFF94s local energy is lower than a certain energy cutoff, *E*MMFF94s *E*cut. Because we are not interested in the quantitative representation of high-energy rotational barriers in ICFF applications, we can focus the benchmarking on the range of conformations with low to moderate energy stress. Moreover, we have to exclude highly stressed conformations from the RMSD analysis anyway, because the reference MMFF94s calculations also lose accuracy when torsion deformations are high. The energy cutoff *E*cut - 5 kcal is applied in Figure 2, with datapoints below *E*cut shown in open circles.

Table 2 summarizes the ICFF torsion energy term accuracy and compares it to the accuracy of the local energy function in the rigid MMFF94s model. The ICFF implicit flexibility model reproduces the reference MMFF94s local conformational energies much better than the rigid MMFF94s model. For example, the RMSD of all ICFF energy comparisons in the moderate energy range (*E*MMFF94s 5 kcal) is only 0.26 kcal, vs. RMSD - 0.58 kcal for the rigid model with the Cartesian energy function.

Although the implicit flexibility incorporated in ICFF torsion parameters effectively reduces inappropriate energy deviations in the rigid geometry model, these errors cannot be eliminated completely. These residual errors come from the fact that we are limited to a one-dimensional definition of the torsion energy function. Although relative positions of a torsion fragment atoms depend on a single torsion variable in the rigid model, their interactions with neighbor atoms outside of the torsion fragment depend on two or more torsion variables. These "outside" atoms can indirectly influence interactions between atoms of the torsion fragment through the covalent geometry of the torsion fragment. In some rare cases this indirect influence can amount to about a 1-kcal energy deviation if the torsion fragment bonds are highly flexible and an atom inside the torsion fragment is severely clashed with an atom outside of the fragment (e.g., see Fig. 2f). On average, though, the error due to this effect is rather small, and apparently does not call for such a major complication as the introduction of torsion potentials that depend on two or more torsion variables.

## *Fourier Series Approximation*

We found that a sixfold Fourier series function is sufficient to reproduce torsion profiles for the test set with an accuracy of 0.01 kcal. Although the performance of a fivefold series (RMSD 0.08 kcal) and a fourfold series (RMSD 0.2) is still acceptable, we can easily afford some excessive accuracy at this step. Because small Fourier coefficients (e.g., less than 0.01) can be nullified, these extra terms do not slow down energy evaluations noticeably.

Moreover, we suggest fitting Fourier coefficients with a weighted least square procedure, where the weight of each data point is proportional to the inverse relative energy, W() [*E*() Min(*E*()) 1]<sup>1</sup> . This fitting procedure puts the emphasis on the low-energy part of the torsion profile improving RMSD of fitting to 0.001 kcal for energies below 5-kcal cutoff. At the same time, the algorithm takes into account a limited accuracy of the source force field for high-energy conformations and reduces the effect of occasional deviations in high-energy regions on the resulting Furrier coefficients.

## *Accuracy of the Full ICFF Energy Function*

The adequate treatment of interactions between atoms, separated by four or five covalent bonds ("1-5, 1-6" interactions) is one of the major challenges for the torsion coordinate modeling. Unlike

![](_page_6_Figure_2.jpeg)

**Figure 2.** Six examples of "local" conformational energy torsion profiles: (a) glycine dipeptide analog (NOCOCON); (b) alanine dipeptide analog (NOCOCON); (c) methyl isopropyl ether (COOOCOC); (d) isopropyl formate (COCOOOC); (e) NOOH, NOMe propionamide (COCONOC); (f) vinyl formate (COOOCOH). Conformational energies are calculated in the rigid geometry model with ICFF (dashes lines) and MMFF94s (dotted lines) "local" energy functions at different values of the torsion angle. For the Cartesian model calculations (solid lines), bond lengths and bond angles for each value of the torsion angle are relaxed in Cartesian coordinates with the MMFF94s local energy. Marked with circles are data points included in the calculation of the RMSD (*E*MMFF94s(local) 5 kcal). Torsion profiles are ranked from 1 to 110 according to their ICFF RMSD values (see Supplementary Material Fig. SM-1 for all 110 profiles).

**Table 2.** The Accuracy of Two Rigid Geometry Models with Different Local Energy Functions, Compared to the Local Energy in the MMFF94s Cartesian Model.

| comparisons <sup>a</sup> | RMSD, kcal                                   | RMSD, kcal                                                                 |
|--------------------------|----------------------------------------------|----------------------------------------------------------------------------|
| 1010                     | 0.12                                         | 0.23                                                                       |
| 2125                     | 0.20                                         | 0.42                                                                       |
| 2651                     | 0.26                                         | 0.58                                                                       |
| 2841                     | 0.30                                         | 0.66                                                                       |
| 2956                     | 0.33                                         | 0.80                                                                       |
| 3091                     | 0.42                                         | 1.23                                                                       |
| 3190                     | 0.53                                         | 1.53                                                                       |
|                          | 1010<br>2125<br>2651<br>2841<br>2956<br>3091 | 1010 0.12<br>2125 0.20<br>2651 0.26<br>2841 0.30<br>2956 0.33<br>3091 0.42 |

<sup>&</sup>lt;sup>a</sup>All energy comparisons with relative MMFF94s energy  $E_{\rm MMFF94s} < E_{\rm cut}$  <sup>b</sup>Rigid geometry model with MMFF94s local energy function.

"1-4" contacts, "1-5, 1-6" interactions cannot be included in the torsion term and have to be calculated explicitly. At the same time using the original MMFF94 van der Waal term in fixed geometry models is not appropriate. We found many examples in our molecular set where "1-5, 1-6" van der Waals clashes would lead to an improper balance between two conformers and large distortions of equilibrium geometry. Figure 3 displays an example of a 1-5 clash in the C=C-O-C trans configuration of vinyl formate. The original MMFF94s van der Waals energy term in the rigid model reaches 24 kcal for the O-H contact (black conformation in Fig. 3a and solid line in Fig. 3b). In contrast, the energy never exceeds 1 kcal when the bonded geometry is relaxed in this conformation with the torsion angles constrained (gray conformation in Fig. 3a and dotted line in Fig. 3b). To mimic such bond flexibility in the ICFF rigid geometry model, we introduce a soft polynomial term for "1-5, 1-6" repulsions as described in the ICFF Formulation section. Calculated with this empirical term, the energy of this van der Waals contact is much closer to the flexible model values (black conformation in Fig. 3a and dashed line in Fig. 3b).

Figure 4 presents a set of six typical examples of the conformational energy torsion landscapes calculated for MMFF94s Cartesian models (solid lines) and ICFF models with implicit flexibility (dashed lines), as well as for ICFF models with the original MMFF94 van der Waals term (dotted line). Figure 4a and b display two typical results in a group of 79 high-accuracy profiles with an ICFF vs. MMFF94s energy RMSD < 0.4 kcal. Replacing the "hard" van der Waals repulsion with a soft cubic term dramatically improves agreement between ICFF and Cartesian models in a variety of cases shown in Figure 4a–f, as well as for the whole molecular set (see Supplementary material Fig. SM-2).

The coulomb electrostatics between 1-5 and 1-6 neighbors as a part of the full conformational energy does not introduce substantial errors in the rigid geometry model approximation, at least for the molecular set studied, and therefore, it does not require any special treatment.

Table 3 summarizes the accuracy of full energy torsion profiles for the rigid geometry model with different "1-5,1-6" repulsion functions. The original MMFF94s van der Waals term in the rigid geometry model (column 5) gives very large RMSDs even in the

low energy range. This RMSDs can be apparently reduced by eliminating the van der Waals term for the "1-5,1-6" interactions completely (column 4), though such a "simplification" would be unacceptable in real molecular mechanics simulations.

When the "1-5,1-6" repulsion is calculated with the ICFF cubic polynomial function (column 3) the relative energies for the torsion set are very close to the reference MMFF94s values in the low and moderate energy range. This accuracy was achieved by fitting the only adjustable factor C in the ICFF repulsion term to minimize the RMSD at three different values of energy cutoff  $E_{\rm cut}=3,5$  and 10 kcal. Yet another fitting of the C factor was performed for the equilibrium conformational energy set, described above. For all four objective functions, the optimal RMSD is achieved in the range of 0.50 < C < 0.60, and all four show a very weak dependence on C factor in this range. Thus, a value C=0.55 was chosen as the consensus for the four optimizations, yielding RMSD = 0.37 kcal at the 5-kcal MMFF94s energy cutoff.

In general, the analysis of our results for "1-5,1-6" repulsion lead to the following conclusions. First, the bond relaxation dra-

![](_page_7_Figure_12.jpeg)

Figure 3. Steric clash between vinyl formate oxygen and hydrogen separated by five covalent bonds ("1-6" interaction). (a) Geometry of the molecule in C—C—O—C cis conformation in the rigid geometry model (black lines) and with covalent geometry relaxed by Cartesian minimization (gray lines). (b) Repulsion energy between the two atoms in the rigid model with the MMFF94s van der Waals repulsion term is shown by the solid line; the rigid ICFF model with the soft repulsion [see eq. (2)] is represented by the dashed line; and the MMFF94s Cartesian model is shown by the dotted line. We varied the distance between O and H by changing the torsion angle,  $\Omega$ .

![](_page_8_Figure_2.jpeg)

**Figure 4.** Six examples of full conformational energy torsion profiles: (a) propenoic acid (HOCOCOO); (b) 2-methylpropenamide (COCOCON); (c) isopropyl formate (COCOOOC); (d) alanine dipeptide (COCONOH); (e) 2-methyl-but-1-ene-3-one (COCOCOC); (f) 4-oxobutanal [HOCOCOC(AO)]. Conformational energy with soft empirical (dashes lines) and MMFF94s van der Waals (dotted lines) "1-5,1-6" repulsion potential is calculated at different values of the torsion angle in the ICFF rigid geometry model. For the Cartesian MMFF94s model calculations (solid lines), bond lengths and bond angles for each value of the torsion angle are relaxed in Cartesian coordinates. Marked with circles are data points included in the calculation of the RMSD (*E*MMFF94s 5 kcal). Torsion profiles are ranked from 1 to 110 according to their ICFF RMSD values (see Supplementary Material Fig. SM-2 for all 110 profiles).

**Table 3.** The Accuracy of the ICFF Models with Different "1-5,1-6" Repulsion Functions, Compared to the MMFF94s Cartesian Model.

| $E_{\rm cut}$ , <sup>a</sup> kcal | Number of comparisons <sup>a</sup> | ICFF<br>RMSD, kcal       | ICFFno1-56 <sup>b</sup><br>RMSD, kcal | ICFFvw1-56°<br>RMSD, kcal |
|-----------------------------------|------------------------------------|--------------------------|---------------------------------------|---------------------------|
| 1.                                | 796                                | 0.21                     | 0.33                                  | 0.87                      |
| 3.                                | 1869                               | <b>0.31</b> <sup>d</sup> | 0.54                                  | 1.67                      |
| 5.                                | 2414                               | 0.37                     | 0.62                                  | 3.49                      |
| 7.                                | 2652                               | 0.45                     | 0.73                                  | 12.37                     |
| 10.                               | 2851                               | 0.56                     | 0.87                                  | 17.73                     |
| 20.                               | 3079                               | 1.02                     | 1.22                                  | 61.13                     |
| All                               | 3190                               | 1.40                     | 1.31                                  | 60.59                     |

 $<sup>^{\</sup>rm a}{\rm All}$  energy comparisons with relative MMFF94s energy  $E_{\rm MMFF94s} < E_{\rm cut}.$ 

matically reduces the conformational energy of "1-5,1-6" clashes in the flexible model, so it must be accounted for implicitly in the rigid model. Second, the functional form we have chosen for the empirical "1-5,1-6" repulsion potential mimics the bond flexibility reasonably well, while utilizing the same parameters  $\varepsilon_{IJ}$ ,  $R_{IJ}^*$ , and  $R_{IJ}^0$  as the original van der Waals term. Let us note here, that this energy function mimics only an *average* flexibility of four or five covalent bonds. A more elaborated "1-5,1-6" repulsion potential may be considered, which implicitly accounts for individual bond rigidity parameters, although it is impractical at this accuracy level.

Note also that while "1-5,1-6" repulsion in ICFF is calculated with a "soft" potential to compensate for improper steric clashes, the interactions between atoms, connected by six or more bonds are treated with the original MMFF94 van der Waals term. This discrimination is based on the assumption that six-bond (and longer) connections have enough freedom in torsion coordinates to easily avoid steric hindrances in torsion models. This assumption

is correct for most chemical structures and for all the molecules in our test suite. In some rare cases though (e.g., in conjugated ring systems) all the torsions can be very stiff even between atoms separated by six or seven bonds, and the ICFF soft repulsion potential would be more appropriate for the interactions between these atoms.

## Effect of the Covalent Geometry on ICFF Accuracy

Presented in Table 4 are results illustrating the accuracy of ICFF models with different covalent geometries, generated by minimizations of "local" and full MMFF94 energy functions, as well as by quantum mechanical (MP2) geometry optimization. All three methods yield acceptable covalent geometries for ICFF models, with the difference between corresponding energy RMSD values not exceeding 25% in the whole range of the energy cutoffs. This relatively low sensitivity to the bonded geometry variations is an

**Table 4.** The Accuracy of the ICFF Models with Different Covalent Geometries, Compared to the MMFF94s Cartesian Model.

| $E_{\rm cut}^{a}$ kcal | Number of comparisons <sup>a</sup> | ICFF <sup>b</sup><br>RMSD, kcal | ICFF//MMFF94s_full <sup>a</sup><br>RMSD, kcal | ICFF//MP2 <sup>d</sup><br>RMSD, kcal |
|------------------------|------------------------------------|---------------------------------|-----------------------------------------------|--------------------------------------|
| 1.                     | 796                                | 0.21                            | 0.23                                          | 0.23                                 |
| 3.                     | 1869                               | 0.31                            | 0.34                                          | 0.36                                 |
| 5.                     | 2414                               | 0.37                            | 0.43                                          | 0.49                                 |
| 7.                     | 2652                               | 0.45                            | 0.51                                          | 0.56                                 |
| 10.                    | 2851                               | 0.56                            | 0.66                                          | 0.63                                 |
| 20.                    | 3079                               | 1.02                            | 1.10                                          | 0.99                                 |
| All                    | 3190                               | 1.40                            | 1.55                                          | 1.23                                 |

 $<sup>^{\</sup>rm a}{\rm All}$  energy comparisons with relative MMFF94s energy  $E_{\rm MMFF94s} < E_{\rm cut}$ 

<sup>&</sup>lt;sup>b</sup>ICFF model without "1-5,1-6" repulsion.

cICFF model with "1-5,1-6" repulsion calculated with MMFF94s van der Waals term.

<sup>&</sup>lt;sup>d</sup>Shown in bold are data used in the fitting of factor C.

<sup>&</sup>lt;sup>b</sup>Covalent geometries for ICFF models obtained by Cartesian optimization of the "local" part of the MMFF94s energy function.

<sup>&</sup>lt;sup>c</sup>Covalent geometries obtained with the full MMFF94s energy function, including van der Waals and electrostatic terms.

<sup>&</sup>lt;sup>d</sup>Covalent geometries obtained by MP2 geometry optimization.<sup>33</sup>

intrinsic feature of the ICFF energy function, achieved by the special treatment of "1-4" and "1-5,1-6" interactions.

Still, the results clearly favor the "local energy" Cartesian optimization, where the local energy includes only bonded terms and "1-4" van der Waals interactions, as the best method to generate covalent geometry for ICFF models. The better accuracy of this approach can be rationalized by the following considerations. The bond lengths and bond angles are always optimized at a specific conformational minimum of the molecule; the freezing of the covalent geometry in this minimum results the model becoming somewhat biased towards this conformation. When the covalent geometry is optimized with a full Cartesian force field, specific electrostatic and van der Waals interactions would significantly distort the result. Omitting these terms from the optimization procedure makes the resulting covalent geometry more "generic," i.e., less biased. (The other possibility to reduce this bias could be the "averaging" of covalent geometries over all torsional configurations of the molecule, but this procedure would be too computationally expensive and does not guarantee better results.)

The other advantage of using "local energy" in the covalent geometry optimizations is the fast and reliable convergence of the procedure even for large organic molecules. The Cartesian minimization problem for the molecule is essentially divided into simple local minimizations, almost independent from each other, and can easily be solved. Choosing starting conformations with torsion variables corresponding to the best energy minima in the torsion profiles is a smart initial guess for this optimization, which makes the procedure deterministic, i.e., it assures a resulting unique and predictable covalent geometry.

## *ICFF Performance*

The ICFF algorithm of a Cartesian force field projection into torsion coordinate space does not require expensive calculations and can be automated, and parameters for the torsion force field can be generated "on the fly." The method, implemented within the ICM molecular modeling package (MolSoft), precalculates fixed geometry and torsion parameters for an average drug-like compound on a 0.1 second scale. The same procedure can be used to generate geometry and torsion parameter libraries for proteins and other biopolymers. This opens a new possibility to the use of a consistent and widely transferable potential for torsion coordinate modeling of a variety of compounds and biopolymers.

The speed performance of ICFF is improved by about 50% compared to the other torsion force field, ECEPP. This is a result of a reduced number of pairwise interactions, due to the implicit treatment of "1-4" van der Waals term and simplified "1-5,1-6" repulsion term.

# **Conclusions**

Using the ICFF approach, a conformational energy function for the fixed covalent geometry model can be derived from a source Cartesian force field directly. Using MMFF94s as a source, we show that its ICFF projection into torsion coordinates closely mimics the Cartesian force field with RMSD of about 0.3– 0.5 kcal. Although the quality of the derived torsion coordinate force field is always limited by the quality of the source energy function, it can approach this limit rather closely.

Our data demonstrates that "implicit flexibility" built into the internal coordinate parameters is critical for the performance of the torsion coordinate model. In a representative test, based on the comparison of conformational energy profiles for a diverse molecular set, the energy RMSD is reduced from 3.5 kcal in the rigid model to 0.37 kcal in the ICFF model. Flexibility is incorporated into ICFF in two ways. First, the composite torsion term is parametrized to account implicitly for bond stretching, bond bending, out-of-plane, and 1-4 van der Waals interactions. Second, an empirical soft potential for "1-5,1-6" repulsion is derived to mimic the flexibility of covalent bonds, connecting these atoms. Fixed covalent geometry for the ICFF models can be obtained by fast and deterministic geometry minimizations with the "local" part of the source Cartesian force filed.

# **Acknowledgments**

The authors thank Dr. Alexander MacKerell for a number of useful suggestions.

# **References**

- 1. Friesner, R. A.; Gunn, J. R. Ann Rev Biophys Biomol Struct 1996, 25, 315.
- 2. Lazaridis, T.; Karplus, M. Curr Opin Struct Biol 2000, 10, 139.
- 3. Hao, M. H.; Scheraga, H. A. Curr Opin Struct Biol 1999, 9, 184.
- 4. Abagyan, R.; Batalov, S.; Cardozo, T.; Totrov, M.; Webber, J.; Zhou, Y. Y. Proteins 1997, Supplement 1, 29.
- 5. Maiorov, V.; Abagyan, R. Proteins Struct Funct Genet 1997, 27, 410.
- 6. Kitao, A.; Go, N. Curr Opin Struct Biol 1999, 9, 164.
- 7. Lengauer, T.; Rarey, M. Curr Opin Struct Biol 1996, 6, 402.
- 8. Sternberg, M. J. E.; Gabb, H. A.; Jackson, R. M. Curr Opin Struct Biol 1998, 8, 250.
- 9. Elcock, A. H.; Sept, D.; McCammon, J. A. J Phys Chem B 2001, 105, 1504.
- 10. Totrov, M.; Abagyan, R. Protein-Ligand Docking as an Energy Optimization Problem; Raffa, R. B., Ed.; John Wiley & Sons, Ltd: London, 2001; p 603, vol. 1.
- 11. Kramer, B.; Metz, G.; Rarey, M.; Lengauer, T. Med Chem Res 1999, 9, 463.
- 12. Abagyan, R.; Totrov, M. Curr Opin Struct Biol 2001, 5, 375.
- 13. Abagyan, R. A.; Totrov, M. M. Ame Chem Soc Abstr 1996, 211, 35-COMP.
- 14. Halgren, T. A. Curr Opin Struct Biol 1995, 5, 205.
- 15. Brooks, B. R.; Bruccoleri, R. E.; Olafson, B. D.; States, D. J.; Swaminathan, S.; Karplus, M. J Comput Chem 1983, 4, 187.
- 16. MacKerell, A. D.; Bashford, D.; Bellott, M.; Dunbrack, R. L.; Evanseck, J. D.; Field, M. J.; Fischer, S.; Gao, J.; Guo, H.; Ha, S.; Joseph–McCarthy, D.; Kuchnir, L.; Kuczera, K.; Lau, F. T. K.; Mattos, C.; Michnick, S.; Ngo, T.; Nguyen, D. T.; Prodhom, B.; Reiher, W. E.; Roux, B.; Schlenkrich, M.; Smith, J. C.; Stote, R.; Straub, J.; Watanabe, M.; Wiorkiewicz–Kuczera, J.; Yin, D.; Karplus, M. J Phys Chem B 1998, 102, 3586.
- 17. MacKerell, A. D. J.; Banavali, N. K. J Comput Chem 2000, 21, 105.
- 18. Allinger, N. L.; Zhou, X. F.; Bergsma, J. J Mol Struct (Theochem) 1994, 118, 69.
- 19. Allinger, N. L.; Durkin, K. A. J Comput Chem 2000, 21, 1229.

- 20. Chen, K. H.; Walker, G. A.; Allinger, N. L. J Mol Struct (Theochem) 1999, 490, 87.
- 21. Weiner, S. J.; Kollman, P. A.; Case, D. A.; Chandra Singth, U.; Ghio, C.; Alagona, G.; Prefeta, S. J.; Wiener, P. J Am Chem Soc 1984, 106, 765.
- 22. Weiner, S. J.; Kollman, P. A.; Nguyen, D. T.; Case, D. A. J Comput Chem 1986, 7, 230.
- 23. Pearlman, D. A.; Case, D. A.; Caldwell, J. W.; Ross, W. S.; Cheatham, T. E.; Debolt, S.; Ferguson, D.; Seibel, G.; Kollman, P. Comput Phys Commun 1995, 91, 1.
- 24. Maple, J. R.; Hwang, M. J.; Jalkanen, K. J.; Stockfisch, T. P.; Hagler, A. T. J Comp Chem 1998, 19, 430.
- 25. Halgren, T. A.; Murphy, R. B.; Jorgensen, W. L.; Friesner, R. A. Am Chem Soc Abstr 2000, 220, 2-COMP.
- 26. Kaminski, G. A.; Friesner, R. A.; Tirado–Rives, J.; Jorgensen, W. L. Am Chem Soc Abstr 2000, 220, 14-COMP.
- 27. Jorgensen, W. L.; Tirado–Rives, J. Am Chem Soc Abstr 1998, 216, 043-COMP.
- 28. Halgren, T. A. J Comput Chem 1996, 17, 490.
- 29. Halgren, T. A. J Comput Chem 1996, 17, 520.
- 30. Halgren, T. A. J Comput Chem 1996, 17, 553.
- 31. Halgren, T. A. J Comput Chem 1996, 17, 616.
- 32. Halgren, T. A.; Nachbar, R. B. J Comput Chem 1996, 17, 587.
- 33. Halgren, T. A. J Comput Chem 1999, 20, 730.
- 34. Halgren, T. A. J Comput Chem 1999, 20, 720.
- 35. Cheng, A.; Best, S. A.; Merz, K. M.; Reynolds, C. H. J Mol Graphics Model 2000, 18, 273.
- 36. Jorgensen, W. L.; Maxwell, D. S.; Tirado–Rives, J. J Am Chem Soc 1996, 118, 11225.
- 37. Coleman, R. S.; McCary, J. L. Bioorg Med Chem Lett 1998, 8, 3039.
- 38. Foloppe, N. D.; MacKerell, A Jr. J Comput Chem 2000, 21, 86.
- 39. Feller, S. E.; MacKerell, A. D. J Phys Chem B 2000, 104, 7510.
- 40. Abagyan, R. A.; Mazur, A. K. J Biomol Struct Dyn 1989, 6, 833.
- 41. Mazur, A. K.; Abagyan, R. A. J Biomol Struct Dyn 1989, 6, 815.
- 42. Horiuchi, T.; Go, N. Proteins Struct Funct Genet 1991, 10, 106.
- 43. Vasquez, M. N. G.; Scheraga, H. A. Chem Rev 1994, 94, 2183.
- 44. Abagyan, R.; Totrov, M.; Kuznetsov, D. J Comput Chem 1994, 15, 488.
- 45. Lee, B.; Kurochkina, N.; Kang, H. S. FASEB J 1996, 10, 119.
- 46. Schaumann, T.; Braun, W.; Wuthrich, K. Biopolymers 1990, 29, 679.

- 47. Cardozo, T.; Totrov, M.; Abagyan, R. Proteins Struct Func Genet 1995, 23, 403.
- 48. Norledge, B. V.; Lambeir, A. M.; Abagyan, R. A.; Rottmann, A.; Fernandez, A. M.; Filimonov, V. V.; Peter, M. G.; Wierenga, R. K. Proteins Struct Funct Genet 2001, 42, 383.
- 49. Totrov, M.; Abagyan, R. Proteins 1997, Supplement 1, 215.
- 50. Trosset, J. Y.; Scheraga, H. A. J Comput Chem 1999, 20, 244.
- 51. Li, Z.; Scheraga, H. A. Proc Natl Acad Sci USA 1987, 84, 6611.
- 52. Abagyan, R.; Totrov, M. J Mol Biol 1994, 235, 983.
- 53. Abagyan, R. A.; Totrov, M. J Comput Phys 1999, 151, 402.
- 54. Momany, F. A.; McGuire, R. F.; Burgess, A. W.; Scheraga, H. A. J Phys Chem 1975, 79, 2361.
- 55. Nemethy, G. P. M. S.; Scheraga, H. A. J Phys Chem 1983, 87, 1883.
- 56. Nemethy, G.; Gibson, K. D.; Palmer, K. A.; Yoon, C. N.; Paterlini, G.; Zagari, A. J Phys Chem 1992, 96, 6472.
- 57. Strynadka, N. C. J.; Eisenstein, M.; KatchalskiKatzir, E.; Shoichet, B. K.; Kuntz, I. D.; Abagyan, R.; Totrov, M.; Janin, J.; Cherfils, J.; Zimmerman, F.; Olson, A.; Duncan, B.; Rao, M.; Jackson, R.; Sternberg, M.; James, M. N. G. Nat Struct Biol 1996, 3, 233.
- 58. Schapira, M.; Totrov, M.; Abagyan, R. J Mol Recognit 1999, 12, 177.
- 59. Stigler, R. D.; Hoffmann, B.; Abagyan, R.; Schneider–Mergener, J. Struct Fold Design 1999, 7, 663.
- 60. Trosset, J. Y.; Scheraga, H. A. J Comput Chem 1999, 20, 412.
- 61. Hinsen, K.; Kneller, G. R. Phys Rev E 1995, 52, 6868.
- 62. Dunfield, L. G. B. A. W.; Scheraga, H. A. J Phys Chem 1978, 82, 2609.
- 63. http://ccl.net/cca/data/MMFF94/ and http://ccl.net/cca/data/MMFF94s/.
- 64. Abagyan, R.; Frishman, D.; Argos, P. Proteins Struct Funct Genet 1994, 19, 132.
- 65. Abagyan, R.; Totrov, M.; Maiorov, V. Am Chem Soc Abstr 1996, 212, 69-COMP.
- 66. This molecular set is also available in computer-readable format at Computational Chemistry List (CCL) website (http://ccl.net/cca/data/ ff\_evaluation\_suite/).
- 67. Ring deformations involve dramatic changes in bonded geometry, and thus the usual torsion potential approach is not applicable to them in torsion coordinates. Rings can be treated in torsion coordinates according to one of the three models: (i) a single fixed ring geometry, which completely ignores deformations, (ii) a set of alternative fixed geometries and corresponding conformational energies, (iii) continuous ring deformations using pseudotorsion coordinates.