# Overcoming potential energy distortions in constrained internal coordinate molecular dynamics simulations

Saugat Kandel, Romelia Salomon-Ferrer, Adrien B. Larsen, Abhinandan Jain, Nagarajan Vaidehi. J. Chem. Phys. 144, 044112 (2016). doi:10.1063/1.4939532

## Abstract

The Internal Coordinate Molecular Dynamics (ICMD) method describes bonded systems (proteins, polymers) in bond/angle/torsion (BAT) coordinates and enables coarsening of the dynamics model at multiple hierarchical levels (e.g. treating helical domains as rigid bodies, loops as flexible torsions). Using constraints to treat bond lengths and bond angles as rigid distorts the potential energy landscape and reduces the number of dihedral transitions and conformational sampling. This paper presents a two-part solution:

1. To alleviate the **intrinsic distortion** (from the reduced phase space of torsional MD), use the **Fixman compensating potential**.
2. To alleviate the **extrinsic distortion** (from coupling between dihedral angles and bond angles in the forcefield), use a proposed **hybrid ICMD method** that selectively relaxes chosen bond angles.

Together these remove the potential energy distortions in constrained ICMD simulations of peptides. The hybrid ICMD method bridges the gap between all-atom MD and torsional MD.

## I. Introduction

All-atom MD solves Newton's equations in absolute Cartesian coordinates with all DOF movable. It is simple but not suited to coarsening the dynamics model for long-time-scale processes (e.g. domain motion). Internal coordinates provide a natural way to apply holonomic constraints directly to the dynamics model: constraints are applied by simply excluding selected DOF, while retaining ODE (not DAE) form for the equations of motion for non-loop systems. The ICMD model with all bonds and angles frozen and only torsions free is the well-known **torsional MD (TMD)** model.

A longstanding challenge: holonomic constraints alter the equilibrium statistical properties relative to unconstrained all-atom models. These statistical differences affect the probability density functions (pdfs) and the transition-barrier crossing rates between microstates. The distortions stem from two sources:

1. **Mass-matrix-dependent intrinsic distortion:** inherent to the constrained model, arising from the reduced dimension of phase space. The configuration of the frozen DOF affects the mass-matrix determinant and thus the partition function. Independent of the forcefield potentials.
2. **Forcefield-dependent extrinsic distortion:** a consequence of dynamic coupling among BAT DOF from the external forcefield. In proteins, bond lengths are well separated in frequency from bond angles and torsions (negligible coupling), but bond angles can have overlapping frequencies with torsions, so freezing bond angles distorts torsional distributions and alters barrier heights.

The Fixman potential (a configuration-dependent correction) has been studied to quantify how well it corrects the constraint-induced bias, but prior investigations were limited to small idealized serial chains without extrinsic forcefield effects. The authors' prior work developed a computationally tractable spatial-operator-algebra algorithm to compute the Fixman potential for serial chains and general branched systems, showing that without extrinsic coupling, the constrained model plus Fixman potential recovers unconstrained torsional pdfs even for branched peptides.

This paper shows: the Fixman potential completely removes intrinsic bias from configuration pdfs but only partially for velocity-dependent barrier-crossing rates; and applying it to general branched molecules with all-atom forcefields establishes that it alone does not compensate extrinsic distortions. To reduce extrinsic distortion in the configurational pdf, they develop the **hybrid ICMD method**, which lets the user treat any desired bond-angle DOF as flexible while keeping others rigid (Fixman applied by default). Opening a small subset of key backbone angles is sufficient to retrieve the all-atom Cartesian config pdfs with no adverse impact on time-step size.

Advantages: (1) freeze torsions in rigid parts while allowing bond-angle motion in flexible parts - coarse graining without compromising forcefield accuracy; (2) a fully flexible internal-coordinate MD simulation becomes possible.

## II. Methods

The hybrid ICMD simulations build on the previously developed **Generalized Newton-Euler Inverse Mass Operator (GNEIMO)** method. Since some bond lengths/angles can be treated as rigid, the DOF in the equations of motion become coupled and take the form of eq:1. The dynamics is obtained by solving eq:1 for the acceleration $\ddot\alpha$ and integrating. GNEIMO uses spatial operator algebra to derive an analytical inverse of the mass matrix, giving the closed form eq:2 for $\ddot\alpha$; the $H,\psi,\mathcal{K}$ etc. terms are mass-matrix factorizations (Jain refs 21, 26). The right-hand side is evaluated by cost-effective recursive algorithms. These recursive equations are generic and unchanged whether some or all bond angles are open (hybrid ICMD). The GNEIMO-Fixman method for the Fixman torque also continues to apply when bond angles are open. GNEIMO underlies the *GneimoSim* ICMD software package.

### A. Fixman potential to correct for the intrinsic distortion

In TMD all bond lengths and bond angles are rigid, changing the statistics from the trajectories. Fixman proposed a correction potential depending on the mass-matrix determinant. The authors previously showed it recovers the correct config pdfs for simple and complex branched molecules. This section shows it only partially removes the velocity-dependent barrier-crossing-rate bias and, as expected, does not correct forcefield-induced extrinsic distortion.

#### 1. Derivation of the Fixman potential

### Derivation (not implemented)

A polymer with $n$ atoms and $3n$ Cartesian coordinates can equivalently be described in $3n$ BAT coordinates. For a constrained system, the $3n$ BAT coords partition into $N$ unconstrained coords $\alpha$ and $(3n-N)$ constrained coords $q$.

For the unconstrained system, both $\alpha$ and $q$ vary; the canonical momenta $\mathfrak{p}$ give the unconstrained Hamiltonian eq:3. At temperature $T$ the partition function is eq:4. Substituting eq:3 into eq:4 and integrating over momenta gives eq:5, hence the unconstrained configuration pdf eq:6.

Remarkably, $\det\{\mathcal{M}_B\}$ does not depend on the torsions and is a simple product over individual BAT bond lengths and bond angles (eq:7), so it factorizes as $f_1(\alpha)f_2(q)$ (eq:8).

For a constrained model with $q$ frozen at $q_0$, the partition function is eq:9, giving pdf eq:10. In contrast to $\det\{\mathcal{M}_B\}$, the constrained $\det\{\mathcal{M}\}$ does NOT decompose into a product over individual bond lengths/angles; it introduces coupling between torsions, bond lengths, and bond angles in $\rho_c$ - a systematic bias called the **intrinsic distortion** (its source is the constrained model's mass matrix, not the forcefield).

To compensate, Fixman proposed a modified potential $U'(\alpha)=U(\alpha,q_0)+U_f(\alpha)$ with the Fixman potential $U_f$ given by eq:11. Using eq:11 in eq:10 gives the compensated pdf eq:12.

Comparing eq:6 and eq:12, the compensated constrained pdf $\rho_f(\alpha,q_0)$ equals the unconstrained $\rho_u(\alpha,q)$ under any of:
1. No force potential, $U(\alpha,q)=0$ (the case studied by prior Fixman work).
2. Constrained and unconstrained coords are **separable**: $U(\alpha,q)=U_1(\alpha)+U_2(q)$; then $e^{-U/kT}$ factorizes and, using eq:8, both sides factorize -> $\alpha$ and $q$ statistically independent.
3. $U(\alpha,q)$ is very steep for $q$ around $q_0$, so $q$ barely varies: $U(\alpha,q)\approx U(\alpha,q_0)$ and $\det\{\mathcal{M}_B(\alpha,q)\}\approx\det\{\mathcal{M}_B(\alpha,q_0)\}$.

In general these conditions may not hold, so the Fixman-compensated pdf will not agree with the unconstrained pdf: the Fixman potential alone cannot overcome the statistical biases for such forcefields.

#### 2. Linear C4 chain with separable degrees of freedom

Since $\det\{\mathcal{M}_B\}$ (eq:7) is torsion-independent, when $\alpha$ are just torsions (TMD) eqs 11 and 12 simplify to eq:13, with $c_f$ a constant from the bond-length/bond-angle contributions of $\det\{\mathcal{M}_B\}$. For separable potentials the Fixman potential should exactly compensate the config-pdf difference.

In prior work with $U(\alpha)=0$, the torsion pdf is uniform, eq:14, and Fixman correctly compensates the bias for serial and branched systems.

##### a. Calculation of transition barrier crossing rates

The Fixman potential cannot fully remove biases when the quantity of interest depends on velocity coordinates. Studying a C4 system with a single-barrier torsional potential $U(\alpha)$ (separable from bond length/angle), Fixman corrects the torsion pdf but only partially the barrier-crossing rates - extending Pear and Weiner's results.

The single-barrier harmonic torsional potential is eq:15 with $k_\alpha=0.30$ kcal/mol; $\alpha_0$ is the peak location. Transition-state rate theory assumes that once the torsion crosses the barrier peak it does not return (a reasonable approximation for Langevin dynamics at low viscosity). The generic 1D transition-state rate is eq:16; for the C4 torsion it is eq:17. For the constrained model (see Appendix) it simplifies to eq:18 with $S^{-1}(\alpha_0)$ given by eq:18b. The rate depends on the mass-matrix determinant. When the Fixman potential is applied, $\det\{\mathcal{M}\}$ drops out (eq:19), but the $S^{-1}(\alpha_0)$ term remains - a residual dependency on the barrier-peak location that the unconstrained model does not exhibit. Thus Fixman fully compensates config-variable pdfs but not velocity-dependent quantities.

##### b. Calculation of barrier crossing rates from simulations

Three simulation types: FLEXIBLE (all-atom Cartesian, no constraints); TMD (ICMD with bonds and angles frozen); FIXMAN (TMD plus Fixman potential). For the C4 system: bond angles 90°, bond lengths 1.54 Å, masses 14.01 amu; barrier locations 0°, 45°, 90°, 135°, 180°; three 20 ns Langevin runs, time step 1 fs, damping 0.1/fs. FLEXIBLE bond/angle spring constants 303.1 kcal/Å² and 63.21 kcal.

Results (Fig 2): for $\alpha_0=90°$, $T=800$ K, the TMD torsion pdf is biased (shifting the effective barrier peak); only Fixman recovers the expected pdf. The TMD pdf RMS deviation from FLEXIBLE also depends on the barrier location; Fixman removes both the bias and this dependence. However, for barrier centers at $\alpha_0=90°$ and $0°$, Fixman does not sufficiently compensate the bias in the barrier crossing rates, which retain a barrier-location dependence (independent of temperature). Fixman does reduce the error and allows correctly identifying the location and magnitude of the potential barrier.

#### 3. Linear chain C4 with non-separable degrees of freedom

When torsions are cross-coupled to bond lengths/angles by a non-separable potential $U(\alpha,q)$, the unconstrained pdf (eq:6) equals the corrected constrained pdf (eq:13) only when $U(\alpha,q)$ is a very steep well around $q_0$. Otherwise $q$ varies significantly and the Fixman-compensated pdf differs from the unconstrained pdf - the difference is the intrinsic + extrinsic distortion.

A Coulombic coupling potential (eq:20) is added to C4 terminal atoms ($k_{\text{coul}}=332.06$ kcal·Å/e²); $r$ depends on both torsion and bond angles, making the potential non-separable. Opposite-sign charges make $U_{\text{coul}}$ attractive with $\rho(\alpha)$ maximal at $\alpha=0°$. Harmonic angle springs (eq:21) with constant $k_\theta$ control stiffness.

Applying charges $+0.2e$ and $-0.2e$, bond lengths 1.54 Å, angles 90°, masses 14.01 amu, and varying $k_\theta$: the Fixman potential removes the torsion-pdf bias only for very stiff $k_\theta$. For realistic AMBER99SB values (30-100 kcal), the FLEXIBLE and FIXMAN torsion pdfs differ (Fig 3, $k_\theta=95$ kcal), agreeing only at large $k_\theta$. So the torsion-angle cross-coupling is not negligible: Fixman compensates the intrinsic distortion but not the extrinsic distortion.

### B. Hybrid ICMD to correct for extrinsic distortions

All-atom forcefields take the form eq:22. It is customary to assume the harmonic bond/angle restraining terms dominate, so the torsional-subspace potential is independent of bond/angle values; then Fixman alone suffices. But Echenique et al. showed that for realistic forcefields there is significant torsion-bond-angle coupling and the harmonic terms do not dominate. Then the forcefield is non-separable and torsion statistics depend on bond-angle DOF too; constraining bond angles introduces both intrinsic and extrinsic distortion, and Fixman does not eliminate the latter.

The proposed **hybrid internal coordinate molecular dynamics** method treats any desired subset of bond angles as flexible at any point during the simulation, keeping bonds and remaining bond angles rigid - bridging all-atom flexible and torsional MD. Implemented in *GneimoSim*. Motivation: Berkholz et al. showed protein backbone covalent geometry is a function of backbone torsions; Hinsen et al. showed including a subset of backbone angles in the unconstrained set accurately represents protein structure. Postulate: opening only specific backbone bond angles compensates the extrinsic distortion and retrieves unconstrained-like landscapes.

## III. Simulation results for the hybrid ICMD method

The hybrid ICMD method is applied to small peptides with all-atom forcefields. For non-separable potentials, constraints distort the free energy surface (FES); long-range forces couple torsions to bond lengths/angles, a distortion Fixman alone cannot compensate. With hybrid ICMD + Fixman, opening a few key bond angles removes these distortions. Even with these angles open, stable simulations are possible with time steps up to 5 fs.

### A. Application of hybrid ICMD to alanine dipeptide

TMD, FIXMAN, and FLEXIBLE simulations. Alanine dipeptide is modeled as rigid-body clusters connected by hinges: each cluster is a non-terminal atom with all its terminal neighbors; aromatic rings treated as rigid clusters; proline rings and disulfide bonds broken into tree structures held by stiff harmonic bond parameters. Software: *GneimoSim* with AMBER99SB and GBSA solvation.

Conformations described by backbone dihedrals $\phi$ = C–N–Cα–C and $\psi$ = N–Cα–C–N; FES projected onto these. 20 FLEXIBLE and 20 TMD simulations each at 300 K and 800 K, each 20 ns. FIXMAN simulations at 300 K then systematically opened bond angles, 20 simulations of 20 ns each per angle/combination. All constrained simulations start from an initial alpha-helical conformation $(\phi,\psi)=(-60°,-40°)$, Langevin dynamics, damping 0.1/fs, time step 2 fs.

#### 1. Results and discussion

At 300 K, TMD leads to barriers limiting sampling of the minima in the 1st and 4th FES quadrants; partially overcome at 800 K but still limited vs FLEXIBLE. Adding Fixman still fails to sample the 1st/4th-quadrant minima - the torsion-angle cross-coupling significantly affects the FES (extrinsic bias uncompensated by Fixman). Opening backbone angles **C–N–Cα and N–Cα–C together** recovers the FLEXIBLE features; opening additional angles has little effect, and opening only one of the two does not alleviate the barriers.

### B. Application of hybrid ICMD to other dipeptides

Applied to valine, leucine, isoleucine, methionine, phenylalanine, tryptophan, proline, tyrosine dipeptides. Separate simulations with each backbone bond angle open, and combinations; 20 simulations of 20 ns each, clustering/methods as in Section III A. Start from alpha-helical $(\phi,\psi)=(-60°,-40°)$, Fixman enabled.

#### 1. Results and discussion

Applying Fixman + hybrid ICMD with open angles C–N–Cα and N–Cα–C lets VAL, LEU, MET, PHE, TRP, TYR sample the first quadrant. Quantified via Hellinger distance between FIXMAN/FLEXIBLE and hybrid/FLEXIBLE pdfs. Opening additional backbone angles does not significantly improve PHE, ALA, LEU, TRP. For ILE, VAL, MET, TYR, additionally opening the Cα–C–N backbone angle gives maximal improvement (Hellinger distances 0.14, 0.13, 0.17, 0.16 for VAL, MET, TYR, ILE respectively). Proline (a ring/loop) is broken to a tree; opening the backbone N–Cα–C plus sidechain Cα–Cβ–Cγ gives a minimal Hellinger distance of 0.26 (vs 0.45 for FIXMAN). Opening angles other than the specified ones had little effect.

Opening only N–Cα–C (as suggested by Arnautova et al. and Hinsen et al.) produced pdfs quite distant from FLEXIBLE; in all cases except proline, at least the additional angle C–N–Cα was needed to significantly reduce the extrinsic bias.

### C. Application of hybrid ICMD to longer peptide chains

Applied to the ten-amino-acid peptide **CLN025** (NMR ensemble available, PDB ID 2RVD). Eight simulations of 10 ns each for FIXMAN, hybrid ICMD, FLEXIBLE. Hybrid ICMD opened N–Cα–C for proline and C–N–Cα–C and N–Cα–C backbone angles for other residues. 300 K via Nose-Hoover, time step 1 fs. GBSA solvation, internal dielectric 4.0, external 78.0. All start from the 2RVD NMR structure.

#### 1. Results and discussion

Percentage of snapshots with hydrogen-bond distances within one standard deviation of the NMR mean (Table I). Only a very small percentage of FIXMAN conformations agree with observed NMR distances; opening backbone bond angles in hybrid ICMD produces a much larger proportion matching the NMR ensemble. Per-residue backbone Ramachandran $(\phi,\psi)$ pdfs show hybrid ICMD gives greater conformational sampling than FIXMAN, with Hellinger distances (to FLEXIBLE) much lower than FIXMAN's. Fixing all bond angles introduces an extrinsic distortion limiting sampling; opening only backbone bond angles significantly ameliorates it.

### D. Time step size for the hybrid ICMD method

Opening angle DOF can require smaller time steps. Hybrid ICMD simulations of dipeptides (N–Cα–C and C–N–Cα open for ALA, PHE, LEU, TRP; additional Cα–C–N for ILE, VAL, MET, TYR; N–Cα–C + Cα–Cβ–Cγ for PRO) used time steps of 5 fs, 20 simulations of 50 ns each. Opening backbone bond angles in coarse-grain ICMD still allows stable use of 5 fs time steps (much larger than all-atom Cartesian MD).

## IV. Conclusions and future work

Distortions in the constrained-ICMD potential energy surface stem from (i) the mass matrix's dependence on BAT coordinates (intrinsic) and (ii) forcefield coupling between constrained and unconstrained DOF (extrinsic). The Fixman potential alleviates the intrinsic distortion in the configuration pdf but only partially overcomes biases for velocity-dependent quantities (e.g. barrier-crossing rate). Long-range nonbonded forces couple constrained and unconstrained DOF, introducing an extrinsic distortion Fixman cannot compensate. The hybrid ICMD method lets the user open any bond-angle DOF; for short peptides, Fixman plus opening a small subset of backbone angles removes the extrinsic distortions while permitting large time steps, and recovers the FES of a complex polypeptide. The hybrid ICMD method (in GneimoSim) now supports dynamics models at multiple coarsening levels - from all angles open to arbitrarily large rigid domains - in the same simulation (e.g. rigid helices with flexible loops). Future work: validating hybrid ICMD for larger proteins in explicit solvent, extending to nucleic acids, and extending the Fixman potential to correct velocity-dependent biases such as barrier-crossing rates.

## Appendix: Calculating the barrier-crossing rate

### Derivation (not implemented)

For a C4 system, the transition-state barrier-crossing rate is eq:17 (repeated as A1). For the constrained model the pdf is eq:A2, with kinetic energy eq:A3; the $\det\{\mathcal{M}(\alpha)\}$ prefactor comes from the momentum->velocity coordinate change.

For $\alpha,\gamma$ of dimensions $m,n$, partition the mass matrix as eq:A4, with Schur complements $S=S_0-VW_0^{-1}V^*$ and $W=W_0-V^*S_0^{-1}V$, giving

$$ \mathcal{M}^{-1}(\alpha) = \begin{bmatrix} S^{-1} & -S_0^{-1}VW^{-1} \\ -W_0^{-1}V^{*}S^{-1} & W^{-1} \end{bmatrix} $$

and the determinant factorization eq:A5. Note $S^{-1}(\alpha)=[\mathcal{M}^{-1}(\alpha)]_\alpha$ (eq:A6). The kinetic energy re-expressed:

$$ E_k(\alpha,\dot\alpha,\dot\gamma) = \tfrac12 \dot\alpha^{*} S \dot\alpha + \tfrac12 (\dot\gamma-\beta)^{*} W_0 (\dot\gamma-\beta), \qquad \beta = -W_0^{-1}V^{*}\dot\alpha $$

For C4, $\alpha$ is the scalar torsion and $\gamma$ the three Euler angles, so $S$ is scalar and $\det\{S\}=S$. Substituting eq:A2 and the above into A1 with $\alpha=\alpha_0$ gives a product of Gaussian integrals (eq:A8). Using the standard Gaussian integrals eq:A9a (with $A=W_0(\alpha_0)/kT$, $p=3$) and eq:A9b (with $s=S(\alpha_0)/kT$) and eq:A5, one obtains eq:18. When the Fixman potential (eq:13) is included in $\mathcal{U}$, the $\det\{\mathcal{M}(\alpha_0)\}$ term drops out (eq:19); the residual $S(\alpha_0)$ dependence on barrier location remains. Hence Fixman fully corrects config-dependent averages but only partially corrects velocity-dependent ones.
