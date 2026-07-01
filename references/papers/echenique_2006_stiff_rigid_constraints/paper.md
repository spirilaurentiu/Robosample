# Quantum Mechanical Calculation of the Effects of Stiff and Rigid Constraints in the Conformational Equilibrium of the Alanine Dipeptide

Pablo Echenique, Iván Calvo, J. L. Alonso. *J Comput Chem* 27: 1733–1747, 2006. DOI 10.1002/jcc.20467.

## Abstract
If constraints are imposed on a macromolecule, two inequivalent classical models may be used: the **stiff** and the **rigid** one. This work studies the effects of such constraints on the conformational equilibrium distribution (CED) of the model dipeptide HCO-L-Ala-NH2 without any simplifying assumption. Ab initio quantum mechanics calculations including electron correlation at the MP2 level describe the system, and the conformational dependence of all correcting terms to the naive CED (based on the potential energy surface) that appear when constraints are considered is measured. These terms are related to mass-metric tensor determinants and also occur in Fixman's compensating potential. Some corrections are non-negligible over the whole Ramachandran space; if only the energetically lower region (containing the principal secondary structure elements) is relevant, all correcting terms may be neglected up to peptides of considerable length.

Key words: constraints; alanine dipeptide; Fixman; ab initio; mass metric tensors.

## Introduction (routing summary)
Reducing DOF of macromolecules via constraints is standard. Two classical models exist and are inequivalent: the **classical rigid** model (constraints exact; velocities orthogonal to the constraint hypersurface vanish) and the **classical stiff** model (constraints approximate; a steep potential drives the system to the hypersurface, orthogonal velocities activated as "heat containers"). The paper does not decide which model is physically correct; it studies both on equal footing.

In internal coordinates the constraints are directly imposed on selected hard coordinates (bond lengths, bond angles, some dihedrals). The marginal probability density in the coordinate part of phase space is NOT proportional to the naive `exp[-β V_Σ(q^i)]`; correcting terms involving determinants of mass-metric tensors and of the Hessian of the constraining potential must be added. These same terms define Fixman's compensating potential, used to reproduce the stiff equilibrium distribution from rigid MD.

## Theory

### Setup and notation
The system is n mass points (atoms). Euclidean coordinates `x⃗_α` (α=1..n). Curvilinear coords `q^μ` (μ=1..3n), with `N := 3n`. Split `q^μ = (q^A, q^a)`: first six `q^A` are **external** (overall translation+rotation), the `q^a` (a=7..N) are **internal**. Internals split into soft `q^i` (i=7..M+6) and hard `q^I` (I=M+7..N). All soft coords `q^u = (q^A, q^i)`, u=1..M+6. See `notation.md` (Table 1) for index ranges.

The hypersurface Σ ⊂ I of dimension M is described by L = N-M-6 constraints (eq 1): `q^I = f^I(q^i)`. We seek the probability density on E × Σ rather than the whole space E × I.

### Classical Stiff Model
Split the potential into value on Σ plus a constraining potential (eq 2). Impose on V_c: (i) Σ is the global minimum of V_c wrt hard-coordinate variations; (ii) small hard-coordinate variations Δq^I produce energy changes ≫ thermal energy RT.

Advantages of this formulation: provides a direct prescription to compute V_Σ(q^i) and Σ (the Potential Energy Surface, PES) via geometry optimization at fixed soft coordinates; clearly shows the necessity of the Hessian-determinant correction term.

The soft coordinates `q^i` need NOT be soft energetically (they may be voluntarily chosen, e.g. Ramachandran φ,ψ whose barriers can be ~40 RT); only the hard coordinates must truly be hard. The labels "important"/"unimportant" (Karplus-Kushick) may be more accurate than soft/hard.

### Derivation (not implemented) - stiff partition function reduction
Because of condition (ii), only the vicinity of the equilibrium hard-coordinate values matters, so V_c is Taylor-expanded to second order about Σ (eq 3). The zeroth-order term is zero by definition of V_c, the linear term is zero by condition (i); the first non-zero term is quadratic, defining the partial Hessian H_JK (only hard-coord derivatives). This yields the stiff Hamiltonian (eq 4) with mass-metric tensor G (eq 5) and its inverse (eq 6), and the stiff partition function (eq 7).

Using condition (ii), the hard coords in G are evaluated at equilibrium f^I(q^i) (eq 8). Integrating over the hard coords gives a Gaussian integral producing `det^{-1/2} H`, taken to the exponent (eq 9). Integrating over the momenta produces `det^{-1/2} G` to the exponent, yielding the coordinate-only stiff partition function (eq 10) with prefactor χ_s(T) (eq 11). No Jacobian appears in the phase-space measure because q^μ, p_μ come from Euclidean coordinates via a canonical transformation. The partial Hessian is positive definite so det H > 0 and its log is well defined.

The exponent is read as a free energy: V_Σ is the internal energy, and two conformation-dependent correcting terms are effective entropies (linear in RT). S_s^k (from det G) is a "kinetic entropy" from averaging momenta; S_s^c (from det H) is a true configurational entropy from averaging out DOF. Definitions: F_s (eq 12a), S_s^c (eq 12b), S_s^k (eq 12c). The stiff equilibrium probability is eq 13.

Although S_s^k depends on external coords q^A, det G factorizes (ref 58) into an external-only times an internal-only function, so the external factor integrates out independently (see Factorization section).

### Classical Rigid Model
Treating eq 1 as exact holonomic constraints, the reduced Hamiltonian on E×Σ is eq 14 with the **reduced mass-metric tensor** g (eq 15), the pull-back of the full G through the constraint map f̃ (eq 16).

### Derivation (not implemented) - rigid reduction
Eq 15 follows from the unconstrained Hamiltonian (eq 17) using the constraints (eq 1) and their time derivatives (eq 18), defining reduced momenta η_ν (eq 19). The rigid partition function (eq 20), integrated over momenta, gives the coordinate-only form (eq 21) with prefactor χ_r(T) (eq 22). Free energy F_r (eq 23a), kinetic entropy S_r^k (eq 23b), rigid equilibrium probability (eq 24).

As with G, det g factorizes (ref 58) into external-only × internal-only functions; the external part integrates out.

### Fixman's compensating potential
Defined as F_s - F_r (eq 25):
`V_F(q^u) = (RT/2) ln[ det G / (det H · det g) ]`.
Performing rigid MD (equilibrium ∝ exp[-β F_r]) and adding V_F to V_Σ reproduces the stiff density P_s ∝ exp[-β F_s]. Historically this application (not MC sampling) motivated interest in mass-metric tensor effects.

Table 2 summarizes both models' equilibrium densities and correcting terms.

## Methods

### Factorization of the External Coordinates
Using SASMIC coordinates (ref 59), det G has the closed form eq 26: a product of atomic-mass cubes, one external `sin^2 θ`, bond-length^4 factors, and bond-angle `sin^2 θ_α` factors - independent of dihedral angles explicitly (Go-Scheraga / Volkenstein result). The mass factor is conformation-independent (droppable); the only external-dependent part `sin^2 θ` integrates out. The stiff kinetic entropy reduces to eq 27 (sums of ln r_α^4 and ln sin^2 θ_α).

Similarly det g factorizes as `sin^2 θ · det g_2(q^i)` (eq 28), with g_2 the internal reduced mass-metric matrix (eq 29) built from total mass m_tot, the skew matrix v(R⃗) (eq 31), the inertia tensor J (eq 30), and derivatives of atom positions in the primed (body-fixed) frame wrt soft internals. Here `m_tot = Σ_α m_α`, `R⃗ = m_tot^{-1} Σ_α m_α x'_α`. After integrating out sin^2 θ, the rigid kinetic entropy depends only on soft internals (eq 32). Because sin^2 θ divides out in eq 25, V_F is independent of the external coordinates too.

### Computational Methods
For HCO-L-Ala-NH2 (Fig 1): M=2 soft internals q^i = (φ,ψ) (Ramachandran angles, Table 3), N=48, L=40, n=16 atoms. Side-chain χ and peptide-bond ω_0, ω_1 dihedrals are treated as hard (barriers ~6-12 RT for χ, larger for ω); the χ integral is approximated by three Gaussian integrals (threefold symmetry), adding only a T- and conformation-independent constant to S_s^c.

Ab initio calculations with GAMESS (ref 79). SASMIC coordinates (ref 59), converted to Delocalized Coordinates (ref 80) for optimization. PES computed on a 12×12 grid of (φ,ψ), each angle from -165° to 165° in 30° steps, via constrained MP2/6-31++G(d,p) optimizations freezing φ,ψ; OPTTOL=1e-5, CONV=1e-6. Result: 144 conformations defining Σ and V_Σ(φ,ψ) (~100 days CPU). Hessian H(φ,ψ) computed at each grid point removing the φ,ψ rows/columns (~140 days CPU).

Eqs 27, 32 give the kinetic entropy terms. For g_2 (eq 29) the primed-frame Euclidean coords x'_α of the 16 atoms and their derivatives wrt (φ,ψ) are needed; two additional 12×12 grids displaced +2° in φ and +2° in ψ give the derivatives by finite differences (~75 days CPU each). All correcting terms also computed at six special secondary-structure points (~16 days). Total ~406 days CPU at MP2. HF/6-31++G(d,p) repeated at ~1/10 the cost (~40 days total).

## Results
See `checks.md` for all numeric tables (Tables 4-8).

- Table 4: max variation / average / std of V_Σ, F_s, F_r and the three correcting terms plus V_F over the grid.
- The conformational dependence of the correcting terms is >1 order of magnitude smaller than that of V_Σ in the worst case, but comparable to chemical accuracy (1 kcal/mol) and to the MP2-vs-HF difference; so they can be relevant.
- Correction-term magnitude ordering: the least important is S_s^k (mass-metric G); the most important is S_s^c (Hessian H) - the latter is persistently underestimated in the literature.

The paper introduces a statistical **distance** d_12 (ref 82) between two potentials: the typical error in energy differences when substituting V_2 for V_1 (linear rescaling allowed). d_12 < RT ⇒ safe substitution. The residue limit for an additive per-residue polypeptide potential is N_res = (RT/d_12)^2 (eq 33). Full-grid results (Table 6, RT≈0.6 kcal/mol at 300 K):
- Dropping both stiff terms (F_s vs V_Σ): d_12=0.74 RT, N_res=1.82 → error above thermal noise already at 2 residues.
- Dropping only S_s^k: d_12=0.11 RT, N_res≈80 → safe up to ~80 residues.
- Dropping S_r^k (F_r vs V_Σ): d_12=0.29 RT, N_res≈12.
- Omitting V_F (F_r vs F_s): d_12=0.67 RT, N_res=2.24 → include Fixman for peptides > 2 residues.

Restricting to the six secondary-structure elements (Table 7, low-energy region) the distances shrink (~4× larger residue limits, Table 8), because σ_2 ≈ 2 kcal/mol there vs ≈ 4 kcal/mol over the whole grid, and d_12 = √2 σ_2 (1-r_12^2)^{1/2} (eq 34).

HF results closely match MP2 (correlations ~0.91-0.98), differing mostly by a linear scale (b_12 > 1) - so future studies can use HF at 1/10 the cost.

## Conclusions
- Stiff MC at room temperature: -TS_s^k (mass-metric G) negligible up to ~80 residues (max variation 0.24 kcal/mol).
- Stiff MC: -TS_s^c (Hessian H) should be included for peptides > 2 residues (max variation 1.67 kcal/mol).
- Rigid MC: -TS_r^k (reduced g) negligible up to ~12 residues (max variation 0.81 kcal/mol).
- Rigid MD aiming at the stiff distribution: Fixman V_F should be included for peptides > 2 residues (max variation 1.68 kcal/mol).
- Restricting to the low-energy secondary-structure region multiplies these residue limits by ≈ 4.
- Errors from dropping the most important corrections are the same order as errors from MP2→HF.
- Caveat: conclusions are for a QM-derived potential on a small dipeptide; simpler force fields and bulkier residues should be re-examined.

## Appendix (Derivation / discussion - not implemented)
Three common simplifying approximations: (i) neglect conformational dependence of det G; (ii) neglect conformational dependence of det H; (iii) assume hard coords are constant (f^I independent of q^i).

Under approximations (i)+(ii) the Fixman potential depends only on det g. Under approximation (iii), g is the soft sub-block of G and Fixman showed det G / det g = 1/det h, with h^{IJ} the hard-coord matrix (eq A1) - sparse because each internal coord involves few atoms.

For realistic force fields (eq A2: harmonic bonds+angles + torsions + long-range), det H is conformation-dependent (long-range terms affect the Hessian on Σ) and the equilibrium hard coords are functions f^I(q^i), NOT the constants r_α^0, θ_α^0. So approximations (i)-(iii) must be checked case by case. One may define "exactly separable" coordinates Q^a (eq A3) with hard coords constant on Σ, but their transformation to Euclidean requires knowing f^I numerically, limiting practical use.
