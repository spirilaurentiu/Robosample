# Equations - ICFF (Katritch, Totrov, Abagyan 2003)

<!-- eq:1 -->
$$ E_{\text{tor}} = C_0 + \sum_{k=0}^{k=6} (A_k \cos k\theta + B_k \sin k\theta) $$
- **what:** ICFF composite torsion energy for a bond, a sixfold (0..6) Fourier series in the torsion angle. This term implicitly absorbs bond stretching, bond bending, out-of-plane, and "1-4" van der Waals contributions of the source Cartesian force field.
- **symbols:** $E_{\text{tor}}$ - torsion energy (kcal/mol); $\theta$ - torsion angle about the specific bond (rad or deg, argument of cos/sin so radians for evaluation); $C_0$ - constant offset (kcal/mol); $A_k, B_k$ - Fourier cosine/sine coefficients for harmonic $k$ (kcal/mol); $k$ - harmonic index, integer $0 \le k \le 6$.
- **note:** $C_0$ and the $k=0$ term ($A_0\cos0 + B_0\sin0 = A_0$) are both constants; the paper writes both explicitly. In practice the constant part is $C_0 + A_0$. Coefficients are fit from a Cartesian energy profile of a "torsion fragment" (see checks.md for grid).

<!-- eq:2 -->
$$ E_{\text{vw}(1\text{-}5,1\text{-}6)} = \varepsilon_{IJ} \left[ C \left( \left( \frac{R_{IJ} - R_{IJ}^*}{R_{IJ}^0 - R_{IJ}^*} \right)^2 - 1 \right) + (1 - C) \left( \left( \frac{R_{IJ} - R_{IJ}^*}{R_{IJ}^0 - R_{IJ}^*} \right)^3 - 1 \right) \right], \quad \text{if } R_{IJ} < R_{IJ}^* $$
- **what:** Soft empirical repulsion term for "1-5" and "1-6" atom pairs (separated by four or five covalent bonds), applied only in the repulsive region $R_{IJ} < R_{IJ}^*$. A blend of a harmonic (square) and cubic polynomial in the reduced distance, mimicking covalent-bond flexibility so rigid-geometry clashes are softened.
- **symbols:** $E_{\text{vw}(1\text{-}5,1\text{-}6)}$ - soft repulsion energy (kcal/mol); $\varepsilon_{IJ}$ - well depth (minimal energy magnitude) of the source vdW surface for atom types $I,J$ (kcal/mol, taken as the vdW minimum value, so the term equals $\varepsilon_{IJ}$ at the minimum); $R_{IJ}$ - actual distance between atoms $I$ and $J$ (Å); $R_{IJ}^*$ - distance at the vdW energy minimum (Å); $R_{IJ}^0$ - distance where the source vdW surface crosses zero (Å); $C$ - blend factor, $0 \le C \le 1$, the only adjustable parameter (chosen $C=0.55$).
- **boundary conditions (from prose, use to verify implementation):** at $R_{IJ}=R_{IJ}^*$, $E = \varepsilon_{IJ}$ with zero derivative; at $R_{IJ}=R_{IJ}^0$, $E = 0$ (matching the original vdW zero crossing). For $R_{IJ} \ge R_{IJ}^*$ the original attractive vdW branch is kept unchanged. The term is finite at $R_{IJ}=0$ (requires a capped electrostatic term to avoid collapse).
<!-- CHECK: paper states E = eps_IJ at the minimum; here eps_IJ denotes the (magnitude of the) minimum energy of the vdW surface. Verify sign convention against the source MMFF94 vdW well depth when porting. -->

<!-- eq:3 -->
$$ E_{\text{MMFF}} = \sum EB_{ij} + \sum EA_{ijk} + \sum EBA_{ijk} + \sum EOOP_{ijk;l} + \sum ET_{ijkl} + \sum EvdW_{ij} + \sum EQ_{ij} $$
- **what:** Full MMFF94(s) Cartesian potential energy used as the "source" force field. Sum of bond stretching, bond bending, stretch-bend, out-of-plane bending, torsion, van der Waals, and electrostatic terms. Full term definitions are in eqs. (1)-(14) of Halgren 1996 (ref. 28); this paper only reproduces the top-level decomposition.
- **symbols:** $E_{\text{MMFF}}$ - total MMFF94 energy (kcal/mol); $EB_{ij}$ - bond stretch (bond $i$-$j$); $EA_{ijk}$ - bond angle bend; $EBA_{ijk}$ - stretch-bend coupling; $EOOP_{ijk;l}$ - out-of-plane bend at center $l$; $ET_{ijkl}$ - torsion; $EvdW_{ij}$ - van der Waals; $EQ_{ij}$ - electrostatic (Coulomb).

<!-- eq:constraint -->
$$ E_r = C_r (\theta - \theta^0)^2, \qquad C_r = 10000 \text{ kcal} $$
- **what:** Harmonic restraint potential used to constrain a torsion angle to a target value $\theta^0$ during the grid-scan Cartesian minimizations that generate the torsion energy profile. With this stiffness the angle deviation from $\theta^0$ typically stays below 0.1 degree.
- **symbols:** $E_r$ - restraint energy (kcal/mol); $\theta$ - current torsion angle; $\theta^0$ - target constrained angle; $C_r$ - restraint force constant, $10000$ kcal/mol (per rad^2 or deg^2 as used).

<!-- eq:weight -->
$$ W(\theta) = \left[ E(\theta) - \min(E(\theta)) + 1 \right]^{-1} $$
- **what:** Per-datapoint weight for the weighted least-squares fit of the Fourier coefficients $\{A_k, B_k\}$ (eq. 1). Weight is inversely proportional to the relative energy, emphasizing the low-energy part of the torsion profile. Improves fit RMSD to ~0.001 kcal below the 5-kcal cutoff and down-weights inaccurate high-energy regions.
- **symbols:** $W(\theta)$ - fit weight at torsion angle $\theta$ (dimensionless); $E(\theta)$ - source (local) energy at $\theta$ (kcal/mol); $\min(E(\theta))$ - minimum energy over the profile (kcal/mol); the "$+1$" is a 1-kcal regularizer keeping the weight finite at the minimum.
<!-- CHECK: OCR gave "W() [E() Min(E()) 1]^1"; reconstructed as [E - min(E) + 1]^-1 from the prose ("weight ... proportional to the inverse relative energy"). -->
