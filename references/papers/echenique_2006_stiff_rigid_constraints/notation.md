# Notation - Echenique 2006 (stiff vs rigid constraints)

## Conventions (READ FIRST)
- **Energy units:** kcal/mol throughout results tables; per-mole energy units, so `β := 1/RT` with R the gas constant (NOT `1/k_B T`). RT ≈ 0.6 kcal/mol at 300 K.
- **Einstein summation** on repeated upper/lower indices.
- **Primed frame** `x'_α` = body-fixed reference frame (translation/rotation removed).
- Mass-metric tensor lowered-index `G_{μν}` is the matrix `G`; raised-index `G^{μν}` is its inverse (eq 6). Same for `g`.
- All correcting terms are "referenced to zero in the grid" in the results tables (only conformational variation matters).

## Index ranges (Table 1)
| Index symbol | Range | Count | Meaning |
|---|---|---|---|
| α, β, γ | 1..n | n | atoms |
| μ, ν, ρ | 1..N | N=3n | all coordinates |
| A, B, C | 1..6 | 6 | external coordinates (overall translation+rotation) |
| a, b, c | 7..N | N-6 | internal coordinates |
| i, j, k | 7..M+6 | M | soft internal coordinates |
| I, J, K | M+7..N | L=N-M-6 | hard internal coordinates |
| u, v, w | 1..M+6 | M+6 | all soft coordinates (external + soft internal) |

## Symbols
| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| n | number of atoms | int | dipeptide: 16 |
| N | total DOF = 3n | int | dipeptide: 48 |
| M | number of soft internal coords | int | dipeptide: 2 (φ,ψ) |
| L | number of hard internal coords = N-M-6 | int | dipeptide: 40 |
| q^μ | generalized (curvilinear) coordinate | R | split into (q^A external, q^a internal) |
| q^A | external coordinates | R^6 | overall position + orientation |
| q^i | soft internal coordinates | R^M | dipeptide: (φ,ψ) Ramachandran angles |
| q^I | hard internal coordinates | R^L | constrained: q^I=f^I(q^i) |
| q^u | all soft coords (q^A, q^i) | R^{M+6} | |
| x^σ, x⃗_α | Euclidean (Cartesian) coordinates of atoms | R, R^3 | space-fixed frame |
| x'_α | atom position in body-fixed (primed) frame | R^3 | translation/rotation removed |
| m_α, m_σ | atomic mass | mass | m_σ = mass of atom owning Euclidean coord σ |
| m_tot | total mass Σ m_α | mass | |
| R⃗ | center of mass in primed frame | R^3 | R⃗ = m_tot^{-1} Σ m_α x'_α |
| f^I(q^i) | constraint function (hard from soft) | R | eq 1 |
| f̃^μ | full parametrization map | eq 16 | identity on soft, f^I on hard |
| V(q^a) | full potential energy | kcal/mol | Born-Oppenheimer PES (ab initio MP2/HF) |
| V_Σ(q^i) | potential on constraint surface Σ (the PES) | kcal/mol | from constrained geom optimization |
| V_c | constraining potential | kcal/mol | zero on Σ (eq 2) |
| G_{μν} | mass-metric tensor | mass, N×N | eq 5 |
| g_{vw} | reduced mass-metric tensor on E×Σ | mass, (M+6)×(M+6) | eq 15 |
| g_2 | internal reduced mass-metric matrix | (M+3)×(M+3) block | eq 29 (det g = sin^2θ det g_2) |
| H_{JK} | partial Hessian of V_c over hard coords on Σ | kcal/mol, L×L | positive definite -> det H>0 |
| J | inertia tensor in primed frame | mass·length^2, 3×3 | eq 30 |
| v(R⃗) | skew cross-product matrix of R⃗ | length, 3×3 | v(R)w = R×w (eq 31) |
| p_μ | momentum conjugate to q^μ | | full space |
| η_u | reduced momentum conjugate to soft coord | | rigid space (eq 19) |
| β | inverse thermal energy = 1/RT | mol/kcal | per-mole units |
| R | gas constant | kcal/(mol·K) | |
| T | temperature | K | 300 K used |
| h | Planck constant | | |
| α_QM | quantum indistinguishability factor | dimensionless | 1/N! for N identical particles |
| H_s, H_r | stiff / rigid Hamiltonian | kcal/mol | eqs 4, 14 |
| Z_s, Z_r | stiff / rigid partition function | | |
| Z'_s, Z'_r | configurational normalizers | | eqs 13, 24 |
| χ_s(T), χ_r(T) | temperature prefactors | eqs 11, 22 | |
| F_s, F_r | stiff / rigid effective free energy | kcal/mol | eqs 12a, 23a - the sampling potentials |
| S_s^c | stiff conformational entropy (Hessian) | kcal/(mol·K) | -R/2 ln det H (eq 12b) |
| S_s^k | stiff kinetic entropy (mass-metric G) | kcal/(mol·K) | +R/2 ln det G (eq 12c) |
| S_r^k | rigid kinetic entropy (reduced g) | kcal/(mol·K) | +R/2 ln det g (eq 23b) |
| V_F | Fixman compensating potential | kcal/mol | F_s - F_r (eq 25) |
| P_s, P_r | equilibrium probability densities | | eqs 13, 24 |
| h^{IJ} | Fixman sparse hard-coord matrix | eq A1 | det G/det g = 1/det h under approx (iii) |
| r_α | bond length of atom α | length (Angstrom) | SASMIC Z-matrix |
| θ_α | bond angle of atom α | radians | SASMIC Z-matrix |
| φ, ψ | Ramachandran backbone dihedrals | degrees | soft coords of dipeptide |
| ω_0, ω_1 | peptide-bond dihedrals | degrees | hard (large barriers) |
| χ | side-chain dihedral | degrees | hard (barriers ~6-12 RT) |
| d_12 | statistical energy distance | RT units | eq 34 |
| N_res | residue limit for dropping a term | dimensionless | eq 33 |
| b_12 | slope of linear rescaling V_1 vs V_2 | dimensionless | |
| r_12 | Pearson correlation coefficient | dimensionless | |
| σ_2 | std deviation of V_2 over working set | kcal/mol | |
