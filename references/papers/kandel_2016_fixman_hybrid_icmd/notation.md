# Notation - Kandel et al. 2016

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| $\alpha$ | unconstrained generalized coords (torsions + open bond angles) | rad, $R^N$ | ICMD internal coords |
| $q$ | constrained BAT coords (frozen bonds/angles) | mixed, $R^{3n-N}$ | frozen at $q_0$ in constrained model |
| $q_0$ | frozen values of constrained coords | - | reference geometry |
| $n$ | number of atoms | int | |
| $N$ | number of unconstrained (movable) DOF | int | $\le 3n$ |
| $\mathfrak{p}$, $p$ | canonical momenta | $R^{3n}$ | conjugate to $(\alpha,q)$ |
| $\mathcal{M}(\alpha)$ | constrained-model mass matrix | $R^{N\times N}$, symmetric PD | reduced (articulated) inertia |
| $\mathcal{M}_B(\alpha,q)$ | full BAT mass matrix (unconstrained) | $R^{3n\times 3n}$ | see eq:7 for $\det$ |
| $\det\{\mathcal{M}^{1/2}\}$ | shorthand for $(\det\{\mathcal{M}\})^{1/2}$ | scalar | |
| $\mathcal{U}, U$ | forcefield potential energy | kcal/mol | |
| $\mathcal{U}_f, U_f$ | Fixman compensating potential | kcal/mol | eq:11; added to forcefield |
| $c_f$ | constant term in TMD Fixman potential | kcal/mol | bond/angle contribution |
| $\mathcal{H}_u$ | unconstrained Hamiltonian | kcal/mol | eq:3 |
| $\mathcal{Z}$ | canonical partition function | - | subscripts $u$/$c$ = unconstr/constr |
| $\rho_u,\rho_c,\rho_f$ | config pdf: unconstrained / constrained / Fixman-corrected | probability density | |
| $k$ | Boltzmann constant | kcal/(mol K) | $kT$ appears as a unit of energy |
| $T$ | temperature | K | |
| $kT$ | thermal energy | kcal/mol | inverse-temperature scale |
| $d_i$ | bond lengths | Å | $(n-1)$ of them |
| $\theta_i$, $\theta$ | bond angles | rad (or deg in results) | $(n-2)$ of them |
| $\theta_0$ | equilibrium bond angle | rad/deg | |
| $\gamma=(\gamma_1,\gamma_2,\gamma_3)$ | overall-orientation Euler angles | rad | ZXZ convention |
| $m_i$ | atomic masses | amu | |
| $\alpha_0$ | torsional barrier peak location | deg | |
| $k_\alpha$ | torsional barrier amplitude | kcal/mol | eq:15 |
| $k_\theta$ | bond-angle spring constant | kcal (paper's stated unit) | eq:21, NO 1/2 prefactor in eq:21 |
| $K_{r_i}, K_{\theta_i}$ | forcefield bond/angle spring constants | kcal/Å², kcal | eq:22, WITH 1/2 prefactor |
| $k_{\text{coul}}$ | Coulomb constant | kcal Å / e² | 332.06 |
| $q_1,q_2$ | test point charges | e | |
| $r$ | interbead distance | Å | function of torsion + angles |
| $f_{TS}$ | transition-state barrier-crossing rate | 1/time | eq:16-19 |
| $S, S_0, W, W_0, V$ | mass-matrix Schur blocks | matrices | eq:A4/A5 |
| $S^{-1}(\alpha_0)$ | $(\alpha,\alpha)$ block of $\mathcal{M}^{-1}$ | scalar for scalar $\alpha$ | eq:18b |
| $E_k$ | kinetic energy | kcal/mol | eq:A3 |
| $\phi,\psi$ | backbone dihedrals | deg | $\phi$=C–N–Cα–C, $\psi$=N–Cα–C–N |
| $H,\psi,\mathcal{K},\mathcal{D},\mathcal{P},\mathfrak{a},\mathfrak{b}$ | GNEIMO spatial operators | operators | eq:2, defined in Jain refs 21,26 |

## Sign / unit conventions (implementation-critical)

- `*` superscript = matrix/vector transpose (NOT complex conjugate).
- Fixman potential (eq:11) uses `+1/2 kT ln(detM/detM_B)`; it is ADDED to the physical forcefield. The "Fixman torque" is $-\partial U_f/\partial\alpha$.
- $\det\{\mathcal{M}_B\}$ (full BAT) is independent of torsion angles (eq:7); only the constrained $\det\{\mathcal{M}(\alpha)\}$ carries torsion dependence -> in pure TMD the Fixman potential reduces to $\tfrac12 kT\ln\det\{\mathcal{M}(\alpha)\}$ + const (eq:13).
- Bond-angle spring eq:21 is written WITHOUT a 1/2, while the forcefield form eq:22 bond/angle springs ARE written WITH a 1/2 - reconcile spring-constant conventions when porting.
- Euler angles use the ZXZ convention.
- Energies in kcal/mol, lengths Å, masses amu, angles deg in results / rad in analytics.
