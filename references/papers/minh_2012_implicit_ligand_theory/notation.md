# Notation - Implicit Ligand Theory (Minh 2012)

Global conventions:
- `β = 1/(k_B T)`, inverse thermal energy. Demonstration at T = 300 K.
- All free energies / energies in kcal/mol.
- Angle brackets `⟨...⟩_{X,...}^{r}` = ensemble average over coords `r` (superscript)
  w.r.t. density ∝ `q_{X,...}`; subscript `X` names the effective potential `U(r_X)`,
  extra subscripts are labels. Shorthand: functions implicitly depend on their coords
  (e.g. `Ψ ≡ Ψ(r_{RL})`).
- Hat `Â` = a statistical estimator (equation evaluated on sampled data).
- Standard concentration C° typically 1 M; the `8π²` factor is the external-orientation
  (rotational) volume for the rigid-body external DOF.

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| ΔG° | standard binding free energy | kcal/mol | R+L ⇋ RL |
| β | inverse temperature 1/(k_B T) | (kcal/mol)⁻¹ | T=300 K in demo |
| k_B | Boltzmann constant | kcal/mol/K | |
| C° | standard concentration | M | typically 1 M |
| C_X | equilibrium concentration of species X | M | X∈{R,L,RL} |
| R, L, RL | receptor, ligand, complex | label | |
| Z_{X,N} | configurational partition function in N solvent molecules | scalar | X∈{RL, R, L} |
| Z_N | pure-solvent partition function (N molecules) | scalar | |
| Z_{RL}, Z_Y | implicit-solvent configurational integrals (=Z_{·,N}/Z_N) | scalar | Y∈{R,L} |
| U(r_X, r_S) | full potential energy (solute+solvent) | kcal/mol | |
| U(r_X) | gas-phase potential of species X alone | kcal/mol | |
| U(r_S) | pure-solvent potential energy | kcal/mol | |
| ψ(r_X, r_S) | solute-solvent interaction = U(r_X,r_S)−U(r_X)−U(r_S) | kcal/mol | |
| W(r_X) | solvation PMF (gas→solvent transfer work), eq:8 | kcal/mol | PB/GB electrostatics + surface-area nonpolar |
| U(r_X) (calligraphic 𝒰) | effective potential = U(r_X)+W(r_X) | kcal/mol | implicit-solvent effective energy |
| Ψ(r_{RL}) | effective interaction energy = 𝒰(r_{RL})−𝒰(r_R)−𝒰(r_L) | kcal/mol | |
| r_{RL} | internal coords of complex | R^n | decomposes into r_R, r_L, ξ_L |
| r_R, r_L | internal coords of receptor, ligand | R^n | |
| ξ_L | 6 external DOF: relative translation (3) + rotation (3) | R^6 | |
| r_S | coords of N solvent molecules | R^{3N} | |
| I_ξ ≡ I(ξ_L) | indicator function selecting bound configs | ∈[0,1] | ≈insensitive for tight binders |
| Ω | binding-site volume = ∫ I_ξ dξ_L | (Å³·rad³ type volume) | in demo ≈ (4/3)π(0.75³)(8π²) |
| Ω_c | confined-reference volume = ∫ I_ξ e^{−βU_c} dξ_L | | |
| ΔG_ξ | ligand external-DOF confinement FE = −β⁻¹ ln(ΩC°/8π²) | kcal/mol | |
| B(r_R) | binding PMF at rigid receptor config r_R (eq:10) | kcal/mol | CENTRAL quantity |
| B̂(r_R) | estimated binding PMF | kcal/mol | |
| q_R(r_R) | receptor sampling density = e^{−β𝒰(r_R)} | unnormalized | |
| q_{L,I} | ligand sampling density = I_ξ e^{−β𝒰(r_L)} | unnormalized | |
| q_{RL,I} | complex sampling density = I_ξ e^{−β𝒰(r_{RL})} | unnormalized | |
| q_{ξ,I} | external-DOF density = I_ξ | unnormalized | |
| U_c(ξ_L) | confining potential biasing ligand orientation | kcal/mol | harmonic ⇒ Gaussian orientation |
| O(r_{RL}) | observable (energy, distance, ...) | varies | |
| Θ(r_R) | interaction-weighted rigid-receptor expectation (eq:12) | units of O | ⟨O e^{−βΨ}⟩ |
| w(r) | importance weight q_T(r)/q_S(r) | dimensionless | eq:21 |
| N | number of samples (context: ligand or receptor) | int | |
| B_cpl | vacuum coupling free energy (eq:23) | kcal/mol | MBAR-estimated |
| B_RL, B_L | vacuum→target transfer FE of complex, ligand (eq:23) | kcal/mol | single-step FEP |
| ΔU(r_X) | target−vacuum potential-energy difference = U_T(r_X)−U(r_X) | kcal/mol | |
| λ | alchemical coupling parameter | ∈[0,1] | 0 non-interacting, 1 fully interacting |
| δΨ | interaction-energy fluctuation Ψ−⟨Ψ⟩ | kcal/mol | cumulant expansion (eq:16) |
| r_E, r_I | explicit / implicit solvent coords (hybrid model) | R^{...} | supplemental |
| B'(r_R,r_E) | binding PMF at fixed receptor+explicit solvent | kcal/mol | supplemental |
