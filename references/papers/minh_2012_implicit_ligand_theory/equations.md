# Equations - Implicit Ligand Theory (Minh 2012)

Sign/unit convention: `β = 1/(k_B T)`. Free energies in the demonstration are in
kcal/mol at T = 300 K. Angle brackets `⟨...⟩_{X,...}^{r}` denote an ensemble
average over the coordinates `r` in the superscript with respect to the density
proportional to `q_{X,...}`, where `X` labels the effective potential
`U(r_X)` and `...` are extra labels.

<!-- eq:1 -->
$$ \Delta G^{\circ} = -\beta^{-1} \ln \left( \frac{C^{\circ} C_{RL}}{C_R C_L} \right) $$
- **what:** standard binding free energy from equilibrium concentrations of R, L, RL.
- **symbols:** ΔG° - standard binding free energy (kcal/mol); β=1/(k_B T) - inverse temperature; C° - standard concentration (typically 1 M); C_X - equilibrium concentration of species X∈{R,L,RL}.

<!-- eq:2 -->
$$ \Delta G^{\circ} = -\beta^{-1} \ln \left( \frac{Z_{RL,N} Z_N}{Z_{R,N} Z_{L,N}} \frac{C^{\circ}}{8\pi^2} \right) $$
- **what:** binding free energy as a ratio of configurational partition functions (symmetry numbers and a small PV term omitted).
- **symbols:** Z_{RL,N}, Z_{Y,N}, Z_N - configurational partition functions of complex, species Y, and pure solvent (N solvent molecules); 8π² - rotational volume factor for external DOF.

<!-- eq:3 -->
$$ Z_{RL,N} = \int I_{\xi}\, e^{-\beta U(r_{RL}, r_S)}\, dr_{RL}\, dr_S $$
- **what:** configurational partition function of the complex in N solvent molecules.
- **symbols:** I_ξ ≡ I(ξ_L) - indicator function ∈[0,1] selecting bound configurations; r_{RL} - internal coords of complex; r_S - coords of N solvent molecules; U - potential energy.

<!-- eq:4 -->
$$ Z_{Y,N} = \int e^{-\beta U(r_Y, r_S)}\, dr_Y\, dr_S $$
- **what:** configurational partition function of species Y∈{R,L} in N solvent molecules.
- **symbols:** r_Y - internal coords of receptor or ligand alone.

<!-- eq:5 -->
$$ Z_N = \int e^{-\beta U(r_S)}\, dr_S $$
- **what:** configurational partition function of pure solvent (N molecules).
- **symbols:** U(r_S) - potential energy of solvent by itself.

<!-- eq:6 -->
$$ Z_{RL} \equiv \frac{Z_{RL,N}}{Z_N} = \int I_{\xi}\, e^{-\beta [U(r_{RL}) + W(r_{RL})]}\, dr_{RL} $$
- **what:** implicit-solvent configurational integral for the complex (solvent integrated out).
- **symbols:** W(r_X) - solvation PMF (eq:8); U(r_X) - gas-phase potential energy of species X.

<!-- eq:7 -->
$$ Z_Y \equiv \frac{Z_{Y,N}}{Z_N} = \int e^{-\beta [U(r_Y) + W(r_Y)]}\, dr_Y $$
- **what:** implicit-solvent configurational integral for species Y∈{R,L}.
- **symbols:** as above.

<!-- eq:8 -->
$$ W(r_X) = -\beta^{-1} \ln \left( \frac{\int e^{-\beta \psi(r_X, r_S)}\, e^{-\beta U(r_S)}\, dr_S}{\int e^{-\beta U(r_S)}\, dr_S} \right) $$
- **what:** solvation potential of mean force; reversible work of transferring species X from gas phase into solvent. Estimated in practice as Poisson-Boltzmann (or Generalized Born) electrostatics + surface-area nonpolar term.
- **symbols:** ψ(r_X,r_S)=U(r_X,r_S)−U(r_X)−U(r_S) - solute-solvent interaction energy.

<!-- eq:9 -->
$$ \Delta G^{\circ} = -\beta^{-1} \ln \left( \frac{Z_{RL}}{Z_R Z_L} \frac{C^{\circ}}{8\pi^2} \right) $$
- **what:** binding free energy in terms of implicit-solvent configurational integrals.
- **symbols:** Z_{RL}, Z_R, Z_L - implicit-solvent integrals (eqs 6,7).

<!-- eq:10 -->
$$ B(r_R) = -\beta^{-1} \ln \left( \frac{\int I_{\xi}\, e^{-\beta \Psi(r_{RL})}\, e^{-\beta \mathcal{U}(r_L)}\, dr_L\, d\xi_L}{\int I_{\xi}\, e^{-\beta \mathcal{U}(r_L)}\, dr_L\, d\xi_L} \right) \equiv -\beta^{-1} \ln \left\langle e^{-\beta \Psi} \right\rangle_{L,I}^{r_L,\xi_L} $$
- **what:** binding PMF for a rigid receptor configuration r_R (CENTRAL RESULT). Exponential average of the effective interaction energy over ligand internal + external coords, at fixed receptor.
- **symbols:** U(r_X)=U(r_X)+W(r_X) - effective (implicit-solvent) potential; Ψ(r_{RL})=U(r_{RL})−U(r_R)−U(r_L) - effective interaction energy; ξ_L - 6 external DOF (relative translation+rotation); q_{L,I}(r_L,ξ_L)=I_ξ e^{−βU(r_L)} - sampling density; average taken over r_L,ξ_L.

<!-- eq:11 -->
$$ \Delta G^{\circ} = -\beta^{-1} \ln \left( \frac{\int I_{\xi}\, e^{-\beta \mathcal{U}(r_{RL})}\, dr_{RL}}{\int e^{-\beta \mathcal{U}(r_{R})}\, dr_{R} \int e^{-\beta \mathcal{U}(r_{L})}\, dr_{L}} \frac{C^{\circ}}{8\pi^{2}} \right) = -\beta^{-1} \ln \left( \frac{\int e^{-\beta [B(r_{R}) + \mathcal{U}(r_{R})]}\, dr_{R}}{\int e^{-\beta \mathcal{U}(r_{R})}\, dr_{R}} \frac{\Omega C^{\circ}}{8\pi^{2}} \right) \equiv -\beta^{-1} \ln \left\langle e^{-\beta B} \right\rangle_{R}^{r_{R}} + \Delta G_{\xi} $$
<!-- CHECK: first-form numerator integrand corrected from OCR U(r_R) to U(r_RL) per Eq (11) line 1 in prose -->
- **what:** binding free energy as an exponential average of the binding PMF over receptor configurations sampled from q_R (CENTRAL RESULT). Separates receptor sampling from ligand sampling.
- **symbols:** Ω=∫I_ξ dξ_L - binding-site volume; ΔG_ξ=−β⁻¹ ln(ΩC°/8π²) - free energy of confining ligand external DOF to the binding site; q_R(r_R)=e^{−βU(r_R)} - receptor sampling density; B(r_R) from eq:10.

<!-- eq:12 -->
$$ \Theta(r_R) = \frac{\int I_{\xi}\, O(r_{RL})\, e^{-\beta \Psi(r_{RL})}\, e^{-\beta \mathcal{U}(r_L)}\, dr_L\, d\xi_L}{\int I_{\xi}\, e^{-\beta \mathcal{U}(r_L)}\, dr_L\, d\xi_L} \equiv \left\langle O\, e^{-\beta \Psi} \right\rangle_{L,I}^{r_L,\xi_L} $$
- **what:** rigid-receptor expectation of observable O weighted by the interaction-energy Boltzmann factor. If O depends only on r_R it reduces to O(r_R) e^{−βB(r_R)}.
- **symbols:** O(r_{RL}) - observable (e.g. mean energy, interaction energy, atom-atom distance); Θ(r_R) - interaction-weighted rigid-receptor expectation.

<!-- eq:13 -->
$$ \langle O \rangle_{RL,I}^{r_{RL}} = \frac{\int \Theta(r_R)\, e^{-\beta \mathcal{U}(r_R)}\, dr_R}{\int e^{-\beta [B(r_R) + \mathcal{U}(r_R)]}\, dr_R} = \frac{\langle \Theta \rangle_R^{r_R}}{\langle e^{-\beta B} \rangle_R^{r_R}} = \langle \Theta \rangle_R^{r_R}\, e^{\beta [\Delta G^{\circ} - \Delta G_{\xi}]} $$
- **what:** bound-ensemble expectation of an observable via implicit ligand theory (generalizes implicit ligand sampling).
- **symbols:** q_{RL,I}(r_{RL})=I_ξ e^{−βU(r_{RL})} - bound-complex density; other symbols as above.

<!-- eq:14 -->
$$ B(r_R) = -\beta^{-1} \ln \left( \frac{\int I_{\xi}\, e^{-\beta \mathcal{U}(r_{RL})}\, dr_L\, d\xi_L}{\int I_{\xi}\, e^{-\beta [\mathcal{U}(r_L) + \mathcal{U}(r_R)]}\, dr_L\, d\xi_L} \right) $$
- **what:** binding PMF as a rigid-receptor free energy difference; computable by FEP, TI, or BAR.
- **symbols:** as above; r_R held rigid.

<!-- eq:15 -->
$$ \hat{B}(r_R) = -\beta^{-1} \ln \frac{1}{N} \sum_{n=1}^{N} e^{-\beta \Psi(r_{RL,n})} $$
- **what:** forward FEP sample-mean estimator of the binding PMF. Ligand coords sampled from q_L; external ξ_L resampled from q_{ξ,I}=I_ξ.
- **symbols:** N - number of complex samples; r_{RL,n} - n-th sample; hat denotes estimator.

<!-- eq:16 -->
$$ B(r_R) \approx \langle \Psi \rangle_{L,I}^{r_L,\xi_L} - \frac{\beta}{2!} \left\langle (\delta\Psi)^2 \right\rangle_{L,I}^{r_L,\xi_L} + \frac{\beta^2}{3!} \left\langle (\delta\Psi)^3 \right\rangle_{L,I}^{r_L,\xi_L} - \frac{\beta^3}{4!} \left[ \left\langle (\delta\Psi)^4 \right\rangle_{L,I}^{r_L,\xi_L} - 3 \left( \left\langle (\delta\Psi)^2 \right\rangle_{L,I}^{r_L,\xi_L} \right)^2 \right] $$
<!-- CHECK: OCR "\xi\Psi" is the fluctuation δΨ=Ψ−⟨Ψ⟩ (defined right after Eq 16); rendered here as δΨ -->
- **what:** fourth-order cumulant expansion of the binding PMF exponential average.
- **symbols:** δΨ = Ψ(r_{RL}) − ⟨Ψ⟩_{L,I}^{r_L,ξ_L} - interaction-energy fluctuation.

<!-- eq:17 -->
$$ B(r_R) = -\beta^{-1} \ln \left\langle e^{-\beta[\Psi - U_c]} \right\rangle_{L,I_c}^{r_L,\xi_L} - \beta^{-1} \ln \left( \frac{\Omega_c}{\Omega} \right) $$
- **what:** binding PMF estimated with a confining potential U_c(ξ_L) biasing the ligand orientation toward favorable poses.
- **symbols:** U_c(ξ_L) - confining potential on external DOF; Ω_c=∫I_ξ e^{−βU_c(ξ_L)} dξ_L; q_{L,Ic}=I_ξ e^{−β[U(r_L)+U_c(ξ_L)]} - biased sampling density.

<!-- eq:18 -->
$$ B(r_R) = \beta^{-1} \ln \left( \frac{\int I_{\xi}\, e^{\beta \Psi(r_{RL})}\, e^{-\beta \mathcal{U}(r_{RL})}\, dr_L\, d\xi_L}{\int I_{\xi}\, e^{-\beta \mathcal{U}(r_{RL})}\, dr_L\, d\xi_L} \right) = \beta^{-1} \ln \left\langle e^{\beta \Psi} \right\rangle_{RL,I}^{r_L,\xi_L} $$
- **what:** inverse (reverse-FEP) form of the binding PMF; sample ligand from the fully interacting rigid-receptor complex ensemble.
- **symbols:** q_{RL,I}(r_L,ξ_L)=I_ξ e^{−βU(r_{RL})} - fully-interacting sampling density.

<!-- eq:19 -->
$$ B(r_R) = \beta^{-1} \ln \left\langle e^{\beta[\Psi - U_c]} \right\rangle_{RL,I}^{r_L,\xi_L} - \beta^{-1} \ln \frac{\Omega_c}{\Omega} $$
- **what:** inverse form with a confining reference state, alleviating phase-space overlap problems.
- **symbols:** as eqs 17-18.

<!-- eq:20 -->
$$ \Delta \hat{G}^{\circ} = -\beta^{-1} \ln \frac{1}{N} \sum_{n=1}^{N} e^{-\beta \hat{B}(r_{R,n})} + \Delta G_{\xi} $$
- **what:** binding free energy sample-mean estimator when receptor configs are drawn from q_R.
- **symbols:** B-hat(r_{R,n}) - estimated binding PMF for n-th of N receptor configs.

<!-- eq:21 -->
$$ \langle O \rangle_T = \frac{\int O(r) q_T(r)\, dr}{\int q_T(r)\, dr} = \frac{\int O(r) w(r) q_S(r)\, dr}{\int w(r) q_S(r)\, dr} = \frac{\langle wO \rangle_S}{\langle w \rangle_S} $$
- **what:** generic importance-sampling identity relating a target-distribution expectation to a sampling-distribution expectation.
- **symbols:** w(r)=q_T(r)/q_S(r) - ratio of unnormalized target/sampling densities; subscripts T,S - target and sampling distributions.

<!-- eq:22 -->
$$ \Delta \hat{G}^{\circ} = -\beta^{-1} \ln \frac{\sum_{n=1}^{N} w(r_{R,n})\, e^{-\beta \hat{B}(r_{R,n})}}{\sum_{n=1}^{N} w(r_{R,n})} + \Delta G_{\xi} $$
- **what:** binding free energy estimator when receptor configs are drawn from a biased distribution (importance-weighted).
- **symbols:** w(r_{R,n})=q_R(r_{R,n})/q_{R,w}(r_{R,n}) - importance weight.

<!-- eq:23 -->
$$ B(r_{R}) = B_{cpl} + B_{RL} - B_{L} - \Delta U(r_{R}) $$
$$ B_{cpl} = -\beta^{-1} \ln \left( \frac{\int I_{\xi}\, e^{-\beta U(r_{RL})}\, dr_{L}\, d\xi_{L}}{\int I_{\xi}\, e^{-\beta [U(r_{L}) + U(r_{R})]}\, dr_{L}\, d\xi_{L}} \right) $$
$$ B_{RL} = -\beta^{-1} \ln \left( \frac{\int I_{\xi}\, e^{-\beta \Delta U(r_{RL})}\, e^{-\beta U(r_{RL})}\, dr_{L}\, d\xi_{L}}{\int I_{\xi}\, e^{-\beta U(r_{L})}\, dr_{L}\, d\xi_{L}} \right) $$
<!-- CHECK: B_RL denominator uses U(r_L) with the RL-numerator Boltzmann factor per printed Eq (23); transcribed as printed -->
$$ B_{L} = -\beta^{-1} \ln \left( \frac{\int I_{\xi}\, e^{-\beta \Delta U(r_{L})}\, e^{-\beta U(r_{L})}\, dr_{L}\, d\xi_{L}}{\int I_{\xi}\, e^{-\beta U(r_{L})}\, dr_{L}\, d\xi_{L}} \right) $$
- **what:** decomposition of the binding PMF used because alchemical coupling was run in vacuum: turn on ligand-receptor coupling in vacuum (B_cpl), then transfer complex, ligand, receptor from vacuum to target (implicit-solvent) state.
- **symbols:** ΔU(r_X)=U_T(r_X)−U(r_X) - potential-energy difference between target state and the (vacuum) state configs were sampled from; U - gas-phase potential; B_cpl estimated by MBAR, B_RL and B_L by single-step FEP.

<!-- eq:24 -->
$$ \hat{\Theta}(r_R) = \frac{\sum_{n=1}^{N} w(r_{RL,n})\, O(r_{RL,n})}{\sum_{n=1}^{N} w(r_{RL,n})} $$
$$ w(r_{RL}) = \frac{e^{-\beta(U_{PBSA}(r_L) - U_0(r_L))}}{1 + \frac{N_1}{N_0}\, e^{-\beta(U_1(r_{RL}) - \hat{B}_{cpl} - U_0(r_{RL}))}} $$
<!-- CHECK: numerator arg printed as U_PBSA(R_L)-U_0(R_L); read as U_PBSA(r_L)-U_0(r_L) -->
- **what:** MBAR estimator of the rigid-receptor interaction-weighted expectation, combining snapshots from non-interacting and fully-interacting states.
- **symbols:** U_0, U_1 - potential energies of non-interacting and fully-interacting complex; U_PBSA(r_L) - PBSA energy of ligand only; N_0, N_1 - snapshot counts of non-interacting / fully-interacting complex; B̂_cpl - MBAR estimate of coupling free energy.

<!-- eq:25 -->
$$ \mathrm{RMSE}(m_1, m_2) = \sqrt{\frac{1}{L} \sum_{l=1}^{L} \left( \Delta G_{l,m_1}^{\circ} - \Delta G_{l,m_2}^{\circ} \right)^2} $$
- **what:** root-mean-square error between binding free energies from two methods across the ligand set.
- **symbols:** L - number of ligands; ΔG°_{l,m} - binding FE estimate for ligand l by method m.

## Supplemental: hybrid implicit-explicit solvent

<!-- eq:26 -->
$$ Z_{RL'} = \int I_{\delta}\, e^{-\beta [U(r_{RL}, r_E) + W(r_{RL}, r_E)]}\, dr_{RL}\, dr_E $$
- **what:** complex partition function with a few explicit solvent molecules treated as part of the receptor.
- **symbols:** r_E - explicitly represented solvent coords; r_I - implicitly represented (integrated); I_δ - indicator function on external DOF δ_L.

<!-- eq:27 -->
$$ Z_{R'} = \int e^{-\beta [U(r_R, r_E) + W(r_R, r_E)]}\, dr_R\, dr_E $$
- **what:** receptor+explicit-solvent partition function.
- **symbols:** as above.

<!-- eq:28 -->
$$ \Delta G^{\circ} = -\beta^{-1} \ln \left( \frac{Z_{RL'}}{Z_R' Z_L} \frac{C^{\circ}}{8\pi^2} \right) = -\beta^{-1} \ln \left( \frac{\int e^{-\beta [B'(r_R, r_E) + \mathcal{U}(r_R, r_E)]}\, dr_R\, dr_E}{\int e^{-\beta \mathcal{U}(r_R, r_E)}\, dr_R\, dr_E} \frac{\Omega C^{\circ}}{8\pi^2} \right) = -\beta^{-1} \ln \left\langle e^{-\beta B'} \right\rangle_{R,E}^{r_R,r_E} + \Delta G_{\delta} $$
- **what:** binding free energy in the hybrid implicit-explicit solvent model.
- **symbols:** B'(r_R,r_E) - binding PMF at fixed receptor + explicit-solvent config; q_{R,E}=e^{−βU(r_R,r_E)}; ΔG_δ - external-DOF confinement free energy.
