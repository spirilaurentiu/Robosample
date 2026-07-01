# Notation - Spiridon & Minh 2017

Reduced/convention notes: canonical ensemble, `β = 1/(k_B T)`. Cartesian mass matrix `M` is diagonal and position-independent; the generalized (internal-coordinate) mass metric tensor `M_φ` is position-dependent. Fixman potential enters the acceptance Hamiltonian only, not the guidance dynamics.

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| `φ` (phi) | generalized (internal) coordinates | R^{N_φ} | torsions/bond-angles as flexible DOF |
| `q` | Cartesian coordinates | R^{3N} | 3 per atom |
| `p` | generalized momenta conjugate to φ | R^{N_φ} | |
| `p_q` | Cartesian linear momenta | R^{3N} | |
| `M` | Cartesian mass matrix | R^{3N x 3N}, diagonal | atomic masses on diagonal, position-independent |
| `M_φ`, `M_{N_f}`, `M_{3N}` | mass metric tensor (generalized / flexible-block / full) | R^{N_φ x N_φ} etc. | `M_φ = J^T M J`, position-dependent |
| `J` | Jacobian ∂q/∂φ | R^{3N x N_φ} | `J_kl = ∂q_k/∂φ_l` |
| `|M_φ|` | determinant of mass metric tensor | scalar | O(n) evaluation algorithm exists |
| `H` | Hamiltonian | energy | `H = ½ p^T M_φ^{-1} p + U(φ)` |
| `H'` | modified (acceptance) Hamiltonian | energy | `H' = H + U_F` |
| `U(φ)` | potential energy | kcal/mol or kJ/mol | unmodified force field |
| `U_F(φ_f)` | Fixman compensating potential | energy | Eq. 3 |
| `U'(φ_f)` | modified potential | energy | `U' = U + U_F` |
| `f_p` | velocity-dependent (Coriolis/gyroscopic) forces | force | |
| `β` | inverse temperature | (energy)^{-1} | `β = (k_B T)^{-1}` |
| `k_B` | Boltzmann constant | energy/K | |
| `T` | temperature | K | sims at 300 K (also 450, 625 K) |
| `N` | number of atoms | int | |
| `N_φ` | number of generalized coordinates | int | |
| `N_f` | number of flexible coordinates | int | |
| `ε` (epsilon) | MD integrator time step | fs | 1.5 fs (unconstrained), 4 fs / up to 10 fs (constrained) |
| `η` (eta) | constraint projection tolerance | dimensionless | 10^-4 |
| `proj(·)` | projection onto constraint manifold | operator | enforces holonomic constraints |
| `N` (steps) | MD steps per HMC trajectory | int | 10 for CDHMC, 100 (or 10) for UDHMC |
| `Z(β)` | partition function | | temperature-dependent |
| `Ω(U)` | density of states at energy U | | |
| `π(U|β)` | prob. density of potential energy U | | |
| `ρ(φ)` | marginal configuration density | | Boltzmann target |
| `T(φ_j|φ_k)` | trial-attempt probability (incl. velocity draw) | | proposal kernel |
| `F(Φ,Ψ)` | free energy vs backbone dihedrals | kJ/mol (figures) | `F = -k_B T ln ρ(Φ,Ψ)` |
| `Φ, Ψ` | backbone dihedral angles | degrees | alanine dipeptide |
| `r` | uniform random number | [0,1) | `r ~ U` |
