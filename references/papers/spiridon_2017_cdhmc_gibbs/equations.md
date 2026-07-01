# Equations - Spiridon & Minh 2017, CDHMC as Gibbs Sampling

<!-- eq:boltzmann-joint -->
$$ \rho(\phi, p) \propto \exp\left\{-\beta\, H(\phi, p)\right\} $$
- **what:** unnormalized joint Boltzmann density over generalized coordinates and momenta (canonical ensemble).
- **symbols:** phi - generalized (internal) coordinates (R^{N_phi}); p - conjugate generalized momenta (R^{N_phi}); beta = 1/(k_B T) - inverse temperature; H - Hamiltonian; k_B - Boltzmann constant; T - temperature.

<!-- eq:hamiltonian -->
$$ H(\phi, p) = \tfrac{1}{2}\, p^{T} \mathcal{M}_{\phi}^{-1} p + U(\phi) $$
- **what:** Hamiltonian = kinetic energy (mass-metric-tensor form) plus potential energy.
- **symbols:** M_phi - mass metric tensor (R^{N_phi x N_phi}, position-dependent); U(phi) - potential energy; p - generalized momenta.

<!-- eq:mass-metric-tensor -->
$$ \mathcal{M}_{\phi} = J^{T} M J, \qquad J_{kl} = \frac{\partial q_k}{\partial \phi_l} $$
- **what:** mass metric tensor for generalized coordinates, built from the Cartesian mass matrix and the coordinate Jacobian.
- **symbols:** M - Cartesian mass matrix (R^{3N x 3N}, diagonal, position-independent, atomic masses on diagonal); J - Jacobian (R^{3N x N_phi}), J_kl = d q_k / d phi_l; q - Cartesian coordinates; phi - generalized coordinates; N - number of atoms; N_phi - number of generalized coords.

<!-- eq:kinetic-cartesian -->
$$ K_q = \tfrac{1}{2}\, p_q^{T} M^{-1} p_q $$
- **what:** Cartesian kinetic energy from linear momenta (reference form; appears in intro prose).
- **symbols:** p_q - Cartesian linear momenta (R^{3N}); M - diagonal Cartesian mass matrix (R^{3N x 3N}).

<!-- eq:2 -->
$$ \rho(\phi) \propto |\mathcal{M}_{\phi}|^{1/2}\, e^{-\beta U(\phi)} $$
- **what:** marginal configuration density obtained by Gaussian integration of Eq. joint over momenta; the mass-metric-tensor determinant makes constrained and unconstrained marginals differ.
- **symbols:** |M_phi| - determinant of the mass metric tensor; U(phi) - potential energy; beta - inverse temperature.

<!-- eq:3 -->
$$ U_F(\phi_f) = \frac{1}{2}\beta^{-1} \ln \frac{|\mathcal{M}_{N_f}|}{|\mathcal{M}_{3N}|} $$
- **what:** Fixman compensating potential; added to U so constrained-dynamics marginal matches the unconstrained marginal.
- **symbols:** phi_f - N_f flexible coordinates; M_{N_f} - mass metric tensor restricted to flexible coords (R^{N_f x N_f}, sub-block of M_{3N}); M_{3N} - full unconstrained mass metric tensor (R^{3N x 3N}); beta^{-1} = k_B T.

<!-- eq:4 -->
$$ \min\left(1,\; \frac{T(\phi_i|\phi_f)\,\rho(\phi_f)}{T(\phi_f|\phi_i)\,\rho(\phi_i)}\right) $$
<!-- CHECK: OCR printed rho(phi_o)/rho(phi_n) in numerator/denominator; standard Metropolis-Hastings HMC acceptance uses rho(phi_f)/rho(phi_i). Corrected to f (final) and i (initial). -->
- **what:** Metropolis-Hastings acceptance probability for an HMC trial move from phi_i to phi_f.
- **symbols:** T(phi_j|phi_k) - probability of attempting trial config phi_j starting from phi_k (includes probability of generating starting velocities); rho - target Boltzmann marginal density; i/f - initial/final configuration.

<!-- eq:momenta-draw -->
$$ p_0^{*} \sim \mathcal{N}\!\left(\mu = 0,\; \Sigma = k_B T \cdot \mathcal{M}^{-1}(\phi^{t-1})\right) $$
- **what:** initial momenta drawn per the internal-coordinate equipartition principle (Boltzmann-weighted kinetic energy); covariance is k_B T times the inverse mass metric tensor.
- **symbols:** p_0* - initial trial momenta; N - multivariate normal; M^{-1}(phi^{t-1}) - inverse mass metric tensor evaluated at current config phi^{t-1}; k_B T - thermal energy.

<!-- eq:vv-position -->
$$ \tilde{\phi}_{n+1} = \phi_n^{*} + \epsilon\, \nabla_p H\!\left(\phi_n^{*}, p_{n+\frac{1}{2}}^{*}\right) + \tfrac{1}{2}\epsilon\, \mathcal{M}^{-1} \nabla_\phi H\!\left(\phi_n^{*}, p_n^{*}\right) $$
<!-- CHECK: reproduced verbatim from pseudocode step 2(a); the position update mixing grad_p H and grad_phi H with M^{-1} is unusual for standard velocity Verlet - verify against Simbody constrained VV before implementing. -->
- **what:** velocity-Verlet position half-update (pre-projection) for the constrained trajectory.
- **symbols:** phi_n* - trial config at step n; tilde phi_{n+1} - unprojected next config; epsilon - MD time step; grad_p H - momentum gradient of Hamiltonian; grad_phi H - coordinate gradient of Hamiltonian; p_{n+1/2}* - half-step momentum.

<!-- eq:vv-position-proj -->
$$ \phi_{n+1}^{*} = \operatorname{proj}(\tilde{\phi}_{n+1}) $$
- **what:** project the updated configuration back onto the constraint manifold.
- **symbols:** proj(.) - projection enforcing the holonomic constraints.

<!-- eq:vv-momentum -->
$$ \tilde{p}_{n+1} = p_n^{*} - \tfrac{1}{2}\epsilon\left[ \nabla_\phi H\!\left(\phi_n^{*}, p_n^{*}\right) + \nabla_\phi H\!\left(\phi_{n+1}^{*}, p_n^{*}\right) + f_p \right] $$
- **what:** velocity-Verlet momentum update (pre-projection), including velocity-dependent Coriolis/gyroscopic forces f_p.
- **symbols:** tilde p_{n+1} - unprojected next momentum; f_p - velocity-dependent (Coriolis/gyroscopic) forces; epsilon - time step; grad_phi H - coordinate gradient of H.

<!-- eq:vv-momentum-proj -->
$$ p_{n+1}^{*} = \operatorname{proj}(\tilde{p}_{n+1}) $$
- **what:** project momenta onto the constraint manifold.
- **symbols:** proj(.) - constraint projection for momenta.

<!-- eq:vv-tolerance -->
$$ \frac{\|\tilde{p}_{n+1} - p_{n+1}^{*}\|}{\|p_{n+1}^{*}\|} > \eta \;\Rightarrow\; \text{repeat} $$
- **what:** convergence test for the constraint projection iteration; repeat until relative correction below tolerance eta.
- **symbols:** eta - constraint tolerance (10^-4 in this work).

<!-- eq:accept-pseudocode -->
$$ \min\left\{1,\; \exp\left[-\beta\left(H'(\phi_N^{*}, p_N^{*}) - H'(\phi_0^{*}, p_0)\right)\right]\right\} < (r \sim \mathcal{U}) \;\Rightarrow\; \phi_t = \phi_N^{*};\quad \text{else } \phi_t = \phi_0^{*} $$
<!-- CHECK: as printed the inequality "min{...} < r" rejects when the acceptance ratio is below r; conventional HMC accepts when r < min{1,exp(...)}. The direction of the inequality in the OCR looks reversed - verify. -->
- **what:** Metropolis-Hastings accept/reject of the constrained trajectory using the modified (Fixman-including) Hamiltonian.
- **symbols:** H' - modified Hamiltonian (Eq modified-hamiltonian); phi_N*,p_N* - endpoint of trajectory; phi_0*,p_0 - start; r ~ U - uniform(0,1) random number; phi_t - next Markov-chain sample.

<!-- eq:modified-hamiltonian -->
$$ H'(\phi, p) = H(\phi, p) + U_F(\phi) $$
- **what:** modified Hamiltonian used only in the acceptance criterion; guidance dynamics use unmodified H (no Fixman torque needed).
- **symbols:** H - unmodified Hamiltonian; U_F - Fixman potential; H' - acceptance Hamiltonian.

<!-- eq:5 -->
$$ \pi(U|\beta) = Z(\beta)^{-1}\,\Omega(U)\,\exp(-\beta U) $$
- **what:** canonical-ensemble probability of observing potential energy value U (basis of Shirts' validation test).
- **symbols:** pi(U|beta) - probability density of energy U at inverse temp beta; Z(beta) - partition function; Omega(U) - density of states at energy U; U - potential energy.

<!-- eq:6 -->
$$ \ln\left[\frac{\pi(U|\beta_2)}{\pi(U|\beta_1)}\right] = \ln\left[\frac{Z(\beta_1)}{Z(\beta_2)}\right] - (\beta_2 - \beta_1)\, U $$
- **what:** Shirts' test - log-ratio of energy histograms at two temperatures is linear in U with slope -(beta_2 - beta_1); used to validate Boltzmann sampling.
- **symbols:** beta_1, beta_2 - two inverse temperatures; Z - partition function; U - potential energy; slope = -(beta_2 - beta_1).

<!-- eq:free-energy -->
$$ F(\Phi, \Psi) = -k_B T \ln[\rho(\Phi, \Psi)] $$
- **what:** 2D free energy surface as a function of backbone dihedrals from the sampled configuration histogram.
- **symbols:** Phi, Psi - backbone dihedral angles; rho(Phi,Psi) - normalized histogram density; k_B T - thermal energy; F - free energy.
