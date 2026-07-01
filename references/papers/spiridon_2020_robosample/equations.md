# Equations - Robosample (Spiridon et al. 2020)

Constrained / generalized-coordinate statistical mechanics underlying Robosample.
Coordinates are split into free (flexible) coordinates `φ_f` and constrained
(velocity) coordinates `φ_v`. Reduced-coordinate dynamics uses only `φ_f`.

<!-- eq:1 -->
$$ \frac{d\left[\phi_{f},\phi_{v}\right]}{dt} = \frac{\partial \mathcal{H}}{\partial\left[p_{f},p_{v}\right]}, \qquad \frac{d\left[p_{f},p_{v}\right]}{dt} = -\frac{\partial \mathcal{H}}{\partial\left[\phi_{f},\phi_{v}\right]} $$
- **what:** Hamilton's equations of motion for the fully flexible system (limiting case m=0, N_f=3N).
- **symbols:** φ_f - free/flexible generalized coordinates; φ_v - constrained generalized coordinates; p_f, p_v - conjugate momenta; H - Hamiltonian; t - time.

<!-- eq:2 -->
$$ H(\phi_f, \phi_v, p_f, p_v) = [p_f, p_v]^T\, M_{tot}(\phi_f, \phi_v)\, [p_f, p_v] + U(\phi_f, \phi_v) $$
- **what:** Full Hamiltonian: kinetic term (quadratic in momenta via the mass matrix tensor) plus potential energy.
- **symbols:** M_tot - mass matrix tensor; U - potential energy; other symbols as in eq:1.
<!-- CHECK: paper prints "M_tot(φ_f, p_v)" and an unbalanced "U((φ_f,φ_v)"; corrected argument to (φ_f,φ_v). Also the kinetic form is written as p^T M_tot p (no 1/2, no inverse). Physically kinetic energy is (1/2) p^T M_tot^{-1} p; reproduce paper's literal notation but treat M_tot here as the inverse mass-matrix operator acting on momenta. -->

<!-- eq:2-def -->
$$ M_{tot} = J^T M J, \qquad J = \left\{ \frac{d x_{i}}{d \varphi_f} \right\} $$
- **what:** Definition of the mass matrix tensor: Cartesian diagonal mass matrix pushed through the coordinate-transformation Jacobian.
- **symbols:** J - Jacobian of the Cartesian-to-generalized coordinate transformation dx_i/dφ_f; M - diagonal Cartesian mass matrix; x_i - Cartesian position of atom i.
<!-- CHECK: paper writes "M_tot = JMJ" and Jacobian "{dx_ii/dφ_f}"; standard SOA form is J^T M J. -->

<!-- eq:3 -->
$$ \rho(\phi_f, \phi_v, p_f, p_v) \propto e^{-\beta H\left(\phi_f, \phi_v, p_f, p_v\right)} $$
- **what:** Boltzmann probability of a microstate in the canonical ensemble.
- **symbols:** ρ - probability density; β = 1/(kT) - inverse temperature; H - Hamiltonian from eq:2.

<!-- eq:4 -->
$$ \rho(\phi_f, \phi_v) \propto |M_{\text{tot}}(\phi)|^{-\frac{1}{2}}\, e^{-\beta U(\phi_f, \phi_v)} $$
- **what:** Marginal configuration probability for the fully flexible system, obtained by Gaussian integration over momenta. Note the mass-matrix-determinant prefactor.
- **symbols:** |M_tot(φ)| - determinant of the mass matrix tensor; U - potential energy; β as above.

<!-- eq:5 -->
$$ H(\phi_f, p_f) = p_f^T\, M(\phi_f)\, p_f + U(\phi_f) $$
- **what:** Hamiltonian for a system represented in reduced (constrained) coordinates - only free coordinates φ_f evolve.
- **symbols:** M(φ_f) - reduced-coordinate mass matrix tensor; U - potential; p_f - free momenta.
<!-- CHECK: same kinetic-form convention as eq:2 (paper omits 1/2 and inverse). -->

<!-- eq:6 -->
$$ \rho(\phi_f) \propto |M(\phi_f)|^{-\frac{1}{2}}\, e^{-\beta U(\phi_f)} $$
- **what:** Marginal configuration probability in reduced coordinates. Differs from eq:4 by the determinant used - this mismatch motivates the Fixman correction.
- **symbols:** |M(φ_f)| - determinant of reduced-coordinate mass matrix tensor.

<!-- eq:7 -->
$$ U'(\phi_f) = kT \ln \left( \frac{|M_{\text{tot}}(\phi)|}{|M(\phi_f)|} \right)^{\frac{1}{2}} $$
- **what:** Fixman correcting potential that makes the reduced-coordinate marginal (eq:6) match the fully flexible marginal (eq:4). Its gradient is the "Fixman torque" used during MD.
- **symbols:** U' - Fixman potential; kT - thermal energy; |M_tot|, |M| - the two mass-matrix-tensor determinants.

<!-- eq:8 -->
$$ H(\phi_f, p_f) = p_f^T\, M(\phi_f)\, p_f + U(\phi_f) + U'(\phi_f) $$
- **what:** Reduced-coordinate Hamiltonian including the Fixman correction, so that constrained MD/HMC samples the correct Boltzmann marginal of the flexible coordinates.
- **symbols:** U'(φ_f) - Fixman potential from eq:7; other symbols as in eq:5.

## Derivations (not implemented)
- eq:4 and eq:6 follow from performing the Gaussian momentum integral of eq:3 /
  the reduced-coordinate Boltzmann density; the `|M|^{-1/2}` prefactor is the
  determinant of the covariance of the Gaussian momentum distribution.
- eq:7 is obtained by requiring the reduced marginal (eq:6) to equal the flexible
  marginal (eq:4); the log-ratio of determinants is the entropic/kinetic correction.

## Notes / conventions
- N_f = 3N - m generalized coordinates (N atoms, m constraints); limiting case m=0 gives N_f = 3N.
- Directly inverting the mass matrix tensor is O(N^3); Spatial Operator Algebra (SOA, Rodriguez-Kreutz-Jain) gives O(N) recursive algorithms and yields the determinant, its gradient, logarithm, and square root of the mass matrix tensor.
- Metropolis-Hastings acceptance uses these energies; the Fixman potential can either enter the accept/reject criterion (HMC) or its derivative (Fixman torque) can be applied during MD.
