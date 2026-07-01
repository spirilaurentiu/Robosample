# Equations - Betancourt 2016, Optimal Integration Time / XHMC

Implementable equations for building an HMC sampler with dynamic trajectory
termination (NUTS and Exhaustive HMC). Pure differential-geometry constructions
(disintegrations, foliations) are kept in `paper.md`; only the codeable rules are
extracted here.

<!-- eq:hamiltonian -->
$$ H = K + V $$
- **what:** total Hamiltonian energy = kinetic + potential; the target on the
  cotangent bundle is $\pi_H = e^{-H}\Omega$.
- **symbols:** H - Hamiltonian (scalar); K - kinetic energy $K(q,p)$; V - potential
  energy $V(q) = -\log \pi(q) + \text{const}$.

<!-- eq:target_lift -->
$$ \pi_H = e^{-H}\Omega $$
- **what:** the lifted target distribution is the canonical (Boltzmann) measure on
  the cotangent bundle; preserved by the Hamiltonian flow.
- **symbols:** $\pi_H$ - joint density over $(q,p)$; $\Omega$ - symplectic volume form
  $\prod dq^i dp_i$.

<!-- eq:kinetic_local -->
$$ K = -\log\!\left(\frac{\mathrm{d}\xi}{dp_1 \wedge \ldots \wedge dp_n}\right),\qquad V = -\log\!\left(\frac{\mathrm{d}\pi}{dq^1 \wedge \ldots \wedge dq^n}\right) $$
- **what:** kinetic energy comes from the cotangent disintegration $\xi$; potential
  from the target $\pi$. For a Euclidean-Gaussian momentum this is the standard
  $K = \tfrac12 p^T M^{-1} p + \tfrac12\log|M|$.
- **symbols:** $\xi$ - cotangent disintegration; $\pi$ - target distribution; $p$ -
  momentum; $q$ - position.

<!-- eq:metropolis_involution -->
$$ a(z_0, R(z_L)) = \min\!\left[1,\, \exp\!\big(H(z_0) - H(R(z_L))\big)\right] $$
- **what:** Metropolis acceptance probability for a static HMC proposal made valid
  by composing the flow with an involution $R$ (usually momentum negation). Note the
  paper writes $\exp(H\circ R(z_L) - H(z_0))$ but the correct acceptance is
  $\min[1,\exp(-(H(R(z_L))-H(z_0)))]$; use the energy-decrease convention.
  <!-- CHECK: paper eq. has apparent sign as exp(H(z_L)-H(z_0)); standard HMC accept is min[1,exp(H(z_0)-H(z_L))]. Sign corrected to energy-decrease. -->
- **symbols:** $z_0$ - initial phase point $(q,p)$; $z_L$ - final point after $L$
  leapfrog steps; $R$ - involution operator; $H$ - Hamiltonian.

<!-- eq:involution_condition -->
$$ \Phi^{\tilde{H}}_{\epsilon,L\cdot\epsilon}\circ R\circ\Phi^{\tilde{H}}_{\epsilon,L\cdot\epsilon} = \mathrm{Id}_{T^*Q} $$
- **what:** condition on $R$ that makes the numerical flow a valid (reversible)
  Metropolis proposal. Momentum negation $R:(q,p)\mapsto(q,-p)$ satisfies it for
  leapfrog.
- **symbols:** $\Phi^{\tilde{H}}_{\epsilon,L\cdot\epsilon}$ - symplectic integrator, $L$
  steps of size $\epsilon$; $\mathrm{Id}$ - identity map.

<!-- eq:3 -->
$$ \mathbb{P}[z\mid\mathfrak{t}] = \frac{\frac{\mathrm{d}\pi_H}{\mathrm{d}\Omega}(z)}{\sum_{z'\in\mathfrak{t}}\frac{\mathrm{d}\pi_H}{\mathrm{d}\Omega}(z')} = \frac{e^{-H(z)}}{\sum_{z'\in\mathfrak{t}} e^{-H(z')}} $$
- **what:** (Eq. 3) Metropolis probability of selecting state $z$ from a trajectory
  $\mathfrak{t}$; the softmax of $-H$ over the trajectory states. Used for multinomial
  state selection (XHMC) and slice sampling (NUTS).
- **symbols:** $z$ - a phase point in the trajectory; $\mathfrak{t}$ - the numerical
  trajectory (set of states); $H(z)$ - Hamiltonian at $z$.

<!-- eq:4 -->
$$ \mathbb{P}[\mathfrak{t}\mid z_1] = \mathbb{P}[\mathfrak{t}\mid z_2],\quad \forall \mathfrak{t}\in\mathfrak{T}_{z_1,L}\cap\mathfrak{T}_{z_2,L}\equiv\mathfrak{T}_{(z_1,z_2),L} $$
- **what:** (Eq. 4) sufficient condition for detailed balance: the probability of
  sampling a given trajectory must be equal from any state it contains.
- **symbols:** $\mathfrak{T}_{z,L}$ - set of length-$L$ trajectories containing $z$;
  $z_1,z_2$ - two states in a shared trajectory.

<!-- eq:uniform_traj -->
$$ \mathbb{P}[\mathfrak{t}\mid z_0] = \begin{cases} 0, & \mathfrak{t}\notin\mathfrak{T}_{z_0,L}\\ 1/L, & \mathfrak{t}\in\mathfrak{T}_{z_0,L}\end{cases} $$
- **what:** uniform sampling of static-length trajectories; implemented by drawing
  $L' \sim U[0,L]$ then integrating backwards $L'$ steps and forwards $L-L'$ steps.
  Trivially satisfies Eq. 4. Equivalent to Neal's Windowed State Algorithm ($W=L$).
- **symbols:** $L$ - fixed trajectory length (leapfrog steps); $L'$ - random backward
  offset.

<!-- eq:slice -->
$$ \frac{e^{-H(z)}}{\sum_{z'\in\mathfrak{t}} e^{-H(z')}} > u,\qquad u\sim U(0,1) $$
- **what:** slice-sampler state selection: after drawing $u$, uniformly sample among
  trajectory points whose normalized weight exceeds $u$.
- **symbols:** $u$ - slice variable; other symbols as in eq:3.

<!-- eq:modified_H -->
$$ \widetilde{H} = H + \epsilon^k G + \mathcal{O}(\epsilon^{k+2}) $$
- **what:** leading-order modified (shadow) Hamiltonian that a $k$-th order symmetric
  symplectic integrator conserves exactly; explains why symplectic integrators do not
  drift. For leapfrog $k=2$.
- **symbols:** $\widetilde{H}$ - modified Hamiltonian; $\epsilon$ - step size; $k$ -
  integrator order; $G$ - leading error term (a function on phase space); $\mathcal{O}$ -
  asymptotic remainder.

<!-- eq:2 -->
$$ |\kappa(T,z)| \le \delta,\qquad \delta \in \mathbb{R}^+ $$
- **what:** (Eq. 2) generic termination criterion: stop integrating when the
  autocorrelation function falls below threshold $\delta$. Defines
  $T_\kappa(z) = \min\{ t \mid |\kappa(t,z)| \le \delta \}$.
- **symbols:** $\kappa(T,z)$ - autocorrelation function; $\delta$ - termination
  threshold; $z$ - initial phase point.

<!-- eq:kappa_u -->
$$ \kappa_u(T, z) \equiv \frac{1}{T}\int_0^T \mathrm{d}t\, \frac{\mathrm{d}u}{\mathrm{d}t}\circ\phi_t^H(z) = \frac{u\circ\phi_T^H(z) - u(z)}{T} $$
- **what:** autocorrelation function from the temporal average of the time-derivative
  of a scalar $u$; telescopes to a boundary difference over $T$. Vanishes as
  $T\to\infty$ for bounded $u$.
- **symbols:** $u$ - bounded scalar function on phase space; $\phi_t^H$ - Hamiltonian
  flow for time $t$; $T$ - integration time.

<!-- eq:virial -->
$$ G = q^i p_i $$
- **what:** the virial; the canonical scalar (besides $H$) used to build the
  exhaustion termination criterion. Its time-derivative $\mathrm{d}G/\mathrm{d}t$ drives
  the exhaustion.
- **symbols:** $G$ - virial (scalar); $q^i$ - position components; $p_i$ - momentum
  components (Einstein summation over $i$).

<!-- eq:exhaustion -->
$$ \left|\frac{1}{T_\delta}\int_0^{T_\delta} \mathrm{d}t\, \frac{\mathrm{d}G}{\mathrm{d}t}\circ\phi_t^H(z)\right| = \left|\frac{G\circ\phi_{T_\delta}^H(z) - G(z)}{T_\delta}\right| < \delta,\quad \forall z\in T^*Q $$
- **what:** (Definition 1) theoretical exhaustion criterion: terminate when the
  temporal average of the virial rate drops below $\delta$. Reduces tuning to a single
  threshold $\delta$.
- **symbols:** $T_\delta$ - exhaustion integration time; $G$ - virial; $\delta$ -
  threshold.

<!-- eq:numexp -->
$$ \left|\frac{1}{|\mathfrak{t}|}\sum_{z\in\mathfrak{t}} \mathbb{P}[z\mid\mathfrak{t}]\, \frac{\mathrm{d}G}{\mathrm{d}t}(z)\right| < \delta $$
- **what:** (Definition 2) numerical exhaustion criterion, the practical XHMC
  termination check: Metropolis-weighted average of the virial rate over the numerical
  trajectory, bounded by $\delta$. This is the CHECK_TERMINATION test for XHMC.
- **symbols:** $|\mathfrak{t}|$ - number of states in the trajectory;
  $\mathbb{P}[z\mid\mathfrak{t}]$ - Metropolis weight (eq:3); $\mathrm{d}G/\mathrm{d}t$ -
  virial rate at state $z$.

<!-- eq:nuts -->
$$ \kappa_{\mathrm{NUTS}}(T) = g_q^{-1}(p, \rho_T) < 0 $$
- **what:** generalized No-U-Turn termination criterion: stop when the momentum $p$ and
  the running position/momentum integral $\rho_T$ point in opposing directions under
  the metric $g_q^{-1}$. Euclidean metric recovers standard NUTS.
- **symbols:** $g_q^{-1}$ - inverse Riemannian metric at $q$ (Euclidean: identity);
  $p$ - current momentum; $\rho_T$ - running integral (eq:rho).

<!-- eq:rho -->
$$ \rho_T = \frac{1}{T}\int_0^T \mathrm{d}t\, \left(\phi_t^H\right)_* \theta $$
- **what:** time-averaged pushed-forward tautological one-form (position/momentum
  accumulator) used inside the No-U-Turn criterion. Euclidean case:
  $\rho_T \propto q(T) - q(0)$, giving the familiar dot-product NUTS test
  $(q^+ - q^-)\cdot p < 0$.
- **symbols:** $\theta$ - tautological one-form; $\phi_t^H$ - flow; $T$ - integration
  time.

<!-- eq:error_cutoff -->
$$ \text{reject } \mathfrak{t}\ \text{if}\ H(z_0) - H(z) > 1000,\ \forall z\in\mathfrak{t} $$
- **what:** divergence / integrator-error cutoff: reject a trajectory whose energy has
  dropped by more than 1000 (nats) relative to $z_0$ at any state. Standard Stan/NUTS
  divergence guard.
- **symbols:** $H(z_0)$ - initial energy; $H(z)$ - energy along trajectory; threshold
  = 1000.

<!-- eq:virial_iid_decomp -->
$$ \frac{\mathrm{d}G}{\mathrm{d}t} = 2\sum_{n=1}^{N} (T_n - V_n) $$
- **what:** for an IID target the virial rate decomposes per dimension; each term
  oscillates to zero at different times, so the incoherent sum delays exhaustion
  termination (explains XHMC over-integration on IID Gaussians).
- **symbols:** $N$ - dimension; $T_n$ - per-dimension kinetic energy; $V_n$ -
  per-dimension potential energy; $\mathrm{d}G/\mathrm{d}t$ - virial rate.

<!-- eq:eff_potential -->
$$ \widecheck{V}(q) = V(q) + \tfrac{1}{2}\log|g_q| + \text{const} $$
- **what:** effective potential energy for a Riemannian disintegration (log-det metric
  correction).
- **symbols:** $\widecheck{V}$ - effective potential; $V$ - base potential; $g_q$ -
  metric tensor at $q$; $|g_q|$ - determinant.

<!-- eq:eff_kinetic -->
$$ \check{K}(q,p) = A \cdot f\!\left(g_q^{-1}(p,p)\right) $$
- **what:** effective kinetic energy for a Riemannian disintegration; general
  kinetic-energy family. Euclidean-Gaussian: $A=\tfrac12$, $f=\mathrm{id}$,
  $g_q^{-1}=I$ gives $\check{K}=\tfrac12 p^T p$.
- **symbols:** $A$ - constant; $f:\mathbb{R}\to\mathbb{R}$ - scalar function; $g_q^{-1}$
  - inverse metric; $(p,p)$ - metric contraction $g_q^{-1\,ij} p_i p_j$.

<!-- eq:test_gaussian_V -->
$$ \widecheck{V}(q) = \frac{1}{2} q^i q^j \frac{\delta_{ij} - (1-\delta_{ij})\rho}{1-\rho^2} + \text{const} $$
- **what:** 2D correlated-Gaussian test target (effective potential), correlation
  $\rho$. Used in the graphical experiments.
- **symbols:** $q^i$ - position components; $\rho$ - correlation coefficient;
  $\delta_{ij}$ - Kronecker delta (NOT the threshold $\delta$).

<!-- eq:test_gaussian_K -->
$$ \check{K}(q,p) = \frac{1}{2} p_i p_j \delta^{ij} $$
- **what:** Euclidean-Gaussian kinetic energy (unit mass) used in the 2D test.
- **symbols:** $p_i$ - momentum components; $\delta^{ij}$ - Kronecker delta (identity
  mass matrix).
