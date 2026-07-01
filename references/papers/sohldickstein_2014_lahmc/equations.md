# Equations - Sohl-Dickstein 2014, LAHMC

<!-- eq:1 -->
$$p(\mathbf{x}) = \frac{1}{Z} \exp(-E(\mathbf{x})).$$
- **what:** Target (Boltzmann) distribution defined by an energy function.
- **symbols:** x - position/sample (R^N); Z - partition function (normalizer); E(x) - energy function (scalar); p(x) - target density.

<!-- eq:2 -->
$$\int p(\mathbf{x})\, T(\mathbf{x}'|\mathbf{x})\, d\mathbf{x} = p(\mathbf{x}').$$
- **what:** Fixed point (stationarity) equation the transition operator must satisfy.
- **symbols:** T(x'|x) - Markov transition density from x to x'; p - target density.

<!-- eq:3 -->
$$p(\mathbf{x})\, T(\mathbf{x}'|\mathbf{x}) = p(\mathbf{x}')\, T(\mathbf{x}|\mathbf{x}').$$
- **what:** Detailed balance condition (sufficient, not used by LAHMC).
- **symbols:** as above.

<!-- eq:4 -->
$$p(\mathbf{v}) = (2\pi)^{-\frac{N}{2}} \exp\left(-\frac{1}{2}\mathbf{v}^T\mathbf{v}\right).$$
- **what:** Auxiliary momentum prior: identity-covariance Gaussian.
- **symbols:** v - momentum (R^N); N - dimension of state space.

<!-- eq:5 -->
$$p(\zeta) = p(\mathbf{x}, \mathbf{v}) = p(\mathbf{x})\, p(\mathbf{v}) = \frac{(2\pi)^{-\frac{N}{2}}}{Z} \exp(-H(\zeta)).$$
- **what:** Joint distribution over extended state; separable in x and v.
- **symbols:** ζ = {x, v} - extended state; H(ζ) - Hamiltonian (total energy).

<!-- eq:6 -->
$$H(\zeta) = H(\mathbf{x}, \mathbf{v}) = E(\mathbf{x}) + \frac{1}{2}\mathbf{v}^{T}\mathbf{v}.$$
- **what:** Hamiltonian = potential energy + kinetic energy (unit mass).
- **symbols:** E(x) - potential energy; (1/2) v^T v - kinetic energy with mass=1.

<!-- eq:7 -->
$$\mathbf{F}\zeta = \mathbf{F}\{\mathbf{x}, \mathbf{v}\} = \{\mathbf{x}, -\mathbf{v}\}.$$
- **what:** Momentum flip operator negates the momentum.
- **symbols:** F - flip operator; x unchanged, v -> -v.

<!-- eq:8 -->
$$\mathbf{F}^{-1}\zeta = \mathbf{F}\zeta.$$
- **what:** F is its own inverse (involution).
- **symbols:** as above.

<!-- eq:9 -->
$$H(\mathbf{F}\zeta) = H(\zeta).$$
- **what:** Flip leaves total energy (hence density) unchanged.
- **symbols:** as above.

<!-- eq:10 -->
$$\left| \det \left( \frac{\partial \mathbf{F}\zeta}{\partial \zeta^T} \right) \right| = 1.$$
- **what:** Flip preserves phase-space volume (unit Jacobian determinant).
- **symbols:** Jacobian of F w.r.t. ζ.

<!-- eq:11 -->
$$\mathbf{L}\zeta = \text{state after } M \text{ leapfrog steps of Hamiltonian dynamics with step length } \epsilon.$$
- **what:** Leapfrog (Stormer-Verlet) integrator operator applied for M steps.
- **symbols:** L(ε, M) - leapfrog operator; ε - step length; M - number of leapfrog steps.

<!-- eq:12 -->
$$\mathbf{L}^{-1}\zeta = \mathbf{F}\mathbf{L}\mathbf{F}\zeta.$$
- **what:** Inverse leapfrog = flip, integrate forward, flip back (reversibility).
- **symbols:** L, F as above.

<!-- eq:13 -->
$$\left| \det \left( \frac{\partial \mathbf{L}\zeta}{\partial \zeta^T} \right) \right| = 1.$$
- **what:** Leapfrog exactly preserves phase-space volume.
- **symbols:** Jacobian of L w.r.t. ζ.

<!-- eq:14 -->
$$\mathbf{R}(\beta)\zeta = \mathbf{R}(\beta)\{\mathbf{x}, \mathbf{v}\} = \{\mathbf{x}, \mathbf{v}'\}.$$
- **what:** Momentum randomization operator; leaves x, updates v to v'.
- **symbols:** R(β) - randomization operator; β - noise fraction in [0,1].

<!-- eq:15 -->
$$\mathbf{v}' = \mathbf{v}\sqrt{1-\beta} + \mathbf{n}\sqrt{\beta}.$$
- **what:** Partial momentum refresh (Ornstein-Uhlenbeck-like mix of old momentum and noise).
- **symbols:** v' - new momentum; n - Gaussian noise; β - mixing fraction (β=1 -> full refresh, β=0 -> no change).

<!-- eq:16 -->
$$\mathbf{n} \sim N(\mathbf{0}, \mathbf{I}).$$
- **what:** Noise drawn from standard multivariate normal.
- **symbols:** n - noise vector (R^N); I - identity covariance.

<!-- eq:17 -->
$$\zeta' = \mathbf{F}\mathbf{L}\zeta^{(t,0)}.$$
- **what:** Standard HMC proposal: leapfrog then flip (self-inverse composite).
- **symbols:** ζ^(t,0) - state at step t, substep 0; ζ' - proposal.

<!-- eq:18 -->
$$\pi_{accept} = \min\left(1, \frac{p(\zeta')}{p(\zeta)}\right).$$
- **what:** Metropolis-Hastings acceptance probability (proposal probs cancel since FL is self-inverse).
- **symbols:** π_accept - acceptance probability; p(ζ) = exp(-H(ζ))/normalizer.

<!-- eq:19 -->
$$\zeta^{(t,1)} = \begin{cases} \zeta' & \text{with probability } \pi_{accept} \\ \zeta^{(t,0)} & \text{with probability } 1 - \pi_{accept} \end{cases}$$
- **what:** Accept/reject update for standard HMC.
- **symbols:** as above.

<!-- eq:20 -->
$$\zeta^{(t,2)} = \mathbf{F}\zeta^{(t,1)}.$$
- **what:** Standard-HMC momentum flip after accept/reject (causes doubling back on rejection).
- **symbols:** as above.

<!-- eq:21 -->
$$\zeta^{(t+1,0)} = \mathbf{R}(\beta)\zeta^{(t,2)}.$$
- **what:** Momentum corruption between full sampling steps.
- **symbols:** as above.

<!-- eq:22 -->
$$\beta = \alpha^{\frac{1}{\epsilon M}}.$$
- **what:** Heuristic for β so a fixed fraction α of momentum is randomized per unit simulation time.
- **symbols:** α - fraction randomized per unit time; ε - step length; M - leapfrog steps per trajectory.

<!-- eq:23 -->
$$\zeta^{(t,1)} = \begin{cases} \mathbf{L}\zeta^{(t,0)} & \text{w.p. } \pi_{\mathbf{L}^1}(\zeta^{(t,0)}) \\ \mathbf{L}^2\zeta^{(t,0)} & \text{w.p. } \pi_{\mathbf{L}^2}(\zeta^{(t,0)}) \\ \cdots & \\ \mathbf{L}^K\zeta^{(t,0)} & \text{w.p. } \pi_{\mathbf{L}^K}(\zeta^{(t,0)}) \\ \mathbf{F}\zeta^{(t,0)} & \text{w.p. } \pi_{\mathbf{F}}(\zeta^{(t,0)}) \end{cases}$$
- **what:** LAHMC transition operator: apply L 1..K times or flip F, no MH step.
- **symbols:** K - max leapfrog applications (positive integer); π_{L^a} - probability of a-fold leapfrog; π_F - flip probability.

<!-- eq:24 -->
$$\zeta^{(t+1,0)} = \mathbf{R}(\beta)\zeta^{(t,1)}.$$
- **what:** LAHMC momentum corruption between sampling steps (same as Eq 21).
- **symbols:** as above.

<!-- eq:25 -->
$$\pi_{\mathbf{L}^{a}}(\zeta) = \min\left[ 1 - \sum_{b < a} \pi_{\mathbf{L}^{b}}(\zeta),\; \frac{p(\mathbf{F}\mathbf{L}^{a}\zeta)}{p(\zeta)} \left( 1 - \sum_{b < a} \pi_{\mathbf{L}^{b}}(\mathbf{F}\mathbf{L}^{a}\zeta) \right) \right].$$
- **what:** Greedy leapfrog transition probabilities (the core LAHMC rule). First arg caps total outgoing prob at 1; second enforces forward rate <= reverse rate. Ratio p(FL^a ζ)/p(ζ) = exp(H(ζ) - H(FL^a ζ)) = exp(H(ζ) - H(L^a ζ)) since F preserves H.
- **symbols:** π_{L^a}(ζ) - prob of moving ζ -> L^a ζ; b<a - previously assigned leapfrog probs; the inner sum is evaluated at the reverse-partner state FL^a ζ.

<!-- eq:26 -->
$$p(\zeta)\, \pi_{\mathbf{L}^{a}}(\zeta) = p(\mathbf{F}\mathbf{L}^{a}\zeta)\, \pi_{\mathbf{L}^{a}}(\mathbf{F}\mathbf{L}^{a}\zeta).$$
- **what:** Generalized detailed balance relation guaranteeing the fixed point (used in the derivation).
- **symbols:** as above.

<!-- eq:27 -->
$$\pi_{\mathbf{F}}(\zeta) = 1 - \sum_{a} \pi_{\mathbf{L}^{a}}(\zeta).$$
- **what:** Residual probability assigned to the momentum flip.
- **symbols:** π_F(ζ) - flip probability; sum over a=1..K of leapfrog probs.

<!-- eq:34 -->
$$E(\mathbf{x}) = \frac{1}{2\sigma_1^2}\left(x_1^2 + x_2^2\right) + \cos\left(\frac{\pi x_1}{\sigma_2}\right) + \cos\left(\frac{\pi x_2}{\sigma_2}\right).$$
- **what:** "Rough well" 2D test energy: isotropic quadratic plus sinusoidal roughness.
- **symbols:** x1, x2 - coordinates; σ1=100 - quadratic scale; σ2=2 - sinusoid period scale.
