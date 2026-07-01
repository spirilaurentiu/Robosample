# Checks / Fixtures - Betancourt 2016 (XHMC)

Concrete numbers usable as regression tests for an HMC / NUTS / XHMC port.
Most are qualitative ratios from an unseeded MCMC run, so treat "increase"
factors as order-of-magnitude expectations, not exact assertions.

## Algorithm / implementation constants

- **Divergence cutoff:** reject a trajectory around $z_0$ if
  `H(z0) - H(z) > 1000` for any state `z` in the trajectory. (Same cutoff as
  Stan/NUTS.)
- **XHMC nominal exhaustion thresholds tested:** `delta = 0.1` and `delta = 0.01`.
- **Integrator:** second-order symplectic leapfrog ($k=2$), multiplicative
  trajectory expansion (trajectory lengths $L = 2^D$).
- **NUTS state selection:** slice sampler over final trajectory.
  **XHMC state selection:** multinomial over Metropolis weights $e^{-H}/\sum e^{-H}$.
- **Detailed-balance check cost:** additive expansion -> $O(L)$ checks and $O(L)$
  stored states per length-$L$ proposal; multiplicative expansion -> $O(\log L)$
  checks and $O(\log L)$ stored states.
- **Reference implementation:** Stan/CmdStan v2.8.0, exhaustions branch,
  commit `c04d34ee77d831a2817cf3c7671aebc50a3bf825`.

## IID Gaussian target (N=100, rho=0)

- given: 100-dim IID standard Gaussian, Euclidean-Gaussian kinetic energy.
- expect: marginal energy distribution is `chi-squared with 100 dof`.
- expect: momentum-resampling energy variation is `chi-squared with 50 dof`.
- expect: every trajectory oscillates with **period $2\pi$**, independent of level
  set and initial point.
- expect: **optimal maximal integration time $T(z) = 2\pi$**, i.e.
  `L = 2*pi/epsilon ~ 64 leapfrog steps` (for the step size used).
- expect: NUTS integrates to ~**half** the optimal time (~$\pi$, ~32 steps);
  XHMC integrates much longer -> worse effective samples per transition AND per
  leapfrog step.
- virial-rate decomposition fixture: `dG/dt = 2 * sum_{n=1}^{N} (T_n - V_n)`;
  for a harmonic Gaussian each `(T_n - V_n)` oscillates about zero with period
  $2\pi$.

## Correlated Gaussian target (covariance Sigma^{ij} = rho^{|i-j|}, rho=0.95)

- given: 100-dim(?) Gaussian with `Sigma^{ij} = 0.95^{|i-j|}`.
- expect: two convergence phases — **superlinear** growth of effective sample size
  up to `~2^7 = 128 leapfrog steps`, then `sqrt(t)` asymptotic growth beyond.
- Table 1 (XHMC vs NUTS, relative to NUTS):
  - `delta = 0.1`:  total leapfrog steps `~5x`,  median ESS `~2x ~ sqrt(5)x`.
  - `delta = 0.01`: total leapfrog steps `~43x`, median ESS `~7x ~ sqrt(43)x`.
  - interpretation: ESS grows only as sqrt of the extra steps -> asymptotic,
    inefficient regime; NUTS is more efficient here.

## Nonlinear target: 1-PL item response theory model (50 students)

- model:
  - `y_i ~ Bernoulli(logistic(theta - b_i))`
  - `b_i ~ Normal(mean=0, sd=10)`
  - `theta ~ Normal(mean=0, sd=10)`
  - (normals specified by mean and standard deviation)
- property: likelihood is non-identified (data constrain only `theta` sum and
  individual `b_i`); posterior has strong nonlinear correlations.
- expect: NUTS terminates prematurely; XHMC gives larger ESS and higher efficiency.
- Table 2 (XHMC vs NUTS):
  - `delta = 0.1`:  total leapfrog steps `~13x`, median ESS `~20x` (`> 13x` ->
    superlinear, XHMC wins).
  - `delta = 0.01`: total leapfrog steps `~110x`, median ESS `~60x` (`< 110x` ->
    sublinear, diminishing returns).

## 2D graphical Gaussian test (effective potential with correlation rho)

- effective potential:
  `V_check(q) = 0.5 * q^i q^j (delta_ij - (1-delta_ij) rho)/(1-rho^2) + const`
- effective kinetic: `K_check(q,p) = 0.5 * p_i p_j delta^{ij}` (unit mass).
- `rho = 0.99`: NUTS terminates prematurely; exhaustion gives longer times for any
  delta; temporal expectations of K_check and V_check vanish far too early.
- `rho = 0.7`: NUTS no longer premature and beats the exhaustion.
