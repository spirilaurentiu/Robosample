# Crooks (1998) - Check fixtures

This is an analytic paper: it states no result tables or benchmark numbers. The
testable content is a set of exact identities and limits that any correct work-
based free-energy estimator must satisfy. Use these as regression oracles.

## Exact identities

- **Work equality (Eq. 1).** Given a set of forward switching paths sampled with
  a canonical initial ensemble at $\lambda_0$ and a fixed protocol
  $\{\lambda_0,\dots,\lambda_\tau\}$: expect
  $\overline{e^{-\beta W}} = e^{-\beta \Delta F}$, i.e.
  $\Delta F = -\beta^{-1}\ln \overline{e^{-\beta W}}$.
  Regression: on an analytically solvable model (e.g. moving-harmonic-well or a
  two-state Ising switch), the estimator must converge to the exact $\Delta F$.

- **First law (Eq. 7).** Given per-step work (Eq. 5) and heat (Eq. 6):
  expect $Q + W = E(i_\tau,\lambda_\tau) - E(i_0,\lambda_0)$ exactly, path by path.
  Regression: assert this equality every trajectory; sign convention $\Delta E = Q + W$.

- **Generalized detailed balance (Eq. 9).** For any forward path and its time
  reverse under a microscopically reversible Markov kernel:
  expect $P_{\text{fwd}}/P_{\text{rev}} = e^{-\beta Q}$.

- **Path-ensemble ratio (Eq. 10).** With equilibrium endpoints:
  expect forward/reverse joint-path probability ratio $= e^{\beta W_d}$,
  where $W_d = W - \Delta F$.

## Limits

- **Reversible / infinitely-slow limit.** As switching time $\to\infty$ the process
  is reversible: expect $W_d \to 0$, so $W \to \Delta F$ for every path and
  $\overline{W} = \Delta F$ (thermodynamic integration).

- **Infinitely-fast / instantaneous limit (Eq. 2).** With $\tau = 1$ and the
  configuration held fixed during the $\lambda$ switch: expect
  $\langle e^{-\beta W}\rangle_0 = e^{-\beta \Delta F}$ with
  $W = E(i,\lambda_1) - E(i,\lambda_0)$ (Zwanzig thermodynamic perturbation).

- **Jensen inequality (second law).** From Eq. (1) by convexity: expect
  $\overline{W} \ge \Delta F$, equality only in the reversible limit; the average
  dissipative work $\overline{W_d} \ge 0$.

- **NPT generalization (Eq. 13).** In the isothermal-isobaric ensemble the same
  estimator returns Gibbs free energy: expect
  $\overline{e^{-\beta W}} = e^{-\beta \Delta G}$.

## Parameters / constants stated in the paper

- $\beta = 1/k_B T$ (definition; the only "constant" the paper introduces).
- No numeric hyperparameters, energies, or benchmark values are reported.
