# Numeric check fixtures - Ballard & Jarzynski 2009 (RENS)

## Limiting-case invariants (analytic, exact)

- Sudden limit: given $\tau = 0$, expect $w = \Delta h$ and RENS reduces exactly
  to REM ($P_{\text{acc}} = \min\{1, e^{-\Delta h}\}$, Eq. 1).
- Quasi-static limit: given a properly thermalized system and $\tau \to \infty$,
  expect $w_A = \Delta f = f_B - f_A = -w_B$, so $w = 0$ and $P_{\text{acc}} = 1$.
- Ideal gas under Eq. 13 dynamics: given ideal-gas particles switched from $T_A$
  to $T_B$, evolution exactly maps a Maxwell-Boltzmann distribution
  $T_A \to T_B$; expect $w_A + w_B = 0$ and $P_{\text{acc}} = 1$ for any $\tau$.
- Jacobian: given Eq. 13 dynamics with $N$ degrees of freedom, expect
  $J_A = (T_B/T_A)^{N/2}$ and $J_B = J_A^{-1}$, independent of initial conditions.
- Work anti-symmetry: expect $w_A(x \to x') = -w_B(\bar{x}' \to \bar{x})$
  (odd under time reversal).
- Detailed balance: RENS transition satisfies
  $P(y',x'|x,y)\,p_{AB}^{eq}(x,y) = P(\bar{x},\bar{y}|\bar{y}',\bar{x}')\,p_{AB}^{eq}(\bar{y}',\bar{x}')$.
- Scaling motivation: REM replica count grows as $M \sim N^{1/2}$ with system
  size $N$ (phase-space overlap requirement).

## Model-system test setup

- Potential: 4-well single-particle rough landscape $U(x)$ (Frenkel & Smit,
  Chapter 14). Fourth well at $x \ge 1.25$.
- Particles: $n_p = 10$, moving independently.
- Replicas: $M = 2$, $T_A = 0.30$ (primary), $T_B = 2.0$ (arbitrary units).
- Dynamics: Eq. 13 augmented Hamilton equations + Andersen thermostat.
- Work-sim attempt rate: $r = 0.166$, giving $\bar\tau_{eq} = 1/r \approx 6.0$
  (~3x the within-well relaxation rate).
- 25 test runs, $\tau$ ranging 0 to 100.
- Occupation-probability check: relative time a tagged particle spends in each
  of the 4 wells agrees (within statistical error) with the exact single-particle
  Boltzmann distribution.

## Reported results (fixtures)

- REM baseline ($f_{sw} = 0$, i.e. $\tau = 0$): given the $M=2$ setup, observed
  $\langle P_{\text{acc}} \rangle \approx 0.003$.
- REM baseline sample cost: at $f_{sw} = 0$, $t^* > 4000$.
- RENS sample cost minimum: broad minimum $t^* \sim 450-500$ for
  $f_{sw} \sim 0.2-0.6$. Highlighted run at $\tau = 2.0$ (corresponds to
  $f_{sw} \approx 0.25$ since $\tau/(\bar\tau_{eq}+\tau) = 2/(6+2) = 0.25$).
- Diminishing returns for $f_{sw} > 0.6$.
- Work-cost adjustment: observed relative CPU cost $\alpha = 2.9$ for a work step
  vs a sampling step (this model, no particle-particle interactions).
- REM replica-count sweep: $\tau = 0$, $M = 2,3,\ldots,11$, with $T_1 = 0.30$,
  $T_M = 2.0$, intermediate replicas spaced evenly in $T^{-1}$. Smallest sample
  cost $t^* = 706$ achieved at $M = 4$.
- Conclusion fixture: RENS with $M = 2$ (optimal $t^* \sim 450-500$) matches or
  beats REM with $M = 4$ ($t^* = 706$) - fewer replicas, comparable efficiency.
- Trace comparison (Fig. 5): REM and RENS output traces of $n_4(t)$ over roughly
  the same number of attempted replica exchanges ($\approx 1700$) show
  substantially different acceptance rates.

## Correlation-time definition

- $t_c = (1/\sigma^2)\int_{-\infty}^{+\infty} dt\, c(t)$, where $\sigma^2$ and
  $c(t)$ are the variance and autocorrelation of $n_4(t)$; evaluated by
  block-averaging.
