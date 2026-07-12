# Crooks (1998) - Notation and conventions

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $i_t$ | internal microstate of the system at time step $t$ | discrete label (index or config vector) | discrete time, discrete phase space assumed |
| $\lambda_t$ | externally controlled parameter at step $t$ | scalar | protocol $\{\lambda_0,\dots,\lambda_\tau\}$ is fixed and prescribed |
| $\tau$ | number of time steps in the switching protocol | integer | $t = 0,\dots,\tau$ |
| $E(i,\lambda)$ | energy of internal state $i$ at control value $\lambda$ | energy | depends on both state and control parameter |
| $\beta$ | inverse temperature | 1/energy | $\beta = 1/k_B T$ |
| $k_B$ | Boltzmann constant | energy/temperature | |
| $T$ | heat-bath temperature | temperature | held constant (isothermal) |
| $F(\beta,\lambda)$ | Helmholtz free energy | energy | $F = -\beta^{-1}\ln\sum_i e^{-\beta E(i,\lambda)}$ |
| $\Delta F$ | free energy difference | energy | $\Delta F = F(\beta,\lambda_\tau) - F(\beta,\lambda_0)$; equals reversible work $W_r$ |
| $P(A\mid\lambda)$ | canonical equilibrium probability of state $A$ | probability | $\propto e^{-\beta E(A,\lambda)}$ |
| $P(A\xrightarrow{\lambda}B)$ | single-step transition probability $A\to B$ at fixed $\lambda$ | probability | Markovian; obeys detailed balance (Eq. 8) |
| $P(A\xleftarrow{\lambda}B)$ | reverse single-step transition probability ($B\to A$ at fixed $\lambda$) | probability | reverse-time counterpart |
| $W$ | total work done ON the system over the protocol | energy | Eq. (5); path dependent |
| $W_r$ | reversible work | energy | $W_r = \Delta F$; path independent |
| $W_d$ | dissipative work | energy | $W_d = W - W_r = W - \Delta F$; $\ge 0$ on average |
| $Q$ | total heat absorbed by the system from the bath | energy | Eq. (6) |
| $\Delta E$ | total internal energy change | energy | $\Delta E = Q + W$ (Eq. 7) |
| $\overline{(\cdot)}$ | average over forward paths from a canonical initial ensemble at $\lambda_0$ | - | nonequilibrium path average |
| $\langle\cdot\rangle_0$ | equilibrium average with $\lambda$ fixed at its initial value | - | used in the fast-switching limit (Eq. 2) |
| $p$ | pressure (NPT generalization) | pressure | |
| $V(i,\lambda)$ | volume of state $i$ | volume | NPT ensemble |
| $\Delta V$ | volume change | volume | baric analogue of heat |
| $\Delta G$ | Gibbs free energy difference | energy | NPT result $\overline{e^{-\beta W}}=e^{-\beta\Delta G}$ |

## Sign / bookkeeping conventions

- Both work $W$ and heat $Q$ are defined as energy added TO the system;
  the first law reads $\Delta E = Q + W$ (Eq. 7), not $\Delta E = Q - W$.
- Work occurs on the control-parameter update substep (state $i_t$ fixed,
  $\lambda_t \to \lambda_{t+1}$). Heat occurs on the state-evolution substep
  ($\lambda$ fixed, $i_{t-1} \to i_t$). Each time step is this ordered pair.
- Entropy change of the universe per path is $\beta W_d$ (in units of $k_B$);
  entropy change of the bath is $-\beta Q$.
- Reverse-time work/heat/energy/free-energy are the negatives of the forward
  values.
- The two required dynamical properties: Markovian (memoryless) and
  microscopically reversible (single-step detailed balance, Eq. 8).
