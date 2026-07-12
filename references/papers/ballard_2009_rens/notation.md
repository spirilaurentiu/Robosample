# Notation - Ballard & Jarzynski 2009 (RENS)

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $\mathcal{R}_i$ | replica $i$ (canonical ensemble) | index $1..M$ | typically ordered by increasing $T$ |
| $M$ | number of replicas | integer | |
| $H_i(x)$ | Hamiltonian/energy of replica $i$ | energy | may differ per replica or be shared |
| $T_i$ | temperature of replica $i$ | temperature | $T_A < T_B$ in illustrative case |
| $k_B$ | Boltzmann constant | energy/temperature | reduced units in model ($k_B = 1$ implied) |
| $h_i(x) = H_i/k_B T_i$ | reduced (dimensionless) Hamiltonian | dimensionless | $p_i^{eq} \propto e^{-h_i}$ |
| $x, y$ | phase/config points in replicas A, B before swap | $\mathbb{R}^N$ | $x$ in A, $y$ in B |
| $x', y'$ | configs at end of work simulations | $\mathbb{R}^N$ | $x' = x_\tau$, $y' = y_\tau$ |
| $\bar{x}$ | momentum-inverted state | $\mathbb{R}^N$ | $\mathbf{p} \to -\mathbf{p}$ |
| $\Delta h$ | REM reduced-energy swap change | dimensionless | Eq. 1 |
| $P_{\text{acc}}$ | swap acceptance probability | $[0,1]$ | $\min\{1, e^{-\Delta h}\}$ (REM) or $\min\{1,e^{-w}\}$ (RENS) |
| $\lambda$ | switching parameter | $[0,1]$ | $h(x;0)=h_A$, $h(x;1)=h_B$ |
| $\lambda_t^A$ | switching protocol in A | function of time | $\lambda_0^A=0$, $\lambda_\tau^A=1$ |
| $\lambda_t^B$ | protocol in B | | time-reversed: $\lambda_t^B = \lambda_{\tau-t}^A$ |
| $\tau$ | switching (work-simulation) time | time | free parameter; $\tau=0$ recovers REM |
| $w_A, w_B$ | reduced work in replicas A, B | dimensionless | Eqs. 4.1, 4.2 |
| $w = w_A + w_B$ | total reduced work | dimensionless | drives Eq. 2 acceptance |
| $J_A, J_B$ | propagation Jacobians | dimensionless | $J_A=|\partial x_\tau/\partial x_0|$; $\ln J$ = reduced heat |
| $M_A, M_B$ | deterministic final-state maps | | $x_\tau = M_A(x_0)$ |
| $\pi_A(x'|x)$ | work-sim arrival probability | delta function | $\delta(x'-M_A(x))$ |
| $s_\lambda$ | momentum-scaling coefficient | 1/temperature-like | $(1/2T_\lambda)\,dT_\lambda/d\lambda$ |
| $T_\lambda$ | interpolated temperature | temperature | $T_0=T_A$, $T_1=T_B$ |
| $N$ | number of degrees of freedom | integer | in Jacobian exponent $N/2$ |
| $f_i = -\ln\int dx\,e^{-h_i}$ | reduced free energy of replica $i$ | dimensionless | |
| $\bar\tau_{eq}$ | average sampling-interval duration | time | $= 1/r$ |
| $r$ | work-simulation attempt rate | 1/time | random initiation |
| $X = \tau/\bar\tau_{eq}$ | overhead ratio | dimensionless | Eq. 16 |
| $f_{sw} = X/(1+X)$ | fraction of time in work intervals | $[0,1)$ | Eq. 18 |
| $t^* = (1+X)M t_c$ | sample cost (figure of merit) | time | smaller is better; Eq. 17 |
| $t_c$ | correlation time of output trajectory | time | block-averaged |
| $\alpha$ | relative CPU cost of work vs sampling step | dimensionless | replace $X \to \alpha X$ |
| $n_p$ | number of particles in model | integer | $=10$ in tests |
| $n_4(t)$ | number of particles in fourth well at time $t$ | integer | observable |
| $U(x)$ | single-particle rough potential | energy | 4-well landscape |

## Conventions

- Reduced units throughout the model system ("arbitrary units"); temperatures
  $T_A = 0.30$, $T_B = 2.0$ are dimensionless. Treat $k_B = 1$.
- Sign of work: $w > 0$ suppresses acceptance ($e^{-w}$). Quasi-static reversible
  switching gives $w = 0$ (ideal, always accept); sudden switching gives
  $w = \Delta h$ (recovers REM).
- Time-reversal: a bar denotes momentum inversion $\mathbf{p} \to -\mathbf{p}$.
  Requires the Hamiltonian to be time-reversal symmetric $h(x;\lambda)=h(\bar{x};\lambda)$.
- $\ln J$ is the reduced heat term; for standard symplectic (volume-preserving)
  dynamics $J = 1$ and the $\ln J$ terms vanish. They are nonzero only for the
  non-volume-preserving momentum-scaling dynamics of Eq. 13.
- Protocol in B is the exact time-reverse of the protocol in A.
