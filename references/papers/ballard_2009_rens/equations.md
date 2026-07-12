# Equations - Ballard & Jarzynski 2009 (RENS)

Reduced Hamiltonian definition (from prose):

<!-- eq:reduced_hamiltonian -->
$$ h_i(x) = \frac{H_i(x)}{k_B T_i}, \qquad p_i^{\text{eq}}(x) \propto \exp(-h_i(x)) $$
- **what:** reduced (dimensionless) Hamiltonian for replica $i$; its Boltzmann factor is the equilibrium density of that replica.
- **symbols:** $H_i$ - Hamiltonian/energy of replica $i$; $T_i$ - temperature of replica $i$; $k_B$ - Boltzmann constant; $x$ - point in configuration or phase space; $p_i^{\text{eq}}$ - equilibrium distribution.

<!-- eq:1 -->
$$ \Delta h(x,y) \equiv h_B(x) + h_A(y) - h_B(y) - h_A(x) $$
- **what:** reduced-energy change for an ordinary REM swap $x \leftrightarrow y$ between replicas A and B; used in $P_{\text{acc}} = \min\{1, e^{-\Delta h}\}$.
- **symbols:** $x$ - config in replica A before swap; $y$ - config in replica B before swap; $h_A, h_B$ - reduced Hamiltonians of the two replicas.

<!-- eq:2 -->
$$ P_{\text{acc}} = \min\{1, e^{-w}\}, \qquad w = w_A + w_B $$
- **what:** RENS work-based swap acceptance probability. Replaces the instantaneous REM criterion; if accepted, $y'$ goes to replica A and $x'$ to replica B.
- **symbols:** $w$ - total reduced work of the two work simulations; $w_A, w_B$ - reduced work in replicas A and B (Eqs. 4.1, 4.2).

<!-- eq:4.1 -->
$$ w_A(x_0 \to x_\tau) = h_B(x_\tau) - h_A(x_0) - \ln J_A(x_0) $$
- **what:** reduced work for the forward work simulation in replica A (switching $h_A \to h_B$). Analogous to first law: energy change minus a heat term $\ln J$.
- **symbols:** $x_0 = x$ - initial microstate; $x_\tau = x'$ - final microstate; $h_A, h_B$ - reduced Hamiltonians; $J_A = |\partial x_\tau/\partial x_0|$ - Jacobian of the propagation.

<!-- eq:4.2 -->
$$ w_B(y_0 \to y_\tau) = h_A(y_\tau) - h_B(y_0) - \ln J_B(y_0) $$
- **what:** reduced work for the reverse work simulation in replica B (switching $h_B \to h_A$).
- **symbols:** $y_0 = y$ - initial; $y_\tau = y'$ - final; $J_B = |\partial y_\tau/\partial y_0|$ - Jacobian.

<!-- eq:13 -->
$$ \dot{\mathbf{q}}_i = \frac{\partial H}{\partial \mathbf{p}_i}, \qquad \dot{\mathbf{p}}_i = -\frac{\partial H}{\partial \mathbf{q}_i} + \dot\lambda \, s_\lambda \, \mathbf{p}_i $$
- **what:** illustrative augmented Hamilton equations for a work simulation. Same $H$, temperatures differ; the extra momentum-scaling term heats (replica A, $\lambda:0\to1$) or cools (replica B) the system during switching.
- **symbols:** $\mathbf{q}_i, \mathbf{p}_i$ - position and momentum of particle $i$; $H$ - common Hamiltonian; $\lambda$ - switching parameter; $\dot\lambda$ - its rate; $s_\lambda$ - scaling factor (below).

<!-- eq:s_lambda -->
$$ s_\lambda = \frac{1}{2 T_\lambda} \frac{d T_\lambda}{d\lambda} $$
- **what:** momentum-scaling coefficient in Eq. 13. $T_\lambda$ interpolates from $T_0 = T_A$ to $T_1 = T_B$.
- **symbols:** $T_\lambda$ - interpolated temperature at switching parameter $\lambda$.

<!-- eq:14 -->
$$ J_A = \exp\!\left( \int_0^\tau dt \; \nabla \cdot \dot{x} \right) = \left( \frac{T_B}{T_A} \right)^{N/2}, \qquad J_B = J_A^{-1} $$
- **what:** Jacobian (phase-space volume change) for a work simulation under Eq. 13 dynamics, since $\nabla \cdot \dot{x} = N \dot\lambda s_\lambda \neq 0$. Independent of initial conditions for these dynamics.
- **symbols:** $N$ - number of degrees of freedom; $T_A, T_B$ - replica temperatures ($T_A < T_B$); $\tau$ - switching time.

<!-- eq:15 -->
$$ w_A = \Delta f = f_B - f_A = -w_B, \qquad f_i = -\ln \int dx\, e^{-h_i} $$
- **what:** quasi-static limit ($\tau \to \infty$): reduced work equals the reduced free energy difference, so $w = w_A + w_B = 0$ and $P_{\text{acc}} = 1$.
- **symbols:** $f_i$ - reduced free energy of replica $i$; $\Delta f$ - reduced free energy difference.

<!-- eq:16 -->
$$ X \equiv \tau / \bar\tau_{eq} $$
- **what:** ratio of switching time to average sampling-interval duration; overhead parameter.
- **symbols:** $\tau$ - work-simulation switching time; $\bar\tau_{eq}$ - average duration of a sampling interval ($= 1/r$ for random-initiation rate $r$).

<!-- eq:17 -->
$$ t^* = (1 + X)\, M\, t_c $$
- **what:** sample cost; total computational cost (over all M replicas) to produce one statistically independent sample in the primary replica. Figure of merit; smaller is better. For unequal per-step costs, replace $X$ by $\alpha X$.
- **symbols:** $M$ - number of replicas; $t_c$ - correlation time of the primary-replica output trajectory; $X$ - overhead ratio (Eq. 16); $\alpha$ - relative CPU cost of a work vs. sampling step.

<!-- eq:18 -->
$$ f_{sw} = \frac{X}{1 + X} = \frac{\tau}{\bar\tau_{eq} + \tau} $$
- **what:** fraction of simulation time devoted to work intervals. For unequal per-step costs, replace $X$ by $\alpha X$.
- **symbols:** as in Eqs. 16, 17.

## Derivations (proof-only; not implemented)

The detailed-balance proof (Eqs. 5-12) uses intermediate identities not needed
for implementation: the deterministic map $x_\tau = M_A(x_0)$ (Eq. 6); the
Dirac-delta arrival probabilities $\pi_A(x'|x) = \delta(x' - M_A(x))$ (Eq. 7);
the joint probability $\pi = \pi_A \pi_B$ (Eq. 8); the acceptance
$\alpha = \min\{1, e^{-w_A - w_B}\}$ (Eq. 9, same content as Eq. 2); the
time-reversal Jacobian identity $\pi_A(x'|x) = \pi_B(\bar{x}|\bar{x}')/J_A(x)$
(Eq. 10); and the work anti-symmetry $w_A(x \to x') = -w_B(\bar{x}' \to \bar{x})$
(Eq. 11). Combined, they yield the detailed-balance relation Eq. 5 /
Eq. 12. Implementers only need the transition rule $P = \pi \alpha$ and the
acceptance criterion Eq. 2.
