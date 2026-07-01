# Notation - Betancourt 2016 (Optimal Integration Time / XHMC)

Reduced/natural units throughout: energies are in **nats** (probabilities enter as
$e^{-H}$, so $\beta = 1$, i.e. $H = -\log \pi_H$). Momentum mass matrix is identity
for the Euclidean-Gaussian disintegration ($\check{K} = \tfrac12 p^T p$). Einstein
summation is used on repeated up/down indices.

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $Q$ | sample space (target manifold) | smooth $N$-dim manifold | positively oriented |
| $N$, $n$ | dimension of sample space | int | $n=N$ locally |
| $q$, $q^i$ | position | $\mathbb{R}^N$ | up index |
| $p$, $p_i$ | momentum | $\mathbb{R}^N$ | down index (cotangent) |
| $z = (q,p)$ | phase-space point | $\mathbb{R}^{2N}$ | element of $T^*Q$ |
| $T^*Q$ | cotangent bundle | $2N$-dim manifold | phase space |
| $\pi$ | target distribution on $Q$ | density | $\pi = e^{-V}$ (times volume form) |
| $\pi_H$ | lifted target on $T^*Q$ | density | $\pi_H = e^{-H}\Omega$ |
| $\pi_E$ | marginal energy distribution | density over $E$ | $\pi_E = H_*\pi$ |
| $\xi$ | cotangent disintegration | $n$-form | defines the kinetic energy |
| $\Omega$ | symplectic volume form | $2N$-form | $\prod_i dq^i dp_i$ |
| $\theta$ | tautological one-form | one-form | on $T^*Q$ |
| $H$ | Hamiltonian (total energy) | nats | $H = K + V = -\log\pi_H/\Omega$ |
| $K$, $\check{K}$ | kinetic / effective kinetic energy | nats | $\check{K}=A f(g_q^{-1}(p,p))$ |
| $V$, $\widecheck{V}$ | potential / effective potential energy | nats | $V = -\log\pi + \text{const}$ |
| $E$ | energy level-set value | nats | $H^{-1}(E)$ = level set |
| $\widetilde{H}$ | modified (shadow) Hamiltonian | nats | conserved by symplectic integrator |
| $G$ | virial | nats | $G = q^i p_i$ |
| $\mathrm{d}G/\mathrm{d}t$ | virial rate | nats/time | drives exhaustion criterion |
| $\phi_t^H$ | Hamiltonian flow for time $t$ | map $T^*Q\to T^*Q$ | exact flow |
| $\phi^H(z)$ | orbit through $z$ | subset of level set | $\{\phi_t^H(z):t\in\mathbb{R}\}$ |
| $\Phi^{\tilde H}_{\epsilon,L\epsilon}$ | symplectic integrator | map | $L$ steps of size $\epsilon$ |
| $\epsilon$ | integrator step size | time | leapfrog step |
| $L$ | number of integrator steps / trajectory length | int | $L = T/\epsilon$; $2^D$ for multiplicative |
| $D$ | tree depth (multiplicative expansion) | int | $L = 2^D$ |
| $T$, $T(z)$ | integration time | time | $\pi_{T(z)} = U(0,T(z))$ |
| $T_\delta(z)$ | exhaustion integration time | time | Definition 1 |
| $k$ | symplectic integrator order | int | leapfrog $k=2$ |
| $\mathfrak{t}$ | numerical trajectory | set of states | $|\mathfrak{t}|$ = #states |
| $\mathfrak{T}_{z,L}$ | set of length-$L$ trajectories containing $z$ | set | detailed-balance bookkeeping |
| $\mathfrak{T}_\delta$ | numerical exhaustion (set) | set | Definition 2 |
| $R$ | involution operator | map $T^*Q\to T^*Q$ | usually momentum negation $(q,p)\to(q,-p)$ |
| $\mathbb{P}[z\mid\mathfrak{t}]$ | Metropolis state weight | prob | $e^{-H(z)}/\sum e^{-H}$ (softmax of $-H$) |
| $\mathbb{P}[\mathfrak{t}\mid z]$ | trajectory sampling prob | prob | $=1/L$ for uniform static scheme |
| $a(z_0,z_L)$ | acceptance probability | prob | $\min[1,e^{-\Delta H}]$ |
| $\kappa(T,z)$ | autocorrelation / termination function | scalar | terminate when $|\kappa|\le\delta$ |
| $\kappa_u$ | virial-style autocorrelation from scalar $u$ | scalar | $=(u\circ\phi_T^H - u)/T$ |
| $\kappa_{\mathrm{NUTS}}$ | No-U-Turn criterion | scalar | terminate when $<0$ |
| $\rho_T$ | running (time-averaged) one-form integral | vector | NUTS accumulator |
| $\delta$ | exhaustion termination threshold | $\mathbb{R}^+$ | nominal 0.1 or 0.01 |
| $\delta_{ij}$ | Kronecker delta | - | NOT the threshold $\delta$ |
| $g$, $g_q$ | Riemannian metric | $N\times N$ SPD | $g_q^{-1}$ = inverse (mass matrix) |
| $\rho$ | Gaussian correlation coefficient | scalar | test targets |
| $\Sigma^{ij}$ | covariance of correlated Gaussian | $N\times N$ | $\Sigma^{ij}=\rho^{|i-j|}$ |
| $u$ | slice variable / scalar function | $U(0,1)$ or scalar | context dependent |
| $f$ | expectation integrand / kinetic function | scalar | context dependent |
