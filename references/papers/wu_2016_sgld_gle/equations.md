# Equations - Wu, Brooks, Vanden-Eijnden 2016 (SGLD-GLE)

Reduced/CGS-free molecular units. $k$ = Boltzmann constant, $T$ = temperature,
$\gamma$ = collision (friction) frequency, $t_L$ = local averaging time,
$\lambda$ = guiding factor, $\mu$ = GLE guiding parameter with
$\lambda = \mu(2-\mu)$.

<!-- eq:1 -->
$$ \dot{\mathbf{p}}_i = \mathbf{f}_i + \mathbf{g}_i - \gamma \mathbf{p}_i + \mathbf{R}_i $$
- **what:** SGLD equation of motion for particle $i$ (standard LD plus guiding force $\mathbf{g}_i$).
- **symbols:** $\dot{\mathbf{p}}_i$ - time derivative of momentum ($\mathbb{R}^3$); $\mathbf{f}_i$ - interaction (conservative) force; $\mathbf{g}_i$ - guiding force; $\gamma$ - collision frequency (scalar, units 1/time); $\mathbf{p}_i$ - momentum; $\mathbf{R}_i$ - white Gaussian random force.

<!-- eq:2 -->
$$ \langle \mathbf{R}_j(0)\,\mathbf{R}_i(t) \rangle = 2 m_i kT \gamma\, \delta(t)\,\delta_{ij} $$
- **what:** covariance (fluctuation-dissipation) of the SGLD white random force.
- **symbols:** $m_i$ - mass of particle $i$; $kT$ - thermal energy; $\delta(t)$ - Dirac delta in time; $\delta_{ij}$ - Kronecker delta over particles.

<!-- eq:3 -->
$$ \mathbf{g}_i(t) = \lambda \gamma\, \widetilde{\mathbf{p}}_i(t) - \xi \gamma\, \mathbf{p}_i(t) $$
- **what:** SGLD guiding force: accelerate low-frequency ($+\lambda\gamma\tilde{\mathbf p}$), damp high-frequency ($-\xi\gamma\mathbf p$).
- **symbols:** $\lambda$ - guiding factor (dimensionless, strength); $\widetilde{\mathbf{p}}_i$ - low-frequency (locally averaged) momentum; $\xi$ - energy-conservation factor (eq. 6).

<!-- eq:4 -->
$$ \widetilde{P}(t) = \frac{1}{t_L}\int_{-\infty}^{t} P(\tau)\, e^{-\frac{t-\tau}{t_L}}\, d\tau \approx \left(1 - \frac{\delta t}{t_L}\right)\widetilde{P}(t-\delta t) + \frac{\delta t}{t_L} P(t) $$
- **what:** exponential local time-average (low-pass filter) of any quantity $P$; right form is the incremental update used per timestep.
- **symbols:** $t_L$ - local averaging time (units of time); $\delta t$ - integration timestep; $P(\tau)$ - quantity being averaged (e.g. momentum $\mathbf{p}_i$).

<!-- eq:5 -->
$$ \sum_i \mathbf{g}_i \cdot \dot{\mathbf{r}}_i = \lambda\gamma \sum_i \widetilde{\mathbf{p}}_i \cdot \dot{\mathbf{r}}_i - \xi\gamma \sum_i \mathbf{p}_i \cdot \dot{\mathbf{r}}_i = 0 $$
- **what:** zero-net-work condition determining $\xi$; guiding force does no net work.
- **symbols:** $\dot{\mathbf{r}}_i$ - velocity of particle $i$; sum runs over all particles.

<!-- eq:6 -->
$$ \xi = \frac{\lambda \sum_i \widetilde{\mathbf{p}}_i \cdot \dot{\mathbf{r}}_i}{\sum_i \mathbf{p}_i \cdot \dot{\mathbf{r}}_i} $$
- **what:** energy-conservation factor solved from eq. (5); recomputed each step.
- **symbols:** as above; $\xi$ - scalar per configuration.

<!-- eq:7 -->
$$ \dot{\mathbf{p}}_i = \mathbf{f}_i - \gamma \int_{-\infty}^{t} d\tau\, K(t-\tau)\, \mathbf{p}_i(\tau) + \eta_i $$
- **what:** generalized Langevin equation (GLE) form with memory kernel $K$ and colored noise $\eta_i$.
- **symbols:** $K(t-\tau)$ - memory kernel (eq. 9); $\eta_i(t)$ - zero-mean Gaussian (colored) noise (eq. 11).

<!-- eq:8 -->
$$ \langle \eta_i(t)\,\eta_j(t') \rangle = \delta_{ij}\, m_i kT \gamma\, K(t-t') $$
- **what:** fluctuation-dissipation theorem relating colored noise covariance to the kernel.
- **symbols:** $K(t-t')$ - memory kernel; other symbols as above.

<!-- eq:9 -->
$$ K(t) = 2\delta(t) - \frac{\lambda}{t_L} e^{-\frac{t}{t_L}} $$
- **what:** memory kernel that makes the GLE dissipation reproduce SGLD (neglecting the $\xi$ term).
- **symbols:** $\delta(t)$ - Dirac delta; $\lambda$ - guiding factor; $t_L$ - averaging time.

<!-- eq:10 -->
$$ \int_{-\infty}^{t} \delta(t-\tau)\, P(\tau)\, d\tau = \frac{1}{2} P(t) $$
- **what:** convention for the half-sided delta integral (upper limit at the delta location gives 1/2).
- **symbols:** $P(\tau)$ - integrand quantity.

<!-- eq:11 -->
$$ \eta_i(t) = \mathbf{R}_i(t) - \frac{\mu}{t_L}\int_{-\infty}^{t} \mathbf{R}_i(\tau)\, e^{-\frac{t-\tau}{t_L}}\, d\tau = \mathbf{R}_i(t) - \mu\, \widetilde{\mathbf{R}}_i(t) $$
- **what:** colored GLE noise built from the white noise $\mathbf{R}_i$ minus $\mu$ times its local average $\widetilde{\mathbf R}_i$.
- **symbols:** $\mu \ge 0$ - GLE guiding parameter (eq. 13); $\widetilde{\mathbf{R}}_i$ - locally averaged random force (same filter as eq. 4). <!-- CHECK: raw eq.(11) printed the exponent as e^{-\tau/t_L}; corrected to e^{-(t-\tau)/t_L} to match the definition of the local average in eq.(4) and the second equality -->

<!-- eq:13 -->
$$ \lambda = \mu\,(2 - \mu) $$
- **what:** relation between guiding factor $\lambda$ and GLE parameter $\mu$; two roots give statistically equivalent noise.
- **symbols:** $\lambda$ - guiding factor; $\mu$ - GLE parameter. Roots: $\mu = 1 \pm \sqrt{1-\lambda}$.

<!-- eq:14 -->
$$ \dot{\mathbf{p}}_i = \mathbf{f}_i + \mathbf{g}_i^{(GLE)} - \gamma \mathbf{p}_i + \mathbf{R}_i $$
- **what:** SGLD-GLE equation of motion (the one actually integrated); LD plus GLE guiding force.
- **symbols:** $\mathbf{g}_i^{(GLE)}$ - GLE guiding force (eq. 15); $\mathbf{R}_i$ - white random force (eq. 2).

<!-- eq:15 -->
$$ \mathbf{g}_i^{(GLE)} = \lambda \gamma\, \widetilde{\mathbf{p}}_i - \mu\, \widetilde{\mathbf{R}}_i $$
- **what:** SGLD-GLE guiding force: deterministic low-frequency momentum boost minus a random-force averaging term (this term restores detailed balance).
- **symbols:** $\widetilde{\mathbf{p}}_i$ - low-frequency momentum; $\widetilde{\mathbf{R}}_i$ - low-frequency (locally averaged) random force; $\mu$ - GLE parameter.

<!-- eq:16 -->
$$ \dot{\mathbf{g}}_i^{(GLE)} = -\frac{1}{t_L}\mathbf{g}_i^{(GLE)} + \frac{\lambda\gamma}{t_L}\mathbf{p}_i - \frac{\mu}{t_L}\mathbf{R}_i $$
- **what:** Markovian auxiliary ODE for the GLE guiding force; converts the non-Markovian GLE into a closed Markov system with $\mathbf{g}_i$ as auxiliary variable.
- **symbols:** as above; integrated alongside eq. (14).

<!-- eq:20 -->
$$ \rho(\omega) = \int_{-\infty}^{\infty} C(t)\, e^{-i\omega t}\, dt $$
- **what:** spectrum (Fourier transform) of the velocity autocorrelation function; used to diagnose which frequencies are enhanced.
- **symbols:** $C(t)$ - velocity autocorrelation function; $\omega$ - angular frequency; $\rho(\omega)$ - spectral density.

## Stationary distribution (for verification, not integration)

<!-- eq:18 -->
$$ \rho(\{\mathbf{r}_i\},\{\mathbf{p}_i\},\{\mathbf{g}_i\}) = C^{-1}\exp\left(-\frac{1}{kT}\left(E_p + \sum_i \left(\frac{|\mathbf{p}_i|^2}{2m_i} + \frac{t_L\, |\mathbf{g}_i|^2}{2 m_i \gamma \mu^2}\right)\right)\right) $$
- **what:** exact stationary density of the extended (position, momentum, guiding-force) system; guiding force is Gaussian-distributed and independent, so marginal is canonical.
- **symbols:** $C$ - normalization constant; $E_p$ - potential energy; $|\mathbf{g}_i|^2$ - squared guiding force. <!-- CHECK: raw eq.(18) printed the last numerator as t_L|g_i^2|^2; interpreted as t_L|g_i|^2 (a Gaussian in g_i with variance m_i*gamma*mu^2/t_L), consistent with Appendix eq.(A2) w(z) -->

<!-- eq:19 -->
$$ \overline{\rho}(\{\mathbf{r}_i\},\{\mathbf{p}_i\}) = \int \rho(\{\mathbf{r}_i\},\{\mathbf{p}_i\},\{\mathbf{g}_i\})\, d\{\mathbf{g}_i\} = \overline{C}^{-1}\exp\left(-\frac{1}{kT}\left(E_p + \sum_i \frac{|\mathbf{p}_i|^2}{2m_i}\right)\right) $$
- **what:** marginal distribution over positions and momenta = exact canonical (NVT) Maxwell-Boltzmann; proves SGLD-GLE samples the canonical ensemble exactly.
- **symbols:** $\overline{C}$ - marginal normalization constant.

## Derivations (not implemented)

The Fokker-Planck equation eq. (17) and the Appendix proof (eqs. A1-A9,
$w(z)$, $\sigma$, $\omega$, $\kappa$, compact GLE eq. A6, FPE eq. A7) establish
that eq. (18) is the exact stationary solution. These are analysis steps, not
implementation targets; see `paper.md`.
