# Equations - Nilmeier 2011 NCMC

> Eqs. 39, 41, 42 use the published Corrections (doi:10.1073/pnas.1207617109), not the
> original body. The merged source contained both; the corrected forms are authoritative.
> Proof-only steps (Eqs. 4-8, 35-38, 49-53) are summarized in `paper.md`; the implementable
> subset is extracted here.

<!-- eq:1 -->
$$ \pi_{\lambda}(x) = Z_{\lambda}^{-1} e^{-u_{\lambda}(x)}; \qquad Z_{\lambda} = \int_{\Gamma} dx \, e^{-u_{\lambda}(x)} $$
- **what:** Boltzmann probability of microstate x at thermodynamic state lambda, with partition function.
- **symbols:** x - microstate in phase space Gamma (coords, momenta, box dims); lambda - thermodynamic parameters; u_lambda - reduced (dimensionless) potential; Z_lambda - partition function (normalizer).

<!-- eq:2 -->
$$ u_{\lambda}(x) = \beta[H(x) + pV(x)] $$
- **what:** reduced potential in the isothermal-isobaric (NpT) ensemble.
- **symbols:** beta - inverse temperature 1/(k_B T); H(x) - Hamiltonian (may include biasing potential, invariant under momentum inversion); p - external pressure; V(x) - system volume.

<!-- eq:3 -->
$$ \pi(x,\lambda) = \frac{Z_{\lambda}\pi_{\lambda}(x)\omega_{\lambda}}{\sum_{\nu} \int_{\Gamma} dy \, Z_{\nu}\pi_{\nu}(y)\omega_{\nu}} $$
- **what:** expanded-ensemble joint distribution over (x, lambda) as a weighted mixture of thermodynamic states.
- **symbols:** omega_lambda > 0 - externally imposed weight of state lambda; sum over nu in the allowed set G (discrete or continuous). If G is a single lambda, pi(x,lambda) = pi_lambda(x).
<!-- CHECK: original denominator subscript printed as "\nu,\nu,\nu"; interpreted as a single sum over nu plus an integral over y. -->

<!-- eq:9 -->
$$ \frac{\alpha(\tilde{X}|\tilde{\Lambda})}{\alpha(X|\Lambda)} \equiv \prod_{t=1}^{T} \frac{\alpha_t(\tilde{x}_t^*, \tilde{x}_{t-1})}{\alpha_t(x_{t-1}, x_t^*)} $$
- **what:** ratio of reverse-to-forward perturbation-kernel path probabilities.
- **symbols:** alpha_t - perturbation kernel at step t; x_t^* - perturbed config before propagation; tilde denotes momentum inversion; T - number of protocol steps.

<!-- eq:10 -->
$$ e^{-\Delta \mathcal{S}(X|\Lambda)} \equiv \prod_{t=1}^{T} \frac{K_t(\tilde{x}_t, \tilde{x}_t^*)}{K_t(x_t^*, x_t)} $$
- **what:** ratio of reverse-to-forward propagation kernels defines the conditional path-action difference Delta S.
- **symbols:** K_t - propagation kernel at step t; x_t - config after propagation; Delta S(X|Lambda) - conditional path action difference along trajectory X under protocol Lambda.

<!-- eq:11 -->
$$ \frac{A(X|\Lambda)}{A(\tilde{X}|\tilde{\Lambda})} = \frac{\omega_{T}}{\omega_{0}} \frac{P(\tilde{\Lambda}|\tilde{x}_{T}, \lambda_{T})}{P(\Lambda|x_{0}, \lambda_{0})} \frac{\alpha(\tilde{X}|\tilde{\Lambda})}{\alpha(X|\Lambda)} e^{-\Delta \mathcal{S}(X|\Lambda) - \Delta u(X|\Lambda)} $$
- **what:** main result - ratio of forward/reverse NCMC acceptance probabilities required by pathwise detailed balance.
- **symbols:** A(X|Lambda) - acceptance probability of trajectory X; omega_0, omega_T - expanded-ensemble weights at start/end states; P(Lambda|...) - protocol selection probability; Delta u = u_T(x_T) - u_0(x_0) - reduced energy difference.
<!-- CHECK: source line 151 printed the leading factor as "omega_T / alpha_0"; from Eqs. 19-20 the correct factor is omega_T / omega_0. -->

<!-- eq:12 -->
$$ A(X|\Lambda) = \min \left\{ 1, \frac{\omega_T}{\omega_0} \frac{P(\tilde{\Lambda}|\tilde{x}_T, \lambda_T)}{P(\Lambda|x_0, \lambda_0)} \frac{\alpha(\tilde{X}|\tilde{\Lambda})}{\alpha(X|\Lambda)} e^{-\Delta \mathcal{S}(X|\Lambda) - \Delta u(X|\Lambda)} \right\} $$
- **what:** Metropolis-Hastings NCMC acceptance probability satisfying Eq. 11.
- **symbols:** as in Eq. 11. Draw uniform U in [0,1); accept candidate if A > U, else flip momentum and keep initial state.

<!-- eq:13 -->
$$ \eta(\phi'|\phi) = [2\pi I_0(\kappa)]^{-1} e^{\kappa \cos(\phi' - \phi)} $$
- **what:** von Mises circular proposal for a stochastically driven torsion angle.
- **symbols:** phi - current torsion; phi' - proposed torsion; I_0(kappa) - modified Bessel function order zero; kappa > 0 - dimensionless force constant (concentration).

<!-- eq:14 -->
$$ \frac{\alpha_t(\tilde{\mathbf{y}}, \tilde{\mathbf{x}})}{\alpha_t(\mathbf{x}, \mathbf{y})} = \frac{\eta(\phi|\phi')J(\phi)}{\eta(\phi'|\phi)J(\phi')} = 1 $$
- **what:** perturbation-kernel ratio for torsion rotation equals 1 because rotation about a bond preserves Cartesian phase-space volume.
- **symbols:** J(phi) - Jacobian from torsion to Cartesian coords; here J(phi')=J(phi)=1.

<!-- eq:15 -->
$$ \pi_t(x)K_t(x,y) = \pi_t(\tilde{y})K_t(\tilde{y},\tilde{x}) $$
- **what:** detailed-balance condition satisfied by a reversible MCMC propagation kernel at step t.
- **symbols:** pi_t(x) = Z_t^{-1} e^{-u_t(x)} - instantaneous target at step t; K_t - propagation kernel; tilde - momentum inversion.

<!-- eq:16 -->
$$ w(X|\Lambda) = \sum_{t=1}^{T} [u_t(x_t^*) - u_{t-1}(x_{t-1})] $$
- **what:** nonequilibrium work accumulated by the perturbation steps.
- **symbols:** u_t - reduced potential at step t; x_t^* - perturbed config; x_{t-1} - post-propagation config of previous step.

<!-- eq:17 -->
$$ q(X|\Lambda) = \sum_{t=1}^{T} [u_t(x_t) - u_t(x_t^*)] $$
- **what:** heat accumulated by the propagation steps. Satisfies w + q = Delta u (first law).
- **symbols:** u_t(x_t) - reduced potential after propagation; u_t(x_t^*) - after perturbation.

<!-- eq:18 -->
$$ \Delta \mathcal{S}(X|\Lambda) = -\ln \prod_{t=1}^{T} \frac{\pi_t(x_t^*)}{\pi_t(x_t)} = -q(X|\Lambda) $$
- **what:** for detailed-balance propagation kernels the path-action difference equals minus the heat.
- **symbols:** pi_t - instantaneous target; q - heat (Eq. 17).

<!-- eq:19 -->
$$ \frac{A(X|\Lambda)}{A(\tilde{X}|\tilde{\Lambda})} = \frac{\omega_T}{\omega_0} \frac{P(\tilde{\Lambda}|\tilde{x}_T, \lambda_T)}{P(\Lambda|x_0, \lambda_0)} \frac{\alpha(\tilde{X}|\tilde{\Lambda})}{\alpha(X|\Lambda)} e^{-w(X|\Lambda)} $$
- **what:** acceptance ratio with reversible-MCMC propagation - the work w replaces the instantaneous energy difference.
- **symbols:** w(X|Lambda) - nonequilibrium work (Eq. 16); other symbols as Eq. 11.

<!-- eq:20 -->
$$ \frac{A(X|\Lambda)}{A(\tilde{X}|\tilde{\Lambda})} = \frac{\omega_T}{\omega_0} \frac{P(\tilde{\Lambda}|\tilde{x}_T, \lambda_T)}{P(\Lambda|x_0, \lambda_0)} \frac{\alpha(\tilde{X}|\tilde{\Lambda})}{\alpha(X|\Lambda)} e^{-\Delta u(X|\Lambda)} $$
- **what:** acceptance ratio for symplectic (deterministic, volume-preserving) propagation where Delta S = 0; work equals energy difference.
- **symbols:** Delta u = u_T(x_T) - u_0(x_0) - reduced energy difference; other symbols as Eq. 11.

<!-- eq:21 -->
$$ C(t) = C_{\infty} + (C_0 - C_{\infty})e^{-t/\tau} $$
- **what:** single-exponential model of the dimer-extension autocorrelation function used to extract the correlation time.
- **symbols:** C(t) = <r(0)r(t)>; C_0 = <r^2>; C_inf = <r>^2; tau - integrated correlation time (iterations).

<!-- eq:g -->
$$ g = 1 + 2\tau $$
- **what:** statistical inefficiency - number of iterations per effectively uncorrelated sample.
- **symbols:** g - statistical inefficiency; tau - correlation time (iterations).

<!-- eq:22 -->
$$ \tau_{\rm eff} = \tau_{\rm MD} \left[ \frac{\tau_{\rm NCMC}}{\tau_{\rm MD} + \tau_{\rm NCMC}} \right] $$
- **what:** effective correlation time for mixed MD+NCMC iterations (equivalent to harmonic-style combination, Eq. 54).
- **symbols:** tau_MD - correlation time for MD-only iterations; tau_NCMC - correlation time for NCMC-only iterations.

<!-- eq:tau_ncmc -->
$$ \tau_{\rm NCMC} \approx -\frac{1}{\ln(1 - 2\gamma)} $$
- **what:** NCMC correlation time estimated from the average acceptance probability.
- **symbols:** gamma - average NCMC acceptance probability (assumed symmetric between the two states).

<!-- eq:23 -->
$$ E \equiv \frac{g_{\rm MD} T_{\rm MD}}{g_{\rm NCMC} (T_{\rm MD} + T_{\rm NCMC})} $$
- **what:** statistical efficiency gain of MD+NCMC relative to MD alone at fixed compute (force evaluations).
- **symbols:** g_MD, g_NCMC - statistical inefficiencies; T_MD, T_NCMC - force evaluations per iteration for the MD and NCMC parts.

<!-- eq:24 -->
$$ U_{\text{bond}}(r) = h \left[ 1 - \frac{(r - r_0 - s)^2}{s^2} \right]^2 $$
- **what:** double-well bonded potential of the bistable dimer (minima at r0 and 2r0).
- **symbols:** r - interparticle distance; h - barrier height (= 5 k_B T); r_0 = r_WCA - compact minimum; s = r_WCA/2 - well-width parameter.

<!-- eq:25 -->
$$ U_{\text{WCA}}(r) = \begin{cases} 4\epsilon \left[ \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^{6} \right] + \epsilon, & r < r_{\text{WCA}} \\ 0, & r \ge r_{\text{WCA}} \end{cases} $$
- **what:** Weeks-Chandler-Andersen purely repulsive nonbonded potential (LJ shifted and truncated at its minimum).
- **symbols:** epsilon - LJ well depth; sigma - LJ diameter; r_WCA = 2^(1/6) sigma - cutoff at LJ minimum.

<!-- eq:26 -->
$$ \Delta r = \begin{cases} +r_0 & \text{if } r < 1.5r_0 \\ -r_0 & \text{if } 1.5r_0 \le r \le 3r_0 \\ 0 & \text{otherwise} \end{cases} $$
- **what:** deterministic extension/contraction rule for the instantaneous MC dimer move.
- **symbols:** r - current dimer extension; r_0 - compact-minimum distance; Delta r - proposed change to extension.

<!-- eq:27 -->
$$ A(x_{\text{new}}|x_{\text{old}}) = \min\{1, e^{-\beta[U(x_{\text{new}}) - U(x_{\text{old}})]}J_r(x_{\text{old}}, x_{\text{new}})\} $$
- **what:** Metropolis-Hastings acceptance for the instantaneous MC dimer move with a radial Jacobian.
- **symbols:** J_r(x_old,x_new) = (r_new/r_old)^2 - Jacobian for radial expansion/contraction; beta - inverse temperature; U - total potential energy.

<!-- eq:28 -->
$$ A(X) = \min\{1, e^{-\beta[H(x_T) - H(x_0)]}J_r(x_0, x_T)\} $$
- **what:** NCMC acceptance for the T-step symplectic dimer move (Hamiltonian difference plus radial Jacobian).
- **symbols:** H - Hamiltonian; J_r(x_0,x_T) = (r_new/r_old)^2; x_0, x_T - start/end microstates of the switching trajectory.

<!-- eq:29 -->
$$ \alpha_t(x,y) = \left[\frac{r(y)}{r(x)}\right]^2 \delta([r(y) - r(x)] - [\Delta r/T]) $$
- **what:** perturbation kernel driving the dimer separation by a fixed increment Delta r / T each step.
- **symbols:** r(x) - dimer separation of config x; Delta r - total extension change; T - number of switching steps; delta - Dirac delta.

<!-- eq:30 -->
$$ \langle A \rangle_{\tau} \approx \frac{1}{N} \sum_{n=1}^{N} A(X_n) $$
- **what:** sample-mean estimator of the mean acceptance probability at switching time tau.
- **symbols:** N - number of trials; A(X_n) - acceptance probability of trial n; tau here indexes the switching length.

<!-- eq:31 -->
$$ \ln \langle A \rangle_{\tau} \approx -\ln N + b + \ln \sum_{n=1}^{N} e^{a_n - b} $$
- **what:** numerically stable log-sum-exp estimator of the log mean acceptance from stored log-acceptances.
- **symbols:** a_n = ln A(X_n); b = max_n a_n - shift for stability; N - number of trials.
<!-- CHECK: source line 342 printed "ln b + sum e^{a_n-b}" without the outer ln on the sum; the standard log-sum-exp form is ln(sum) + b, transcribed here. -->

<!-- eq:32 -->
$$ \mathscr{P}_{\text{vac}}(r) \propto 4\pi r^2 e^{-\beta U_{\text{bond}}(r)} $$
- **what:** analytic reference radial distribution of dimer extension in vacuum (spherical Jacobian x Boltzmann).
- **symbols:** r - dimer extension; U_bond - Eq. 24; beta - inverse temperature.

<!-- eq:33 -->
$$ U_{\rm umbrella}(r) = k_{\rm B} T \ln r^2 + \theta(r_{\rm min} - r)(K/2)[r - r_{\rm min}]^2 + \theta(r - r_{\rm max})(K/2)[r - r_{\rm max}]^2 $$
- **what:** flat-bottomed umbrella potential removing the barrier to estimate the solvated distribution.
- **symbols:** theta - Heaviside step (1 for arg >= 0); r_min = r_0; r_max = 2.05 r_0; K = k_B T / eta^2 with eta = 0.3 Angstrom; the k_B T ln r^2 term cancels the radial Jacobian.

<!-- eq:34 -->
$$ \mathcal{P}_{\text{sol}}(r) \propto \frac{\sum_{n=1}^{N} \delta(r - r_n) e^{-\beta [U_{\text{bond}}(r_n) - U_{\text{umbrella}}(r_n)]}}{\sum_{n=1}^{N} e^{-\beta [U_{\text{bond}}(r_n) - U_{\text{umbrella}}(r_n)]}} $$
- **what:** reweighting estimator of the true solvated distribution from umbrella-sampled data.
- **symbols:** r_n - bond separation of sample n; a finite-width histogram bin replaces delta(r) in practice.

<!-- eq:39 -->
$$ x_t = x_t^* + \frac{\Delta t}{\gamma m} F_t(x_t^*) + \sqrt{2} \left(\frac{\Delta t}{\gamma m}\right)^{1/2} \xi_t $$
- **what:** corrected Ermak-Yeh overdamped-Langevin (Brownian) propagation step (coordinates only).
- **symbols:** x_t^* - perturbed coord; F_t(x) = -dH_t/dx - systematic force; gamma - collision frequency (inverse time); m - mass; Delta t - timestep; xi_t - Gaussian noise, mean 0, variance beta^-1.
<!-- CHECK: force sign is "+" per the published Corrections (Eq. 39); the original body printed "-". -->

<!-- eq:40 -->
$$ \phi(\xi_t) = \frac{1}{\sqrt{2\pi\beta^{-1}}} \exp\left[-\frac{\beta}{2}\xi_t^2\right] $$
- **what:** distribution of the Ermak-Yeh noise variate (zero mean, variance beta^-1).
- **symbols:** xi_t - per-DOF noise; beta - inverse temperature.

<!-- eq:41 -->
$$ x_{t}^{*} = x_{t} + \frac{\Delta t}{\gamma m} F_{t}(x_{t}) + \sqrt{2} \left(\frac{\Delta t}{\gamma m}\right)^{1/2} \tilde{\xi}_{t} $$
- **what:** corrected reverse Ermak-Yeh step (from x_t back to x_t^*) defining the reverse noise.
- **symbols:** tilde xi_t - reverse noise generating the reverse transition; other symbols as Eq. 39.
<!-- CHECK: force sign "+" per published Corrections (Eq. 41). -->

<!-- eq:42 -->
$$ \tilde{\xi}_{t} = -\frac{1}{\sqrt{2}} \left( \frac{\Delta t}{\gamma m} \right)^{1/2} [F_{t}(x_{t}) + F_{t}(x_{t}^{*})] - \xi_{t} $$
- **what:** corrected reverse-noise history for the Ermak-Yeh kernel in terms of forward noise and forces.
- **symbols:** tilde xi_t - reverse noise; xi_t - forward noise; F_t - systematic force at x_t and x_t^*.
<!-- CHECK: leading "-1/sqrt(2)" and trailing "- xi_t" per published Corrections (Eq. 42). -->

<!-- eq:43 -->
$$ \Delta\mathcal{S}(X) = -\ln\prod_{t=1}^T \frac{K_t(x_t, x_t^*)}{K_t(x_t^*, x_t)} = \frac{\beta}{2}\sum_{t=1}^T (\tilde{\xi}_t^2 - \xi_t^2) $$
- **what:** conditional path-action difference for the Ermak-Yeh integrator; Jacobians cancel, leaving squared-noise difference.
- **symbols:** tilde xi_t - reverse noise (Eq. 42); xi_t - forward noise; beta - inverse temperature.

<!-- eq:44 -->
$$ \begin{aligned} v_t' &= v_t^* + \frac{\Delta t}{2m} \left( F_t(r_t^*) - \gamma m v_t^* + \sqrt{\frac{2\gamma m}{\Delta t}} \xi_t \right) \\ r_t &= r_t^* + \Delta t\, v_t' \\ v_t &= \frac{1}{1 + \frac{\gamma \Delta t}{2}} \left[ v_t' + \frac{\Delta t}{2m} \left( F_t(r_t) + \sqrt{\frac{2\gamma m}{\Delta t}} \xi_t' \right) \right] \end{aligned} $$
- **what:** velocity-Verlet discretization of the Brünger-Brooks-Karplus (BBK) Langevin integrator, one step.
- **symbols:** r_t, v_t - Cartesian position/velocity of microstate x_t; v_t' - auxiliary velocity; gamma - collision frequency; m - mass; xi_t, xi_t' - two Gaussian noise variates per DOF per step (mean 0, variance beta^-1).

<!-- eq:45 -->
$$ \psi(\xi_t, \xi_t') = \frac{1}{2\pi\beta^{-1}} \exp\left[-\frac{\beta}{2}(\xi_t^2 + \xi_t'^2)\right] $$
- **what:** joint distribution of the two BBK noise variates.
- **symbols:** xi_t, xi_t' - per-DOF noise; beta - inverse temperature.

<!-- eq:46 -->
$$ \tilde{\xi}_t = \xi_t' - \sqrt{2\gamma m \Delta t}\; v_t, \qquad \tilde{\xi}_t' = \xi_t - \sqrt{2\gamma m \Delta t}\; v_t^* $$
- **what:** reverse BBK noise variables generating the time-reversed step (r_t,-v_t) -> (r_t^*,-v_t^*).
- **symbols:** v_t, v_t^* - post/pre-step velocities; gamma - collision frequency; m - mass; Delta t - timestep.
<!-- CHECK: source printed the sqrt as "sqrt(2 gamma m Delta t v_t)"; grouped here as sqrt(2 gamma m Delta t) * v_t for dimensional consistency with the noise variance. -->

<!-- eq:48 -->
$$ \Delta S(X) = -\ln \prod_{t=1}^{T} \frac{K_{t}(\tilde{x}_{t}, \tilde{x}_{t}^{*})}{K_{t}(x_{t}^{*}, x_{t})} = \frac{\beta}{2} \sum_{t=1}^{T} [(\tilde{\xi}_{t}^{2} + \tilde{\xi}_{t}^{\prime 2}) - (\xi_{t}^{2} + \xi_{t}^{\prime 2})] $$
- **what:** conditional path-action difference for the BBK integrator; Jacobian ratio cancels (noise-independent).
- **symbols:** tilde xi_t, tilde xi_t' - reverse noise (Eq. 46); xi_t, xi_t' - forward noise; beta - inverse temperature.
