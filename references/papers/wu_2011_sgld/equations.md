# Equations — Wu & Brooks 2011, SGLD canonical ensemble

Reduced/CHARMM conventions: `k` = Boltzmann constant, `T` in Kelvin, energies in kcal/mol, time in
ps, momentum `p_i = m_i * ṙ_i`. Tilde (`~`) marks a local/evolving average (low-frequency property);
`x - x̃` is the high-frequency part.

<!-- eq:1 -->
$$ \langle \mathbf{p} \rangle_{L} = \frac{1}{L} \sum_{i=n-L+1}^{n} \mathbf{p}_{i} = \frac{1}{t_{L}} \int_{t-t_{L}}^{t} \mathbf{p}(\tau)\, d\tau \approx \left( 1 - \frac{1}{L} \right) \tilde{\mathbf{p}}_{n-1} + \frac{1}{L} \mathbf{p}_{n} = \left( 1 - \frac{\delta t}{t_{L}} \right) \tilde{\mathbf{p}}(t - \delta t) + \frac{\delta t}{t_{L}} \mathbf{p}(t) = \tilde{\mathbf{p}} $$
- **what:** Local average = mean over the most recent L points; the practical form is the recursive evolving (exponential moving) average with weight `δt/t_L` on the current value.
- **symbols:** `⟨p⟩_L` / `p̃` - low-frequency (locally averaged) momentum; L - local averaging size (# points); `t_L = L·δt` - local average time (ps); `δt` - time step (ps); `p_n = p(t)` - instantaneous value.

<!-- eq:1b -->
$$ \frac{d\tilde{\mathbf{p}}(t)}{dt} = \frac{\mathbf{p}(t) - \tilde{\mathbf{p}}(t)}{t_L} $$
- **what:** Continuous-time (`δt -> 0`) limit of the evolving average: first-order relaxation of the average toward the instantaneous value with time constant `t_L`.
- **symbols:** as above.

<!-- eq:2 -->
$$ \tilde{\mathbf{p}}(t) = \frac{1}{t_L} \int_0^t \mathbf{p}(\tau)\, e^{-\frac{t-\tau}{t_L}}\, d\tau $$
- **what:** Solution of eq:1b — evolving average is an exponentially-decaying-weighted history of the property, decay rate `1/t_L`.
- **symbols:** as above.

<!-- eq:3 -->
$$ \tilde{q}(t) = \frac{2\pi \varpi t_L (e^{-t/t_L} - \cos(2\pi \varpi t)) + \sin(2\pi \varpi t)}{1 + 4\pi^2 t_L^2 \varpi^2} $$
- **what:** Evolving average of the test signal `q(t)=sin(2π ϖ t)`; used to show frequency response of the averaging (test fixture).
- **symbols:** ϖ - signal frequency; `q(t)=sin(2π ϖ t)` - test function; other symbols as above.

<!-- eq:4 -->
$$ q(t) - \tilde{q}(t) = \frac{-2\pi \varpi t_L (e^{-t/t_L} - \cos(2\pi \varpi t)) + 4\pi^2 t_L^2 \varpi^2 \sin(2\pi \varpi t)}{1 + 4\pi^2 t_L^2 \varpi^2} $$
- **what:** High-frequency portion of the test signal (instantaneous minus evolving average).
- **symbols:** as eq:3.

<!-- eq:5 -->
$$ \tilde{E}_k = \frac{1}{2} \sum_i \frac{\tilde{p}_i^2}{m_i} $$
- **what:** Low-frequency kinetic energy from low-frequency momenta.
- **symbols:** `Ẽ_k` - low-frequency kinetic energy (kcal/mol); `p̃_i` - low-frequency momentum of atom i; `m_i` - mass.

<!-- eq:6 -->
$$ \tilde{T} = \frac{\tilde{E}_k}{N_{\rm DF}\, k} $$
- **what:** Low-frequency temperature from low-frequency kinetic energy.
- **symbols:** `T̃` - low-frequency temperature (K); `N_DF` - number of degrees of freedom; `k` - Boltzmann constant.

<!-- eq:7 -->
$$ \dot{\mathbf{p}}_i = \mathbf{f}_i - \gamma_i \mathbf{p}_i + \mathbf{R}_i $$
- **what:** Langevin dynamics equation of motion.
- **symbols:** `ṗ_i` - time derivative of momentum; `f_i` - interaction (potential) force; `γ_i` - collision (friction) frequency (1/ps); `R_i` - random force.

<!-- eq:8 -->
$$ \langle \mathbf{R}_i(0)\mathbf{R}_i(t)\rangle = 2 m_i k T \gamma_i \delta(t) $$
- **what:** Fluctuation-dissipation relation defining the random-force autocorrelation.
- **symbols:** `δ(t)` - Dirac delta; T - target temperature (K); other symbols as eq:7.

<!-- eq:9 -->
$$ \dot{\mathbf{p}}_i = \mathbf{f}_i + \mathbf{g}_i - \gamma_i \mathbf{p}_i + \mathbf{R}_i $$
- **what:** SGLD equation of motion — Langevin plus guiding force `g_i`.
- **symbols:** `g_i` - guiding force (eq:10); others as eq:7.

<!-- eq:10 -->
$$ \mathbf{g}_{i}(t) = \lambda_{i} \gamma_{i} (\tilde{\mathbf{p}}_{i}(t) - \xi \mathbf{p}_{i}(t)) $$
- **what:** Guiding force from low-frequency momentum, with energy-conservation subtraction `ξ p_i`.
- **symbols:** `λ_i` - guiding factor (dimensionless, input); ξ - energy-conservation factor (eq:12); `p̃_i` - low-frequency momentum; `p_i` - instantaneous momentum.

<!-- eq:11 -->
$$ \sum_{i} \mathbf{g}_{i} \cdot \dot{\mathbf{r}}_{i} = \sum_{i} \lambda_{i} \gamma_{i} \tilde{\mathbf{p}}_{i} \cdot \dot{\mathbf{r}}_{i} - \xi \sum_{i} \lambda_{i} \gamma_{i} \mathbf{p}_{i} \cdot \dot{\mathbf{r}}_{i} = 0 $$
- **what:** Energy-conservation constraint: total guiding-force power is zero.
- **symbols:** `ṙ_i` - velocity of atom i; others as eq:10.

<!-- eq:12 -->
$$ \xi = \frac{\sum_{i} \lambda_{i} \gamma_{i} \tilde{\mathbf{p}}_{i} \cdot \dot{\mathbf{r}}_{i}}{\sum_{i} \lambda_{i} \gamma_{i} \mathbf{p}_{i} \cdot \dot{\mathbf{r}}_{i}} $$
- **what:** Energy-conservation factor solved from eq:11.
- **symbols:** as eq:11.

<!-- eq:13 -->
$$ \Theta_{\text{SGLD}} = \sum \Omega \exp\left(-\frac{\lambda_{\text{lf}}\tilde{E}_p}{kT_{\text{lf}}} - \frac{\lambda_{\text{hf}}(E_p - \tilde{E}_p)}{kT_{\text{hf}}}\right) $$
- **what:** SGLD partition function split into low- and high-frequency energy surfaces.
- **symbols:** `Θ_SGLD` - SGLD partition function; Ω - density of states; `λ_lf`,`λ_hf` - low/high-frequency energy factors; `Ẽ_p` - low-frequency potential energy; `E_p` - potential energy; `T_lf`,`T_hf` - effective temperatures in low/high-frequency space.

<!-- eq:15 -->
$$ \lambda_{\rm lf} = \frac{\left\langle \sum_{i} (\tilde{\mathbf{f}}_{i} + \tilde{\mathbf{g}}_{i} - \gamma_{i} \tilde{\mathbf{p}}_{i}) \cdot \tilde{\mathbf{f}}_{i} \right\rangle}{\left\langle \sum_{i} \tilde{\mathbf{f}}_{i} \cdot \tilde{\mathbf{f}}_{i} \right\rangle} $$
- **what:** Low-frequency energy factor = projection of the total low-frequency force onto the low-frequency force direction.
- **symbols:** `f̃_i` - low-frequency force; `g̃_i` - low-frequency guiding force; `p̃_i` - low-frequency momentum; `⟨·⟩` - ensemble/time average.

<!-- eq:17 -->
$$ \lambda_{\rm hf} = \frac{\left\langle \sum_{i} \left( \mathbf{f}_{i} - \tilde{\mathbf{f}}_{i} + \mathbf{g}_{i} - \tilde{\mathbf{g}}_{i} - \gamma_{i} (\mathbf{p}_{i} - \tilde{\mathbf{p}}_{i}) \right) \cdot (\mathbf{f}_{i} - \tilde{\mathbf{f}}_{i}) \right\rangle}{\left\langle \sum_{i} \left( \mathbf{f}_{i} - \tilde{\mathbf{f}}_{i} \right) \cdot (\mathbf{f}_{i} - \tilde{\mathbf{f}}_{i}) \right\rangle} $$
- **what:** High-frequency energy factor = projection of total high-frequency force onto the high-frequency force direction.
- **symbols:** high-frequency parts `x - x̃` of force, guiding force, momentum.

<!-- eq:Tlf_Thf -->
$$ T_{\rm lf} = C_{\rm lf}\tilde{T}, \qquad T_{\rm hf} = C_{\rm hf}(T - \tilde{T}), \qquad C_{\rm lf} = \frac{T}{\tilde{T}_0}, \qquad C_{\rm hf} = \frac{T}{T - \tilde{T}_0} $$
- **what:** Effective temperatures proportional to low/high-frequency temperatures; constants fixed by the `λ=0` reference where `T_lf=T_hf=T`.
- **symbols:** `T̃_0` - reference low-frequency temperature (`T̃` at `λ=0`); `C_lf`,`C_hf` - proportionality constants; `T̃` - low-frequency temperature.

<!-- eq:18 -->
$$ \Theta_{\text{SGLD}} = \sum \Omega \exp\left(-\frac{\lambda_{\text{lf}}\tilde{T}_0}{\tilde{T}}\frac{\tilde{E}_p}{kT} - \frac{\lambda_{\text{hf}}(T - \tilde{T}_0)}{T - \tilde{T}}\frac{E_p - \tilde{E}_p}{kT}\right) $$
- **what:** SGLD partition function expressed with temperature ratios (substituting the `T_lf`,`T_hf` forms into eq:13).
- **symbols:** as eq:13 plus `T̃_0`, `T̃`, T.

<!-- eq:19 -->
$$ \dot{\tilde{\mathbf{p}}}_i = \tilde{\mathbf{f}}_i - \chi_{\rm lf} \gamma_i \tilde{\mathbf{p}}_i + \tilde{\mathbf{R}}_i $$
- **what:** Low-frequency equation of motion rewritten as Langevin dynamics with effective collision frequency `χ_lf γ_i`.
- **symbols:** `χ_lf` - low-frequency collision factor (eq:20); `R̃_i` - low-frequency random force.

<!-- eq:20 -->
$$ \chi_{\text{lf}} = \frac{\sum_{i} (\gamma_{i} \tilde{\mathbf{p}}_{i} - \tilde{\mathbf{g}}_{i}) \cdot \gamma_{i} \tilde{\mathbf{p}}_{i}}{\sum_{i} \gamma_{i}^{2} \tilde{\mathbf{p}}_{i} \cdot \tilde{\mathbf{p}}_{i}} $$
- **what:** Low-frequency collision factor from projecting the effective low-frequency friction.
- **symbols:** as eq:19.

<!-- eq:22 -->
$$ \tilde{T}_0 = \tilde{T} \, \chi_{\rm lf} $$
- **what:** Reference low-frequency temperature estimated from the same SGLD run (fluctuation-dissipation, guiding force does not affect random force).
- **symbols:** as eq:19, eq:6.

<!-- eq:23 -->
$$ \Theta_{\text{SGLD}} \approx \sum \Omega \exp \left( -\lambda_{\text{lf}} \chi_{\text{lf}} \frac{\tilde{E}_p}{kT} - \lambda_{\text{hf}} \frac{T - \chi_{\text{lf}} \tilde{T}}{T - \tilde{T}} \frac{E_p - \tilde{E}_p}{kT} \right) $$
- **what:** SGLD partition function using `T̃_0 = T̃ χ_lf` (self-contained; no separate `λ=0` run needed).
- **symbols:** as eq:18 plus `χ_lf`.

<!-- eq:24 -->
$$ \Theta_{LD} = \sum \Omega \exp\left(-\frac{\tilde{E}_{p}}{kT} - \frac{E_{p} - \tilde{E}_{p}}{kT}\right) = \Theta_{\text{SGLD}} \langle w_{\text{SGLD}} \rangle_{\text{SGLD}} $$
- **what:** Canonical (LD) partition function equals the SGLD partition function times the average SGLD weighting factor.
- **symbols:** `Θ_LD` - LD/canonical partition function; `w_SGLD` - SGLD weighting factor (eq:25); `⟨·⟩_SGLD` - average over the SGLD ensemble.

<!-- eq:25 -->
$$ w_{\text{SGLD}} \approx \exp\left(\left(\lambda_{\text{lf}}\chi_{\text{lf}} - 1\right)\frac{\tilde{E}_{p}}{kT} + \left(\lambda_{\text{hf}}\frac{T - \chi_{\text{lf}}\tilde{T}}{T - \tilde{T}} - 1\right)\frac{E_{p} - \tilde{E}_{p}}{kT}\right) $$
- **what:** SGLD reweighting factor converting a SGLD sample to canonical (LD) weight. KEY implementable output. Exact form uses `λ_lf T̃_0/T̃` and `λ_hf (T-T̃_0)/(T-T̃)`; the approximation substitutes `T̃_0 = χ_lf T̃`.
- **symbols:** `λ_lf`,`λ_hf`,`χ_lf`,`Ẽ_p`,`E_p`,`T̃`,`T`,`k` as above.

<!-- eq:26 -->
$$ \langle A \rangle = \frac{\langle A\, w_{\text{SGLD}} \rangle_{\text{SGLD}}}{\langle w_{\text{SGLD}} \rangle_{\text{SGLD}}} $$
- **what:** Reweighted canonical ensemble average of any observable A from a SGLD trajectory.
- **symbols:** A - observable; `w_SGLD` - eq:25.

<!-- eq:27 -->
$$ T_{\rm sg} = \frac{T_{\rm lf}}{T_{\rm hf}} T = \frac{\tilde{T}(T - \tilde{T}_0)}{\tilde{T}_0(T - \tilde{T})} T $$
- **what:** Self-guiding temperature — a temperature-unit measure of conformational search ability. `T_sg=T` for LD.
- **symbols:** `T_sg` - self-guiding temperature (K); others as eqs 6, 22.

## Appendix — leap-frog Verlet SGLD algorithm

<!-- eq:A1 -->
$$ \rho(\mathbf{R}_i) = \frac{1}{\sqrt{4\pi \gamma_i m_i kT}}\, e^{-\frac{\mathbf{R}_i^2}{4\gamma_i m_i kT}} $$
- **what:** Gaussian distribution of the random force (zero mean). NOTE: variance implied here is `2 γ_i m_i kT` (per-step), consistent with eq:8 up to the `δt` discretization.
- **symbols:** `ρ` - probability density of `R_i`; other symbols as eq:8. <!-- CHECK: normalization written with 4π γ m kT; verify against per-step variance 2 m k T γ / δt used in the integrator -->

<!-- eq:A2 -->
$$ \tilde{\mathbf{p}}_{i}(t) = \left(1 - \frac{\delta t}{t_{L}}\right)\tilde{\mathbf{p}}_{i}(t - \delta t) + \frac{\delta t}{t_{L}}\mathbf{p}_{i}\left(t - \frac{\delta t}{2}\right) $$
- **what:** Evolving update of low-frequency momentum using the previous half-step momentum (leap-frog).
- **symbols:** as eq:1.

<!-- eq:A3 -->
$$ \dot{\mathbf{r}}_{i}(t) = \dot{\mathbf{r}}_{i}\!\left( t - \tfrac{\delta t}{2} \right) + \frac{\delta t}{2 m_{i}} (\mathbf{f}_{i}(t) + \mathbf{g}'_{i}(t) + \mathbf{R}_{i}(t)) - \frac{\delta t}{2} (\gamma_{i} + \xi \lambda_{i} \gamma_{i}) \dot{\mathbf{r}}_{i}(t) $$
- **what:** Half-step velocity including friction and energy-conservation friction.
- **symbols:** `g'_i = λ_i γ_i p̃_i` - uncorrected guiding force; others as before.

<!-- eq:A4 -->
$$ \dot{\mathbf{r}}'_{i}(t) = \dot{\mathbf{r}}_{i}\!\left( t - \tfrac{\delta t}{2} \right) + \frac{\delta t}{2 m_{i}} (\mathbf{f}_{i}(t) + \mathbf{g}'_{i}(t) + \mathbf{R}_{i}(t)) $$
- **what:** Friction-free half-step velocity (auxiliary quantity `ṙ'_i`).
- **symbols:** `ṙ'_i` - friction-free half-step velocity; others as eq:A3.

<!-- eq:A5 -->
$$ \dot{\mathbf{r}}_{i}(t) = \frac{\dot{\mathbf{r}}'_{i}(t)}{1 + \frac{(1 + \xi \lambda_{i})\gamma_{i}\delta t}{2}} \approx \frac{\dot{\mathbf{r}}'_{i}(t)}{1 + \frac{\gamma_{i}\delta t}{2}} - \frac{\dot{\mathbf{r}}'_{i}(t)}{\left(1 + \frac{\gamma_{i}\delta t}{2}\right)^{2}} \frac{\xi \lambda_{i}\gamma_{i}\delta t}{2} $$
- **what:** Velocity expressed from the friction-free velocity, linearized in ξ.
- **symbols:** as eq:A4.

<!-- eq:A6 -->
$$ \xi = \frac{\sum_{i}^{N} \lambda_{i} \gamma_{i} \tilde{\mathbf{p}}_{i}(t) \cdot \dot{\mathbf{r}}'_{i}(t) \left(1 + \frac{\gamma_{i} \delta t}{2}\right)^{-1}}{\sum_{i}^{N} \lambda_{i} \gamma_{i} m_{i} \dot{\mathbf{r}}'^{2}_{i}(t) \left(1 + \frac{\gamma_{i} \delta t}{2}\right)^{-2} + \frac{\delta t}{2} \sum_{i}^{N} \lambda_{i}^{2} \gamma_{i}^{2} \tilde{\mathbf{p}}_{i}(t) \cdot \dot{\mathbf{r}}'_{i}(t) \left(1 + \frac{\gamma_{i} \delta t}{2}\right)^{-2}} $$
- **what:** Discrete-integrator energy-conservation factor (neglecting higher powers of ξ).
- **symbols:** N - number of atoms; others as above.

<!-- eq:A7 -->
$$ \mathbf{g}_{i}(t) = \lambda_{i} \gamma_{i} \tilde{\mathbf{p}}_{i}(t) - \xi \mathbf{p}_{i}(t) = \lambda_{i} \gamma_{i} \tilde{\mathbf{p}}_{i}(t) - \frac{\xi m_{i} \dot{\mathbf{r}}'_{i}(t)}{1 + \frac{(1 + \xi \lambda_{i}) \gamma_{i} \delta t}{2}} $$
- **what:** Actual (energy-corrected) guiding force used in the integrator.
- **symbols:** as eq:10, eq:A4.

<!-- eq:A8 -->
$$ \lambda_{lf} = 1 + \frac{GLF}{FLF}, \qquad \lambda_{hf} = 1 + \frac{GHF}{FHF}, \qquad \chi_{lf} = \frac{\tilde{T}_0}{\tilde{T}} = 1 - \frac{GPLF}{PPLF} $$
- **what:** Collision/energy factors from the running accumulators.
- **symbols:** FLF, FHF, GLF, GHF, PPLF, GPLF - accumulators (see below); others as before.

Accumulators (summed over steps t and atoms i):
$$ FLF = \sum_t \sum_i \tilde{\mathbf{f}}_i \cdot \tilde{\mathbf{f}}_i, \quad FHF = \sum_t \sum_i (\mathbf{f}_i - \tilde{\mathbf{f}}_i)\cdot(\mathbf{f}_i - \tilde{\mathbf{f}}_i) $$
$$ GLF = \sum_t \sum_i (\tilde{\mathbf{g}}_i - \gamma_i \tilde{\mathbf{p}}_i)\cdot \tilde{\mathbf{f}}_i, \quad GHF = \sum_t \sum_i (\mathbf{g}_i - \tilde{\mathbf{g}}_i - \gamma_i(\mathbf{p}_i - \tilde{\mathbf{p}}_i))\cdot(\mathbf{f}_i - \tilde{\mathbf{f}}_i) $$
$$ PPLF = \sum_t \sum_i \gamma_i^2 \tilde{\mathbf{p}}_i \cdot \tilde{\mathbf{p}}_i, \quad GPLF = \sum_t \sum_i \tilde{\mathbf{g}}_i \cdot \gamma_i \tilde{\mathbf{p}}_i $$

<!-- eq:A9 -->
$$ w_{\text{SGLD}} = \exp\left(\left(\lambda_{\text{lf}}\chi_{\text{lf}} - 1\right) \frac{\tilde{E}_p - \bar{E}_p}{kT} + \left(\lambda_{\text{hf}} \frac{T - \chi_{\text{lf}}\tilde{T}}{T - \tilde{T}} - 1\right) \frac{E_p - \tilde{E}_p}{kT}\right) $$
- **what:** Practical weighting factor; the mean potential `Ē_p` is subtracted from the low-frequency energy to prevent exponential overflow.
- **symbols:** `Ē_p` - average potential energy over the trajectory; others as eq:25.

<!-- eq:A10 -->
$$ \dot{\mathbf{r}}_{i}\!\left(t + \tfrac{\delta t}{2}\right) = (2\chi_{i} - 1)\dot{\mathbf{r}}_{i}\!\left(t - \tfrac{\delta t}{2}\right) + \chi_{i}\frac{\delta t}{m_{i}}(\mathbf{f}_{i}(t) + \mathbf{g}_{i}(t) + \mathbf{R}_{i}(t)) $$
- **what:** Leap-frog velocity update to the next half step.
- **symbols:** `χ_i` - per-atom scaling parameter (eq:A11).

<!-- eq:A11 -->
$$ \chi_i = \left(1 + \frac{(1 + \xi \lambda_i)\gamma_i \delta t}{2}\right)^{-1} $$
- **what:** Velocity scaling parameter for the leap-frog update.
- **symbols:** as before.

<!-- eq:A12 -->
$$ \mathbf{r}_{i}(t+\delta t) = \mathbf{r}_{i}(t) + \dot{\mathbf{r}}_{i}\!\left(t+\tfrac{\delta t}{2}\right)\delta t $$
- **what:** Leap-frog position update.
- **symbols:** `r_i` - position of atom i.

<!-- eq:A13 -->
$$ \mathbf{f}_{i}^{\text{CON}}(t+\delta t) = \frac{2 m_{i}}{\delta t^{2}} (\mathbf{r}_{i}^{\text{CON}}(t+\delta t) - \mathbf{r}_{i}(t+\delta t)) $$
- **what:** Constraint force recovered from the SHAKE position correction; must be added to the low-frequency force.
- **symbols:** `r_i^CON` - constrained position; `f_i^CON` - constraint force.

## Derivations (not implemented)

The proportionality assumptions (`T_lf = C_lf T̃`, `T_hf = C_hf (T-T̃)`) and the algebra combining
eq:13 -> eq:18 -> eq:23 -> eq:24/25 are first-order-perturbation derivations; only the final
weighting factor (eq:25 / eq:A9), factor definitions (eq:15, eq:17, eq:20, eq:A8), and integrator
(eq:A1–A13) are implemented. Eq:3, eq:4 are test-signal derivations, kept for the frequency-response
check.
