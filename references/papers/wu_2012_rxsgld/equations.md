# Equations - Wu 2012 RXSGLD

Reduced/physical convention: `k` is Boltzmann's constant, `T` temperature,
`β = 1/(kT)`. Bold symbols are 3-vectors per particle. A tilde `~` denotes the
low-frequency (locally time-averaged) part of a quantity.

<!-- eq:1 -->
$$ \dot{\mathbf{p}}_i = \mathbf{f}_i + \mathbf{g}_i - \gamma_i \mathbf{p}_i + \mathbf{R}_i $$
- **what:** SGLD equation of motion for particle i (Langevin + guiding force).
- **symbols:** $\dot{\mathbf{p}}_i$ - time derivative of momentum ($R^3$); $\mathbf{f}_i$ - interaction force; $\mathbf{g}_i$ - guiding force (eq:3); $\gamma_i$ - collision frequency (1/time); $\mathbf{p}_i$ - momentum; $\mathbf{R}_i$ - random force (eq:2).

<!-- eq:2 -->
$$ \langle \mathbf{R}_i(0)\mathbf{R}_i(t)\rangle = 2 m_i k T \gamma_i \delta(t) $$
- **what:** fluctuation-dissipation relation for the Langevin random force.
- **symbols:** $m_i$ - particle mass; $k$ - Boltzmann constant; $T$ - simulation temperature; $\gamma_i$ - collision frequency; $\delta(t)$ - Dirac delta.

<!-- eq:3 -->
$$ \mathbf{g}_i(t) = \lambda_i \gamma_i \left( \tilde{\mathbf{p}}_i(t) - \xi \mathbf{p}_i(t) \right) $$
- **what:** guiding force from low-frequency momentum minus a conserving fraction of the momentum.
- **symbols:** $\lambda_i$ - guiding factor (strength; $\lambda_i=0$ recovers plain Langevin); $\tilde{\mathbf{p}}_i$ - low-frequency momentum (eq:6 applied to $\mathbf{p}$); $\xi$ - energy conservation factor (eq:5); $\gamma_i$ - collision frequency.

<!-- eq:4 -->
$$ \sum_i \mathbf{g}_i \cdot \dot{\mathbf{r}}_i = \sum_i \lambda_i \gamma_i \tilde{\mathbf{p}}_i \cdot \dot{\mathbf{r}}_i - \xi \sum_i \lambda_i \gamma_i \mathbf{p}_i \cdot \dot{\mathbf{r}}_i = 0 $$
- **what:** energy-conservation constraint: guiding force does zero net work over the whole system.
- **symbols:** $\dot{\mathbf{r}}_i$ - velocity of particle i; sum runs over all particles.

<!-- eq:5 -->
$$ \xi = \frac{\sum_i \lambda_i \gamma_i \tilde{\mathbf{p}}_i \cdot \dot{\mathbf{r}}_i}{\sum_i \lambda_i \gamma_i \mathbf{p}_i \cdot \dot{\mathbf{r}}_i} $$
- **what:** energy conservation factor determined each time step (solves eq:4).
- **symbols:** as above; ratio evaluated at the current step.

<!-- eq:6 -->
$$ \tilde{P}(t) = \left(1 - \frac{\delta t}{t_L}\right)\tilde{P}(t - \delta t) + \frac{\delta t}{t_L}P(t) $$
- **what:** progressive local (exponential moving) average defining the low-frequency part of any property P; a low-pass filter.
- **symbols:** $\tilde{P}$ - low-frequency part of P; $P$ - instantaneous value; $\delta t$ - time step; $t_L$ - local averaging time (e.g. 0.2 ps). High-frequency part is $P - \tilde{P}$.

<!-- eq:7 -->
$$ \Theta_{\text{SGLD}} \approx \sum \exp\left(-\frac{\lambda_{\text{lf}}\chi_{\text{lf}}\tilde{E}_{\text{p}}}{kT} - \frac{\lambda_{\text{hf}}\chi_{\text{hf}}(E_{\text{p}} - \tilde{E}_{\text{p}})}{kT}\right) $$
- **what:** configurational partition function of the SGLD ensemble; sum over all microscopic states.
- **symbols:** $E_p$ - potential energy; $\tilde{E}_p$ - low-frequency potential energy; $\lambda_{lf},\lambda_{hf}$ - low/high frequency energy factors (eq:8, eq:9); $\chi_{lf},\chi_{hf}$ - low/high frequency collision factors (eq:10, eq:11).

<!-- eq:8 -->
$$ \lambda_{\rm lf} = \frac{\left\langle \sum_i (\tilde{\mathbf{f}}_i + \tilde{\mathbf{g}}_i - \gamma_i \tilde{\mathbf{p}}_i) \tilde{\mathbf{f}}_i \right\rangle}{\left\langle \sum_i \tilde{\mathbf{f}}_i \tilde{\mathbf{f}}_i \right\rangle} $$
- **what:** low-frequency energy factor; average projection of total low-frequency force onto interaction force.
- **symbols:** $\tilde{\mathbf{f}}_i$ - low-frequency interaction force; $\tilde{\mathbf{g}}_i$ - low-frequency guiding force; $\langle\cdot\rangle$ - ensemble average. LD value = 1.

<!-- eq:9 -->
$$ \lambda_{\rm hf} = \frac{\left\langle \sum_i \left( \mathbf{f}_i - \tilde{\mathbf{f}}_i + \mathbf{g}_i - \tilde{\mathbf{g}}_i - \gamma_i (\mathbf{p}_i - \tilde{\mathbf{p}}_i) \right) (\mathbf{f}_i - \tilde{\mathbf{f}}_i) \right\rangle}{\left\langle \sum_i \left( \mathbf{f}_i - \tilde{\mathbf{f}}_i \right) (\mathbf{f}_i - \tilde{\mathbf{f}}_i) \right\rangle} $$
- **what:** high-frequency energy factor; projection of total high-frequency force onto high-frequency interaction force.
- **symbols:** $\mathbf{f}_i - \tilde{\mathbf{f}}_i$ - high-frequency interaction force; other terms high-frequency guiding/friction. LD value = 1.

<!-- eq:10 -->
$$ \chi_{\text{lf}} = \frac{\tilde{T}_0}{\tilde{T}} = 1 - \frac{\left\langle \sum_i \tilde{\mathbf{g}}_i \gamma_i \tilde{\mathbf{p}}_i \right\rangle}{\left\langle \sum_i \gamma_i^2 \tilde{\mathbf{p}}_i \tilde{\mathbf{p}}_i \right\rangle} $$
- **what:** low-frequency collision factor; ratio of reference to actual low-frequency temperature.
- **symbols:** $\tilde{T}$ - low-frequency temperature (eq:12); $\tilde{T}_0$ - reference low-frequency temperature (value at zero guiding). LD value = 1.

<!-- eq:11 -->
$$ \chi_{\rm hf} = \frac{T - \tilde{T}_0}{T - \tilde{T}} = \frac{T - \chi_{\rm lf} \tilde{T}}{T - \tilde{T}} = 1 - \frac{\left\langle \sum_i \gamma_i (\mathbf{g}_i - \tilde{\mathbf{g}}_i) \cdot (\mathbf{p}_i - \tilde{\mathbf{p}}_i) \right\rangle}{\left\langle \sum_i \gamma_i^2 (\mathbf{p}_i - \tilde{\mathbf{p}}_i) \cdot (\mathbf{p}_i - \tilde{\mathbf{p}}_i) \right\rangle} $$
- **what:** high-frequency collision factor. LD value = 1.
- **symbols:** $T$ - target temperature; $\tilde{T}$ - low-frequency temperature; $\tilde{T}_0$ - reference low-frequency temperature; $\chi_{lf}$ from eq:10.

<!-- eq:12 -->
$$ \tilde{T} = \frac{1}{N_{\rm DF}k} \left\langle \sum_i \frac{\tilde{\mathbf{p}}_i^2}{m_i} \right\rangle $$
- **what:** low-frequency temperature from the low-frequency momentum.
- **symbols:** $N_{DF}$ - number of degrees of freedom; $\tilde{\mathbf{p}}_i$ - low-frequency momentum; $m_i$ - mass; $k$ - Boltzmann constant.

<!-- eq:13 -->
$$ \Theta_{\rm LD} = \sum \exp\left(-\frac{E_{\rm p}}{kT}\right) $$
- **what:** canonical (Langevin) partition function; the SGLD partition function (eq:7) with all four factors = 1.
- **symbols:** $E_p$ - potential energy; sum over all microstates.

<!-- eq:14 -->
$$ w_{\text{SGLD}} = \exp\left((\lambda_{\text{lf}}\chi_{\text{lf}} - 1)\frac{\tilde{E}_{\text{p}}}{kT} + (\lambda_{\text{hf}}\chi_{\text{hf}} - 1)\frac{E_{\text{p}} - \tilde{E}_{\text{p}}}{kT}\right) $$
- **what:** SGLD reweighting factor connecting SGLD ensemble to canonical ensemble ($\Theta_{LD} = \Theta_{SGLD}\langle w_{SGLD}\rangle_{SGLD}$).
- **symbols:** $\tilde{E}_p$ - low-frequency potential energy; $E_p - \tilde{E}_p$ - high-frequency potential energy; factors from eq:8-11.

<!-- eq:15 -->
$$ \langle P \rangle_{\rm LD} = \frac{\langle P w_{\rm SGLD} \rangle_{\rm SGLD}}{\langle w_{\rm SGLD} \rangle_{\rm SGLD}} $$
- **what:** canonical (LD) ensemble average of any property P recovered by reweighting an SGLD run.
- **symbols:** $P$ - any observable; $w_{SGLD}$ - reweighting factor (eq:14); $\langle\cdot\rangle_{SGLD}$ - average over SGLD trajectory.

<!-- eq:16 -->
$$ T_{\rm SG} = \frac{\chi_{\rm hf}}{\chi_{\rm lf}} T = \frac{\tilde{T}(T - \tilde{T}_0)}{\tilde{T}_0(T - \tilde{T})} T $$
- **what:** self-guiding temperature; a measure of conformational searching ability in temperature units.
- **symbols:** $\chi_{lf},\chi_{hf}$ - collision factors (eq:10, eq:11); $T$ - temperature; $\tilde{T},\tilde{T}_0$ - low-frequency and reference low-frequency temperatures.

<!-- eq:17a -->
$$ T^{(i)} = T^{(0)} \left(\frac{T^{(k)}}{T^{(0)}}\right)^{\frac{i}{k}} $$
- **what:** exponential temperature ladder for stage i (of k+1 stages).
- **symbols:** $T^{(i)}$ - temperature at stage i; $T^{(0)}$ - base stage temperature; $T^{(k)}$ - top stage temperature; $k$ - number of stages above base. <!-- CHECK: raw printed the exponent as 1/k; the standard geometric ladder uses i/k so that stage 0 -> T^(0) and stage k -> T^(k). -->

<!-- eq:17b -->
$$ T_{\rm SG}^{(i)} = T_{\rm SG}^{(0)} \left(\frac{T_{\rm SG}^{(k)}}{T_{\rm SG}^{(0)}}\right)^{\frac{i}{k}} = T^{(0)} \left(\frac{T_{\rm SG}^{(k)}}{T^{(0)}}\right)^{\frac{i}{k}} $$
- **what:** exponential self-guiding-temperature ladder for stage i; right form uses $T_{SG}^{(0)} = T^{(0)}$.
- **symbols:** $T_{SG}^{(i)}$ - self-guiding temperature at stage i; $T_{SG}^{(0)}=T^{(0)}$ - base; $T_{SG}^{(k)}$ - top. <!-- CHECK: exponent i/k (see eq:17a note). -->

<!-- eq:18 -->
$$ \rho_{\text{SGLD}}(X_m^{(i)}) = \frac{1}{\Theta_{\text{SGLD}}^{(m)}} \exp\left(-\frac{\lambda_{\text{lf}}^{(m)} \chi_{\text{lf}}^{(m)} \tilde{E}_{p}^{(i)}}{k T_m} - \frac{\lambda_{\text{hf}}^{(m)} \chi_{\text{hf}}^{(m)} (E_{p}^{(i)} - \tilde{E}_{p}^{(i)})}{k T_m}\right) = \frac{1}{\Theta_{\text{SGLD}}^{(m)}} \exp\left(-\tilde{\mu}_m \tilde{E}_{p}^{(i)} - \mu_m E_{p}^{(i)}\right) $$
- **what:** distribution probability of conformation i at stage m under the SGLD ensemble.
- **symbols:** $X_m^{(i)}$ - conformation i at stage m; $T_m$ - stage-m temperature; superscript (m) - stage-m factors; $\tilde{\mu}_m$ (eq:19a), $\mu_m$ (eq:19b); $\Theta_{SGLD}^{(m)}$ - stage-m partition function.

<!-- eq:19a -->
$$ \tilde{\mu}_m = \frac{\lambda_{\text{lf}}^{(m)} \chi_{\text{lf}}^{(m)} - \lambda_{\text{hf}}^{(m)} \chi_{\text{hf}}^{(m)}}{kT_m} $$
- **what:** low-frequency exchange coefficient. Equivalently $\tilde{\mu}_m = \beta_m(\lambda_{lf}^{(m)}\chi_{lf}^{(m)} - \lambda_{hf}^{(m)}\chi_{hf}^{(m)})$.
- **symbols:** $\beta_m = 1/(kT_m)$; stage-m factors. <!-- CHECK: raw eq:19a wrote kT in denominator; consistent with eq:18/19b it should be kT_m. -->

<!-- eq:19b -->
$$ \mu_m = \frac{\lambda_{\text{hf}}^{(m)} \chi_{\text{hf}}^{(m)}}{kT_m} $$
- **what:** high-frequency exchange coefficient. Equivalently $\mu_m = \beta_m \lambda_{hf}^{(m)}\chi_{hf}^{(m)}$.
- **symbols:** $\beta_m = 1/(kT_m)$; stage-m factors.

<!-- eq:20 -->
$$ \pi_{\text{RX}}\left(\left\{\mathbf{X}_m^{[i]}, \mathbf{X}_n^{[j]}\right\} \to \left\{\mathbf{X}_m^{[j]}, \mathbf{X}_n^{[i]}\right\}\right) \approx \exp\left(-\left(\tilde{\mu}_m - \tilde{\mu}_n\right)\left(\tilde{E}_p\left(\mathbf{X}_n^{[j]}\right) - \tilde{E}_p\left(\mathbf{X}_m^{[i]}\right)\right) - \left(\mu_m - \mu_n\right)\left(E_p\left(\mathbf{X}_n^{[j]}\right) - E_p\left(\mathbf{X}_m^{[i]}\right)\right)\right) $$
- **what:** RXSGLD exchange probability between stages m and n (approximated with $T_m=T_n$ so low-frequency energies match per conformation). Accept via Metropolis $\min\{1,\pi_{RX}\}$.
- **symbols:** $\mathbf{X}_m^{[i]}$ - replica i at stage m; $\tilde{\mu},\mu$ - exchange coefficients (eq:19a,b); $\tilde{E}_p,E_p$ - low-frequency and total potential energies.

<!-- eq:21 -->
$$ \pi_{\text{TRXLD}}\left(\left\{\mathbf{X}_m^{[i]}, \mathbf{X}_n^{[j]}\right\} \to \left\{\mathbf{X}_m^{[j]}, \mathbf{X}_n^{[i]}\right\}\right) = \exp\left(-\left(\beta_m - \beta_n\right)\left(E_p\left(\mathbf{X}_n^{[j]}\right) - E_p\left(\mathbf{X}_m^{[i]}\right)\right)\right) $$
- **what:** standard temperature replica-exchange (TRXLD) probability; the special case of eq:20 with all SGLD factors = 1 ($\tilde{\mu}=0$, $\mu=\beta$).
- **symbols:** $\beta_m = 1/(kT_m)$, $\beta_n = 1/(kT_n)$; $E_p$ - potential energy.

<!-- eq:22 -->
$$ \bar{\alpha}(t) = \left(1 - \frac{\delta t}{t_{\text{est}}}\right) \bar{\alpha}(t - \delta t) + \frac{\delta t}{t_{\text{est}}} \alpha(t) $$
- **what:** evolving (exponential moving) average used to estimate $\lambda_{lf},\lambda_{hf},\chi_{lf},\chi_{hf}$ on the fly during simulation.
- **symbols:** $\alpha(t)$ - instantaneous value; $\bar{\alpha}(t)$ - estimated average; $t_{est}$ - estimation time (typically $10\,t_L$); $\delta t$ - time step.

<!-- eq:23a -->
$$ \mathbf{p}'^{(n)}_i = s_{mn}\mathbf{p}_i^{(m)} $$
- **what:** momentum rescaling of a replica moving from stage m to stage n on accepted exchange.
- **symbols:** $\mathbf{p}_i^{(m)}$ - momentum at stage m; $s_{mn}$ - temperature-scaling factor (eq:24).

<!-- eq:23b -->
$$ \mathbf{p}'^{(m)}_i = s_{nm}\mathbf{p}_i^{(n)} $$
- **what:** momentum rescaling of the paired replica moving from stage n to stage m.
- **symbols:** $s_{nm} = 1/s_{mn}$.

<!-- eq:24 -->
$$ s_{mn} = \sqrt{\frac{T_n}{T_m}} = \frac{1}{s_{nm}} $$
- **what:** momentum temperature-scaling factor between stages m and n.
- **symbols:** $T_m,T_n$ - stage temperatures.

<!-- eq:25a -->
$$ \tilde{\mathbf{p}}'^{(n)}_i = \tilde{s}_{mn}\tilde{\mathbf{p}}_i^{(m)} $$
- **what:** low-frequency-momentum rescaling on accepted exchange (m -> n).
- **symbols:** $\tilde{\mathbf{p}}_i$ - low-frequency momentum; $\tilde{s}_{mn}$ - low-frequency scaling factor (eq:26).

<!-- eq:25b -->
$$ \tilde{\mathbf{p}}'^{(m)}_i = \tilde{s}_{nm}\tilde{\mathbf{p}}_i^{(n)} $$
- **what:** low-frequency-momentum rescaling for the paired replica (n -> m).
- **symbols:** $\tilde{s}_{nm} = 1/\tilde{s}_{mn}$.

<!-- eq:26 -->
$$ \tilde{s}_{mn} = \sqrt{\frac{\tilde{T}_n}{\tilde{T}_m}} = \frac{1}{\tilde{s}_{nm}} $$
- **what:** low-frequency momentum temperature-scaling factor between stages.
- **symbols:** $\tilde{T}_m,\tilde{T}_n$ - stage low-frequency temperatures (eq:12).

<!-- eq:27 -->
$$ \varepsilon_{\rm p}(x,y,z) = \frac{a}{w^2}(x^2 + z^2) + \frac{b}{w^4}y^2(y-w)^2 + \frac{s}{w}y $$
- **what:** skewed double-well test potential (single particle). Harmonic in x,z; double well in y (wells at y=0 and y=w); skew s tilts the wells.
- **symbols:** $a$ - x,z stiffness; $b$ - double-well depth; $w$ - well separation; $s$ - skew (energy offset). Test values: $a=20000kT_0$, $b=160kT_0$, $w=2\,\text{Å}$, $s\in\{0,kT_0,2kT_0\}$, $T_0=50\,\text{K}$.

<!-- eq:28a -->
$$ \Theta_{xz} = \frac{\pi k T w^2}{a} $$
- **what:** configurational partition function for the (harmonic) x,z degrees of freedom of eq:27.
- **symbols:** as eq:27.

<!-- eq:28b -->
$$ \Theta_y = \int_{-\infty}^{\infty} \exp\left(-\frac{\frac{b}{w^4}y^2(w-y)^2 + \frac{s}{w}y}{kT}\right) dy $$
- **what:** configurational partition function for the y degree of freedom (numerically integrated).
- **symbols:** as eq:27.

<!-- eq:29a -->
$$ E_{xz} = kT $$
- **what:** analytic average energy of the x,z (2 harmonic DOF) subsystem.
- **symbols:** $k$ - Boltzmann constant; $T$ - temperature.

<!-- eq:29b -->
$$ E_y = \frac{1}{\Theta_y} \int_{-\infty}^{\infty} \left( \frac{b}{w^4} y^2 (w - y)^2 + \frac{s}{w} y \right) \exp\left(-\frac{\frac{b}{w^4}y^2(w-y)^2 + \frac{s}{w}y}{kT}\right) dy $$
- **what:** analytic (numerically integrated) average energy of the y degree of freedom.
- **symbols:** $\Theta_y$ from eq:28b; other symbols per eq:27. <!-- CHECK: raw dropped the 1/kT inside the exponent of the Boltzmann weight; restored to match eq:28b/eq:30b. -->

<!-- eq:30a -->
$$ \rho_{xz}(r_{xz}) = \frac{2a}{kTw^2} e^{-\frac{a r_{xz}^2}{kTw^2}} $$
- **what:** radial distribution in the x-z plane, $r_{xz}=\sqrt{x^2+z^2}$.
- **symbols:** $r_{xz}$ - radial coordinate; other symbols per eq:27.

<!-- eq:30b -->
$$ \rho_y(y) = \frac{1}{\Theta_y} \exp\left(-\frac{\frac{b}{w^4}y^2(w-y)^2 + \frac{s}{w}y}{kT}\right) $$
- **what:** Boltzmann distribution along y.
- **symbols:** $\Theta_y$ from eq:28b; other symbols per eq:27.

<!-- eq:31 -->
$$ I_i = \begin{cases} 1 & s_i \in R_i(1) \\ 2 & s_i \in R_i(2) \\ \dots & \dots \\ k_i & s_i \in R_i(k_i) \end{cases} $$
- **what:** region index of subset variable $s_i$ for subset-indexing clustering (SIC).
- **symbols:** $s_i$ - subset variable (e.g. a dihedral angle); $R_i(1..k_i)$ - the $k_i$ regions of subset i; $I_i$ - resulting index. A cluster is $\{I_1,...,I_m\}$; total clusters $N_c=\prod_{i=1}^m k_i$.
