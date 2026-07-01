# Equations - Minh 2019 (AlGDock)

Reduced units: energies expressed as reduced potential $u = U/(k_B T)$; free energies $f$ in units of $k_B T$. $\beta = (k_B T)^{-1}$.

<!-- eq:1 -->
$$B(r_R) = -\beta^{-1} \ln \left( \frac{\int I(\xi) J(\xi) e^{-\beta U(r_{RL})} \, dr_L \, d\xi}{\int I(\xi) J(\xi) e^{-\beta [U(r_L) + U(r_R)]} \, dr_L \, d\xi} \right)$$
- **what:** Binding potential of mean force (BPMF): a ratio of configurational integrals over the bound complex vs. separated receptor and ligand.
- **symbols:** $B(r_R)$ - BPMF for receptor conformation $r_R$ (energy, $k_B T /\beta$); $r_{RL}$ - complex internal coords; $r_R$ - receptor internal coords; $r_L$ - ligand internal coords; $\xi$ - relative translation+rotation; $I(\xi)$ - bound indicator (1 bound, 0 unbound); $J(\xi)$ - Jacobian of the coordinate transform; $U(\cdot)$ - potential energy of a species in solvent; $\beta=(k_B T)^{-1}$.

<!-- eq:2 -->
$$u_I(d) = \begin{cases} 0 & \text{if } d \le d_0 \\ \frac{1}{2}\beta k(d - d_0)^2 & \text{if } d > d_0 \end{cases}$$
- **what:** Flat-bottom harmonic restraint confining the ligand center of mass to the binding site (reduced units).
- **symbols:** $u_I$ - reduced restraint energy ($k_B T$); $d$ - distance from ligand COM to binding-site center (nm); $d_0 = 6.0$ Å = 0.6 nm - binding-site radius; $k = 10000$ kJ/(mol nm$^2$) - spring constant; $\beta=(k_B T)^{-1}$. <!-- CHECK: original writes d_o once; treated as d_0 throughout -->

<!-- eq:3 -->
$$\Psi_g(r_{RL}) = \Psi_{PBSA}(r_{RL}) + \Psi_{vdW}(r_{RL})$$
- **what:** Total grid-interpolated receptor-ligand interaction energy = electrostatic (Poisson-Boltzmann) + van der Waals.
- **symbols:** $\Psi_g$ - grid interaction energy (kJ/mol); $\Psi_{PBSA}$ - electrostatic (partial charges times interpolated PB potential); $\Psi_{vdW}$ - van der Waals (transformed/interpolated grid).

<!-- eq:4 -->
$$u_{\alpha}(r_{RL}) = \frac{1}{k_B T(\alpha)} \left[ U_s(r_L) + U_s(r_R) + \alpha_{sg}(\alpha)\,\Psi_{sg}(r_{RL}) + \alpha_g(\alpha)\,\Psi_g(r_{RL}) \right]$$
$$\alpha_{sg}(\alpha) = -(2\alpha - 1)^2 + 1$$
$$\alpha_g(\alpha) = \frac{(2\alpha - 1)^2}{1 + \exp\left[-1000(\alpha - \tfrac{1}{2})\right]}$$
$$T(\alpha) = (T_T - T_H)\alpha + T_H$$
- **what:** Reduced-potential switching protocol for states CD (progress variable $\alpha$): turns on soft grids first, then unperturbed grids, while ramping temperature. Consistent with milestone C at $\alpha=0$, milestone D at $\alpha=1$.
- **symbols:** $u_\alpha$ - reduced potential at $\alpha$ ($k_B T$); $\alpha \in [0,1]$ - progress variable; $U_s$ - sampling force-field potential (ligand $r_L$ and receptor $r_R$); $\Psi_{sg}$ - soft-grid interaction energy; $\Psi_g$ - unperturbed grid energy; $\alpha_{sg}$ - soft-grid scaling (parabola, peaks 1 at $\alpha=1/2$); $\alpha_g$ - unperturbed-grid scaling (sigmoid-gated parabola); $T(\alpha)$ - linear temperature ramp; $T_T=300$ K (target), $T_H=600$ K (high). Note $T(0)=T_H$, $T(1)=T_T$.

<!-- eq:5 -->
$$p_{acc} = \min \left[ 1, \; e^{-u_a(y) - u_b(x) + u_a(x) + u_b(y)} \right]$$
- **what:** Hamiltonian replica-exchange (Metropolis) acceptance probability for swapping configurations $x$ (state $a$) and $y$ (state $b$). Preserves Boltzmann distribution in both states.
- **symbols:** $p_{acc}$ - acceptance probability; $u_a, u_b$ - reduced energies of states $a, b$; $x$ - original config in state $a$; $y$ - original config in state $b$. Exponent $= -[u_a(y)+u_b(x)] + [u_a(x)+u_b(y)]$ = negative change in total reduced energy upon swap.

<!-- eq:6 -->
$$\mathcal{L} \equiv \int_0^1 \sqrt{\sum_{i,j} \frac{\partial \gamma^i}{\partial \alpha} \, g(\gamma)_{ij} \, \frac{\partial \gamma^j}{\partial \alpha}} \; d\alpha$$
- **what:** Thermodynamic length along a protocol path $\gamma(\alpha)$ through parameter space; states equidistant in $\mathcal{L}$ minimize statistical error.
- **symbols:** $\mathcal{L}$ - thermodynamic length (dimensionless); $\alpha \in [0,1]$ - path parameter; $\gamma^i(\alpha)$ - $i$-th parameter along path ($\gamma(0)$ initial, $\gamma(1)$ final); $g(\gamma)_{ij}$ - metric tensor (Fisher information, eq:7).

<!-- eq:7 -->
$$g(\gamma)_{ij} \equiv \sigma_\lambda^2 \left[ \partial^i l_\lambda, \; \partial^j l_\lambda \right]$$
- **what:** Fisher-information metric tensor: covariance of derivatives of the log-probability w.r.t. parameters.
- **symbols:** $g_{ij}$ - metric element; $\sigma_\lambda^2[\cdot,\cdot]$ - covariance in state $\lambda$; $l_\lambda(x) = -u_\lambda(x) - \ln Z_\lambda$ - normalized log probability; $u_\lambda(x)=U_\lambda(x)/(k_B T_\lambda)$ - reduced potential; $Z_\lambda=\int e^{-u_\lambda(x)}dx$ - partition function; $\partial^i \equiv \partial/\partial\lambda^i$.

<!-- eq:8 -->
$$\Delta\lambda^i = \frac{s}{\sigma_0\left[\partial^i l_\lambda\right]}$$
- **what:** Single-parameter step rule to keep thermodynamic length constant between adjacent states.
- **symbols:** $\Delta\lambda^i$ - change in parameter $\lambda^i$ between adjacent states; $s$ - thermodynamic speed (adjustable); $\sigma_0[\partial^i l_\lambda]$ - standard deviation of $\partial^i l_\lambda$ in the initial state.

<!-- eq:9 -->
$$\Delta\lambda^i = -\frac{s_{bc}\,k T^2}{\sigma_\lambda[U_S]}$$
- **what:** Temperature-step rule for states BC (Full pathway), specializing eq:8 with $\lambda^i = T$ and $l_\lambda = -U_S(r_L)/(k_B T) - \ln Z_\lambda$.
- **symbols:** $\Delta\lambda^i$ - temperature increment (K); $s_{bc}=20.0$ - thermodynamic speed for BC; $k$ - Boltzmann constant; $T$ - current temperature (K); $\sigma_\lambda[U_S]$ - std dev of sampling potential energy of the ligand in state $\lambda$. <!-- CHECK: sign/units: increment increases T; k here is k_B -->

<!-- eq:10 -->
$$\Delta\lambda^i = s_{cd} \left[ \left| \frac{d\alpha_{sg}}{d\alpha} \right| \frac{\sigma_\lambda[\Psi_{sg}]}{k_B T(\alpha)} + \left| \frac{d\alpha_g}{d\alpha} \right| \frac{\sigma_\lambda[\Psi_g]}{k_B T(\alpha)} + |T_T - T_H| \frac{\sigma_\lambda[u_\alpha(r_{RL})]}{T(\alpha)} \right]^{-1}$$
- **what:** Progress-variable step rule for states CD; accounts for soft-grid, unperturbed-grid, and temperature contributions to the metric.
- **symbols:** $\Delta\lambda^i$ - increment in $\alpha$; $s_{cd}=0.2$ - thermodynamic speed for CD; $\alpha_{sg},\alpha_g$ - scaling functions (eq:4); $\sigma_\lambda[\cdot]$ - std dev in state $\lambda$; $\Psi_{sg},\Psi_g$ - soft/unperturbed grid energies; $u_\alpha$ - reduced potential; $T_T=300$ K, $T_H=600$ K; $T(\alpha)$ - eq:4.

<!-- eq:11 -->
$$f'_{CD} = -\ln \frac{\int I(\xi)J(\xi)e^{-\beta_T[U(r_L)+\Psi_g(r_{RL})]}\,dr_L\,d\xi}{\int I(\xi)J(\xi)e^{-\beta_H U(r_L)}\,dr_L\,d\xi}$$
- **what:** Reduced free energy $f'_{CD}$; used with $f_{BC,L}$ (ligand warming) to avoid computing the receptor internal energy $U(r_R)$. Estimated by MBAR.
- **symbols:** $f'_{CD}$ - reduced free energy ($k_B T$); $\beta_T=(k_B\cdot300\,\text{K})^{-1}$; $\beta_H=(k_B\cdot600\,\text{K})^{-1}$; $\Psi_g$ - grid interaction (eq:3); $I,J,U$ as in eq:1.

<!-- eq:12 -->
$$w_c = \exp\left[-\beta_T \left(U_T(r_{RL}) - U_S(r_L) - \Psi_g(r_L)\right)\right]$$
- **what:** Configuration reweighting factor from milestone D to E (total-energy variant, using sampling ligand energy).
- **symbols:** $w_c$ - unnormalized weight of configuration $c$; $\beta_T=(k_B\cdot300\,\text{K})^{-1}$; $U_T$ - target force-field energy of complex $r_{RL}$; $U_S(r_L)$ - sampling force-field ligand energy; $\Psi_g(r_L)$ - grid interaction energy. <!-- CHECK: paper writes Psi_q; interpreted as Psi_g (grid); argument r_L per source -->

<!-- eq:13 -->
$$w_c = \exp\left[-\beta_T \left(U_T(r_{RL}) - U_T(r_L) - \Psi_g(r_L)\right)\right]$$
- **what:** Reweighting factor assuming only interaction energies change between milestones D and E (uses target ligand energy).
- **symbols:** as eq:12, with $U_T(r_L)$ - target force-field ligand energy.

<!-- eq:14 -->
$$f_{EE_p} = -\ln \frac{\sum_c w_c}{\sum_p \sum_c w_c}$$
- **what:** Reduced free energy of pose $p$ from cumulative cluster weight.
- **symbols:** $f_{EE_p}$ - reduced free energy of pose $p$ ($k_B T$); $\sum_c$ - sum over configurations in the cluster; $\sum_p$ - sum over poses; $w_c$ - weight (eq:12/13).

<!-- eq:15 -->
$$f_{AE,p} = f_{AE} + f_{EE_p}$$
- **what:** Pose-specific BPMF.
- **symbols:** $f_{AE,p}$ - pose-specific reduced BPMF; $f_{AE}$ - overall reduced BPMF (A->E); $f_{EE_p}$ - pose free energy (eq:14).

## Derivations note

The thermodynamic-length material (eq:6-8) is derivation-heavy; the single-parameter approximation $\mathcal{L} = \Delta\lambda^i\,\sigma_0[\partial^i l_\lambda]$ and its specializations (eq:9, eq:10) are the implementable step rules. The intermediate relation $\mathcal{L} = \int_0^1 \frac{\partial\gamma^i}{\partial\alpha}\sigma_\lambda[\partial^i l_\lambda]\,dt$ (single varying parameter) is a proof step, not separately implemented.
