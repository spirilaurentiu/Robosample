# Equations - Kofke 2002 (replica-exchange acceptance probability)

NOTE: Eq. (12) and Eq. (13) are the erratum-corrected forms
(J. Chem. Phys. 120, 10852, 2004). Both are flagged with `<!-- CHECK -->`.

<!-- eq:1 -->
$$ P_{\rm acc}(\text{swap } 0,1) = \min[1, \exp(-(\beta_0 - \beta_1)(U_1 - U_0))] $$
- **what:** Metropolis acceptance for a single replica-exchange (parallel-tempering) swap of two configurations between temperatures.
- **symbols:** $\beta_i = 1/kT_i$ - reciprocal temperature of system $i$ (scalar); $U_i$ - potential energy of system $i$ before the swap (energy); $k$ - Boltzmann constant. Configurations are swapped, temperatures held fixed.

<!-- eq:2 -->
$$ \bar{p}_{\rm acc} \sim \exp(-\Delta S/k) $$
- **what:** Average swap-acceptance probability scales with the inter-replica entropy difference (entropy-model estimate).
- **symbols:** $\bar{p}_{\rm acc}$ - fraction of all swap trials accepted (dimensionless); $\Delta S > 0$ - entropy of the high-temperature system minus that of the low-temperature system; $k$ - Boltzmann constant.

<!-- eq:3 -->
$$ \bar{p}_{\rm acc} \sim \left(\frac{T_0}{T_1}\right)^{N c_V/k} $$
- **what:** Constant-heat-capacity form of the entropy-model estimate; acceptance decays as a power of the temperature ratio, exponent extensive in $N$.
- **symbols:** $T_0 < T_1$ - low/high temperatures; $N$ - number of molecules; $c_V$ - molar constant-volume heat capacity; $k$ - Boltzmann constant. Uses $T\Delta S = N c_V \Delta T$.

<!-- eq:4 -->
$$ \bar{p}_{acc} = \int_{U_m}^{\infty} dU_1\, p_1(U_1) \int_{U_m}^{\infty} dU_0\, p_0(U_0) \times \min[1, \exp(-(\beta_0 - \beta_1)(U_1 - U_0))] $$
- **what:** Exact average acceptance probability as a double integral over the two replicas' energy distributions weighted by the swap acceptance.
- **symbols:** $p_i(U)$ - energy distribution of system at temperature $\beta_i$ (probability density in energy); $U_m$ - lowest possible energy; other symbols as in eq:1.

<!-- eq:5 -->
$$ p_i(U) = \frac{1}{Q(\beta_i)} \Omega(U) \exp(-\beta_i U) $$
- **what:** Canonical decomposition of the energy distribution into density of states times Boltzmann factor over the partition function.
- **symbols:** $\Omega(U)$ - density of states, $\beta$-independent; $Q(\beta_i) = \int \Omega(U) e^{-\beta_i U}\,dU$ - canonical partition function; $p_i(U)$ - normalized energy distribution.

<!-- eq:7 -->
$$ \bar{p}_{\text{acc}} = 2 \int_{U_m}^{\infty} dU_0\, p_0(U_0) \int_{U_m}^{U_0} dU_1\, p_1(U_1) $$
- **what:** Exact acceptance probability rewritten as (twice) the overlap of the two energy distributions; ranges 0 (no overlap) to 1 (identical distributions). Requires $\beta_0 > \beta_1$.
- **symbols:** $p_0$ - low-temperature (cold) distribution; $p_1$ - high-temperature (hot) distribution, peaks at larger $U$; $U_m$ - lowest possible energy. Inner integral runs $U_1 \in [U_m, U_0]$.

<!-- eq:8 -->
$$ \Omega(U) = \left(1 + \frac{1}{C}\beta_r(U - U_r)\right)^C \Omega(U_r) $$
- **what:** Density of states for a system of constant heat capacity across the temperature range.
- **symbols:** $C \equiv C_V/k$ - extensive constant-volume heat capacity in units of $k$ (dimensionless); $U_r, T_r, \beta_r=1/kT_r$ - arbitrary reference state; $\Omega(U_r)$ - density of states at the reference.

<!-- eq:9 -->
$$ p_i(U) = \frac{\beta_i}{C\,\Gamma(C)} [\beta_i(U - U_m)]^C \exp[-\beta_i(U - U_m)] $$
- **what:** Normalized energy distribution for the constant-heat-capacity model (a gamma distribution in $U - U_m$ with shape $C+1$).
- **symbols:** $\Gamma$ - gamma function; $U_m = U_r - C/\beta_r$ - minimum possible energy; $C$, $\beta_i$ as above.

<!-- eq:sub -->
$$ \kappa \equiv \frac{\beta_1(U_1 - U_m)}{\beta_0(U_0 - U_m)}, \qquad \gamma = \beta_1(U_1 - U_m) + \beta_0(U_0 - U_m) $$
- **what:** Change of variables used to reduce the eq:7/eq:9 double integral to the single integral eq:10.
- **symbols:** $\kappa$ - ratio variable (integration variable of eq:10, range $[0, \beta_1/\beta_0]$); $\gamma$ - sum variable (integrated out).

<!-- eq:10 -->
$$ \bar{p}_{acc} = 4 \frac{\Gamma(2C)}{[\Gamma(C)]^2} \frac{2C+1}{C} \int_0^{\beta_1/\beta_0} d\kappa\, \frac{\kappa^C}{(1+\kappa)^{2(C+1)}} $$
- **what:** Exact acceptance probability under the constant-heat-capacity model; a one-dimensional integral evaluated numerically. Depends on temperatures only through the ratio $\beta_1/\beta_0$.
- **symbols:** $C \equiv C_V/k$ - heat capacity in units of $k$; upper limit $\beta_1/\beta_0 = T_0/T_1 \le 1$ (since $\beta_1 < \beta_0$); $\Gamma$ - gamma function.

<!-- eq:11 -->
$$ \Delta S/k = -C \ln(\beta_1/\beta_0) $$
- **what:** Entropy difference in units of $k$ for the constant-heat-capacity model; positive since $\beta_1 < \beta_0$.
- **symbols:** $\Delta S$ - hot-minus-cold entropy difference; $C \equiv C_V/k$; $\beta_1/\beta_0 = T_0/T_1 < 1$.

<!-- eq:12 -->
<!-- CHECK: erratum-corrected (J. Chem. Phys. 120, 10852, 2004). Original 2002 text had (1+B)/(1-2B+B^2)=(1+B)/(1-B)^2; erratum restores a 1/2 power on that denominator, giving (1+B)/(1-B). -->
$$ \bar{p}_{\rm acc} \sim \frac{\exp(-\Delta S/k)}{(\pi C)^{1/2}} \left[ \frac{4}{(1+B)^2} \right]^{C+1} \frac{1+B}{1-B} \left(1 + O(C^{-1/2})\right), \quad C \to \infty $$
- **what:** Large-system asymptotic form of the exact acceptance probability eq:10; the bracketed $C+1$ term (>1) attenuates the exponential decay with size. NOTE (usage): accurate only away from $B \to 1$. For closely-spaced temperatures ($\Delta S/k$ down to a few $k$) this asymptotic form breaks down; use the exact integral eq:10 or the Gaussian eq:13 instead, which Fig. 2 (erratum) shows track the simulation data down to small $\Delta S$. Do not use eq:12 to set a tight temperature ladder near unity.
- **symbols:** $B = \beta_1/\beta_0 < 1$ (equals $T_0/T_1$); $C \equiv C_V/k$; $\Delta S/k = -C\ln B$ from eq:11. Valid for $B$ not too close to 1.

<!-- eq:13 -->
<!-- CHECK: erratum-corrected (J. Chem. Phys. 120, 10852, 2004). Original 2002 text omitted the (1+(beta1/beta0)^2)^{1/2} denominator, giving erfc[(C/2)^{1/2}(1-beta1/beta0)]. -->
$$ \bar{p}_{\text{acc}} = \text{erfc}\left[ \left(\tfrac{1}{2}C\right)^{1/2} \frac{1 - \beta_1/\beta_0}{\left(1 + (\beta_1/\beta_0)^2\right)^{1/2}} \right] $$
- **what:** Acceptance probability when both replicas' energy distributions are approximated as Gaussians. The original 2002 formula was in error (missing denominator term); the erratum correction restores agreement, and the erratum states it "gives more credence to the use of a Gaussian form" — the corrected erfc tracks the data down to small $\Delta S$ (Fig. 2), so it is the preferred formula near unity where eq:12 fails.
- **symbols:** $\text{erfc}$ - complementary error function; $C \equiv C_V/k$; $\beta_1/\beta_0 = T_0/T_1 < 1$.
