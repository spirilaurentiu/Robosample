# Equations - Shirts 2012, Ensemble validation tests

Core idea: for any Boltzmann-form ensemble, the log-ratio of the probability
distributions of an extensive variable between two simulations that differ only
in an intensive parameter is a **known linear function** of that extensive
variable, with slope fixed by the intensive parameters and *independent of the
unknown density of states*. Fitting that slope validates the sampled ensemble.

## Ensemble probability distributions

<!-- eq:1 -->
$$ P(\vec{x}|\beta) \propto \exp(-\beta H(\vec{p}, \vec{q})) $$
- **what:** canonical (NVT) microstate probability.
- **symbols:** x - full microstate; beta = 1/(k_B T) - inverse temperature; H - Hamiltonian; p - momenta; q - positions.

<!-- eq:2 -->
$$ P(\vec{x}, V | \beta, P) \propto \exp(-\beta (H(\vec{p}, \vec{q}) + PV)) $$
- **what:** isobaric-isothermal (NPT) microstate probability.
- **symbols:** V - volume; P - pressure; PV - pressure-volume work term.

<!-- eq:3 -->
$$ P(\vec{x}, \vec{N}|\beta, \vec{\mu}) \propto \exp\left(-\beta\left(H(\vec{p}, \vec{q}) - \sum_{\text{species}} \mu_i N_i\right)\right) $$
- **what:** grand canonical microstate probability.
- **symbols:** N_i - particle number of species i; mu_i - chemical potential of species i.

## Canonical energy distribution and validation relation

<!-- eq:4 -->
$$ P(E|\beta) = Q(\beta)^{-1}\,\Omega(E)\,\exp(-\beta E) $$
- **what:** probability density of total energy E in the canonical ensemble.
- **symbols:** Omega(E) = exp(S(N,V,E)/k_B) - density of states (function of E, NOT beta); Q(beta) = integral of Omega(E) exp(-beta E) dE - canonical partition function (function of beta, NOT E); S - entropy; A = -beta^-1 ln Q - Helmholtz free energy.

<!-- eq:5 -->
$$ \frac{P(E|\beta_2)}{P(E|\beta_1)} = \frac{\exp(-\beta_2 E)/Q(\beta_2)}{\exp(-\beta_1 E)/Q(\beta_1)} = \exp\big([\beta_2 A_2 - \beta_1 A_1] - [\beta_2 - \beta_1]E\big) $$
- **what:** ratio of two canonical energy distributions at different beta; the unknown Omega(E) cancels.
- **symbols:** A_1, A_2 - Helmholtz free energies at beta_1, beta_2.

<!-- eq:6 -->
$$ \ln \frac{P(E|\beta_2)}{P(E|\beta_1)} = [\beta_2 A_2 - \beta_1 A_1] - [\beta_2 - \beta_1]E $$
- **what:** THE core ensemble-consistency test. Linear in E: form alpha_0 + alpha_1 E with intercept alpha_0 = beta_2 A_2 - beta_1 A_1 and slope alpha_1 = -(beta_2 - beta_1). The slope is independent of the (unknown) free energies. A fitted slope deviating from -(beta_2-beta_1) implies the data is NOT canonical.
- **symbols:** alpha_0 - intercept (contains free energies); alpha_1 = -(beta_2 - beta_1) - slope (known from user temperatures).

## Fitting the slope

<!-- eq:7 -->
$$ S_r(\alpha_0, \alpha_1) = \sum_i \left[ \frac{P_1(E_i)}{P_2(E_i)} - \exp(\alpha_0 + \alpha_1 E_i) \right]^2 $$
- **what:** nonlinear least-squares objective: minimize sum of squared residuals between the histogram ratio and the model exp(alpha_0 + alpha_1 E).
- **symbols:** P_1(E_i), P_2(E_i) - histogram probabilities in bin at energy E_i for the two simulations; alpha_0, alpha_1 - fit parameters.
<!-- CHECK: eq:6 uses P2/P1 giving slope -(b2-b1); eq:7 residual is written with P1/P2 = exp(alpha0+alpha1 E), so its alpha1 has opposite sign convention. Keep both conventions straight when implementing. -->

<!-- eq:8 -->
$$ \ln L(\vec{\alpha}|\text{data}) = \sum_{i=1}^{N_1} \ln f(-\alpha_0 - \alpha_1 E_i) + \sum_{j=1}^{N_2} \ln f(\alpha_0 + \alpha_1 E_j) $$
- **what:** maximum-likelihood (logistic-regression) objective for the two fit parameters. Histogram-free. First sum over energies at T_1, second over energies at T_2. Concave with a single maximum -> always converges.
- **symbols:** f(x) = [1 + exp(-x)]^{-1} - Fermi/logistic function; N_1, N_2 - number of samples from sim 1 and 2; alpha_0 = ln(Q_1/Q_2) = beta_2 A_2 - beta_1 A_1; alpha_1 = -(beta_2 - beta_1).

### Reparameterization to two free parameters (canonical)
$$ \alpha_0 + \alpha_1 E = \beta_{\text{ave}}\Delta A - \Delta\beta\, E $$
- with beta_ave = (beta_{1,user}+beta_{2,user})/2 held fixed, DeltaA = A_2 - A_1 (setting A_1 + A_0 = 0), Delta_beta = beta_2 - beta_1.

## Covariance / error estimators

<!-- eq:cov -->
$$ \text{cov}(\vec{\alpha})_{\text{lin}} = (X^T W X)^{-1}, \quad \text{cov}(\vec{\alpha})_{\text{nonlin}} = (J^T W^{-1} J)^{-1}, \quad \text{cov}(\vec{\alpha})_{\text{ML}} = \big(\text{Hess}(\ln L)_{\alpha}\big)^{-1} $$
- **what:** asymptotic covariance of the fit parameters for linear, nonlinear, and maximum-likelihood fits.
- **symbols:** X - (M+1) x N design matrix (first column ones, cols 2..M+1 are the N observations of the M observables); J - Jacobian of the model wrt alpha at the minimum; W - diagonal weight matrix of per-point variances; Hess(ln L) - Hessian of log-likelihood wrt alpha.
<!-- CHECK: main text writes nonlinear cov as (J^T W^-1 J)^-1 but appendix eq (below) writes (J^T W J)^-1; W is defined as inverse-variance weights in the appendix, so (J^T W J)^-1 is the consistent form. -->

## Harmonic-oscillator toy model (canonical)

<!-- eq:ho-energy -->
$$ E = \tfrac{1}{2} K \sum_{i=1}^{D} (x_i - x_{i,0})^2 $$
- **what:** potential energy of a D-dimensional harmonic oscillator, equal spring constant K per dimension.
- **symbols:** K - spring constant; x_{i,0} - equilibrium position in dimension i; D - dimensionality.

<!-- eq:ho-Q -->
$$ Q(\beta) = \left(\frac{2\pi}{\beta K}\right)^{D/2}, \qquad A(\beta) = -\frac{D}{2\beta}\ln\!\left(\frac{2\pi}{\beta K}\right) $$
- **what:** partition function and Helmholtz free energy of the D-dim harmonic oscillator.

<!-- eq:ho-P -->
$$ P(\vec{x}|\beta) = \left(\frac{\beta K}{2\pi}\right)^{D/2} \exp\left(-\frac{\beta K}{2} \sum_{i=1}^{D} |x_i - x_{i,0}|^2\right) $$
- **what:** configurational distribution; each coordinate is Gaussian with variance 1/(beta K).

Average energy: $\langle E\rangle = -\partial \ln Q/\partial\beta = D/(2\beta)$.

## Kinetic / potential separability (NVT)

<!-- eq:sep -->
$$ P(E_{\text{pot}} + E_{\text{kin}}|\beta) = P(E_{\text{pot}}|\beta)\,P(E_{\text{kin}}|\beta) $$
- **what:** in NVT the total-energy distribution factorizes into independent kinetic and potential parts (momenta sample independently of coordinates), so eq:6 holds for E_kin and E_pot separately.

<!-- eq:19 -->
$$ \frac{P(E_{\text{kin}}|\beta_2)}{P(E_{\text{kin}}|\beta_1)} = \frac{Q_{\text{kin}}(\beta_2)}{Q_{\text{kin}}(\beta_1)} \exp(-[\beta_2 - \beta_1]E_{\text{kin}}) $$
- **what:** kinetic-energy validation relation.

<!-- eq:20 -->
$$ \frac{P(E_{\text{pot}}|\beta_2)}{P(E_{\text{pot}}|\beta_1)} = \frac{Q_{\text{pot}}(\beta_2)}{Q_{\text{pot}}(\beta_1)} \exp(-[\beta_2 - \beta_1]E_{\text{pot}}) $$
- **what:** potential-energy validation relation (usable for Monte Carlo, where only E_pot is defined).

<!-- eq:21 -->
$$ \frac{P(E_{\text{kin}}|\beta_2)}{P(E_{\text{kin}}|\beta_1)} = \left(\frac{\beta_2}{\beta_1}\right)^{3N/2} \exp([\beta_1 - \beta_2]E_{\text{kin}}) $$
- **what:** explicit kinetic ratio using Q_kin = prod_i (m_i/(pi beta))^{3/2}; mass terms cancel so it holds for identical and non-identical particles.
- **symbols:** N - number of particles; 3N - kinetic DOF (replace by actual DOF if constraints/COM removal).

<!-- eq:22 -->
$$ = \left(\frac{\beta_{\text{ave}} + \Delta\beta}{\beta_{\text{ave}} - \Delta\beta}\right)^{3N/2} \exp(-\Delta\beta\, E_{\text{kin}}) $$
- **what:** eq:21 rewritten in the single free parameter Delta_beta = beta_2 - beta_1 (beta_ave = (beta_1+beta_2)/2).

<!-- eq:kin-dist -->
$$ P(E_{\text{kin}}) = \frac{\beta}{\sqrt{3N\pi}} \exp\left(-\frac{\left(\beta E_{\text{kin}} - \frac{3N}{2}\right)^2}{3N}\right) $$
- **what:** Gaussian approximation to the chi-squared kinetic-energy distribution (valid for >~60 DOF ~20 particles). Mean 3N/(2 beta), variance sigma^2 = 3N/(2 beta^2).
- **symbols:** mean <E_kin> = 3N/(2 beta) (equipartition); sigma^2 = k_B T^2 C_V with kinetic C_V = 3 N k_B / 2 giving sigma^2 = 3N/(2 beta^2). Replace 3N by actual DOF.

Equipartition DOF estimate: $\langle E_{\text{kin}}\rangle = \tfrac{k_B T}{2}(\#\text{DOF})$.

## Isobaric-isothermal (NPT) validation relations

<!-- eq:npt-H -->
$$ P(H|\beta,P) = \frac{\beta P}{h^{3N}} \Omega'(H,P)\, \Delta(\beta,P)^{-1} \exp(-\beta H) $$
- **what:** enthalpy distribution at fixed P; H = E + PV instantaneous enthalpy.
- **symbols:** Delta(beta,P) - isothermal-isobaric partition function; Omega'(H,P) - density of states counting states with E+PV=H (function of P, not beta); h - Planck constant; beta P prefactor cancels units.

<!-- eq:9 -->
$$ \frac{P(H|\beta_2, P)}{P(H|\beta_1, P)} = \frac{\beta_1 \Delta(\beta_1, P)}{\beta_2 \Delta(\beta_2, P)} \exp(-[\beta_2 - \beta_1]H) $$
- **what:** enthalpy ratio for two sims at same P, different T. Slope -(beta_2-beta_1) on H.

<!-- eq:10 -->
$$ \ln\frac{P(H|\beta_2, P)}{P(H|\beta_1, P)} = \ln(\beta_1/\beta_2) + [\beta_2 G_2 - \beta_1 G_1] - [\beta_2 - \beta_1]H $$
- **what:** log form; Gibbs free energy G replaces A (plus ln(beta_1/beta_2) correction).
- **symbols:** G_1, G_2 - Gibbs free energies.

<!-- eq:11 -->
$$ \ln \frac{P(V|\beta, P_2)}{P(V|\beta, P_1)} = \ln(P_1/P_2) + \beta(G_2 - G_1) - \beta(P_2 - P_1)V $$
- **what:** volume validation: two sims at same T, different P. Slope on V is -beta(P_2 - P_1).

<!-- eq:14 -->
$$ \ln\frac{P(V,E|\beta_2,P_2)}{P(V,E|\beta_1,P_1)} = \ln(\beta_1 P_1/\beta_2 P_2) + [\beta_2 G_2 - \beta_1 G_1] - [\beta_2 - \beta_1]E - [\beta_2 P_2 - \beta_1 P_1]V $$
- **what:** joint (E,V) validation for two sims differing in both T and P; a multilinear fit in E and V.

<!-- eq:15 -->
$$ (\beta_2 P_2 - \beta_1 P_1) = \tfrac{1}{2}\big((\Delta\beta)(P_2 + P_1) + (\beta_2 + \beta_1)(\Delta P)\big) = (\Delta\beta) P_{\text{ave}} + \beta_{\text{ave}}(\Delta P) $$
- **what:** decomposition of the V-coefficient into Delta_beta and Delta_P pieces.
- **symbols:** Delta_beta = beta_2 - beta_1; Delta_P = P_1 - P_2; P_ave = (P_1+P_2)/2; beta_ave = (beta_1+beta_2)/2.

### NPT maximum-likelihood objectives (constants dropped)

<!-- eq:16 -->
$$ \ln \frac{P(H|\beta_2, P)}{P(H|\beta_1, P)} = \beta_{\text{ave}}(\Delta G) - (\Delta\beta)H $$
- **what:** ML model for enthalpy test (constant T-pair, different P... here same P different T). Slope -Delta_beta on H.

<!-- eq:17 -->
$$ \ln \frac{P(V|\beta, P_2)}{P(V|\beta, P_1)} = \beta\big(\Delta G - (\Delta P)V\big) $$
- **what:** ML model for volume test (same T, different P). Slope -beta Delta_P on V.

<!-- eq:18 -->
$$ \ln \frac{P(E, V | \beta_2, P_2)}{P(E, V | \beta_1, P_1)} = \beta_{\text{ave}}(\Delta G) - (\Delta\beta)(E + P_{\text{ave}}V) - \beta_{\text{ave}}(\Delta P)V $$
- **what:** ML model for joint (E,V) test. Coefficient of E is -Delta_beta; coefficient of V is -(Delta_beta P_ave + beta_ave Delta_P).

## NPT harmonic-oscillator toy (volume-dependent spring constant)

Spring constant $K = (a/V)^2$, x_0 = 0.

<!-- eq:npt-toy-Delta -->
$$ \Delta(P,\beta) = (\beta P)^{-2} \sqrt{\frac{2\pi}{a^2 \beta}} $$
- **what:** isothermal-isobaric partition function for the toy (Delta = integral_V Q(V,beta) exp(-beta P V) dV).
- **symbols:** a - constant setting K=(a/V)^2.

<!-- eq:npt-toy-P -->
$$ P(x,V|\beta,P) = a(\beta P)^2 \sqrt{\frac{\beta}{2\pi}} \exp\left(-\frac{\beta a^2 x^2}{2V^2} - \beta PV\right) $$
- **what:** joint (x,V) distribution for the NPT toy.

Gibbs sampling of the toy: at fixed V, x is Gaussian with std $\sigma = V/(a\sqrt{\beta})$ (= $1/\sqrt{\beta K}$).
<!-- CHECK: paper writes "sigma = sqrt(K/beta) = (V/a) beta^{-1/2}"; the middle form sqrt(K/beta) is a typo - the RHS (V/a)beta^{-1/2} = 1/sqrt(beta K) is the correct Gaussian width. -->
Conditional V-sampling via rejection from an exponential envelope: $P(V|x_i) \propto \exp\!\left(-\beta\frac{a^2 x_i^2}{2V^2} - \beta PV\right)$, bounded by $M\exp(-\beta PV)$. Average acceptance efficiency factor $\exp(-\beta/2)$ (uses $\langle x^2\rangle = \beta V^2/a^2$).

## Temperature/pressure gap selection (rule of thumb)

<!-- eq:gap-T -->
$$ \sigma_E = T\sqrt{C_V k_B}, \qquad \frac{\Delta T}{T} = \sqrt{\frac{2k_B}{C_V}} = \frac{2 k_B T}{\sigma_E} $$
- **what:** choose temperatures whose energy-distribution means differ by ~2 sigma_E; optimal error discrimination near 1-2x (molecular systems ~1.5-2x) this gap.
- **symbols:** C_V - constant-volume heat capacity; sigma_E - std of the energy distribution.

<!-- eq:gap-P -->
$$ 2\sigma_V = \Delta P\left(\frac{\partial V}{\partial P}\right)_T, \quad \left(\frac{\partial V}{\partial P}\right)_T = -\frac{\sigma_V^2}{k_B T}, \quad |\Delta P| = \frac{2 k_B T}{\sigma_V} = \sqrt{\frac{2 k_B T}{V \kappa_T}} $$
- **what:** analogous pressure-gap rule for volume tests.
- **symbols:** sigma_V - std of volume distribution; kappa_T = -(1/V)(dV/dP)_T - isothermal compressibility.

Enthalpy-gap variant uses C_P: $\Delta T/T = \sqrt{2 k_B / C_P}$.

## Grand canonical (Appendix A)

<!-- eq:23 -->
$$ P(x,N|\beta,\mu) = \Xi(\beta,\mu)^{-1} \exp(-\beta E + \beta \mu N) $$
- **what:** grand canonical microstate probability; Xi is the grand partition function.

<!-- eq:25 -->
$$ \ln\frac{P_2(N|\beta,\mu_2)}{P(N|\beta,\mu_1)} = \beta\big(-[(PV)_2-(PV)_1] + [\mu_2-\mu_1]N\big) $$
- **what:** particle-number validation at fixed beta, two chemical potentials. Slope beta(mu_2-mu_1) on N.

<!-- eq:27 -->
$$ \ln \frac{P(N,E|\beta_2,\mu_2)}{P(N,E|\beta_1,\mu_1)} = -\big(\beta_2 (PV)_2 - \beta_1 (PV)_1\big) - [\beta_2 - \beta_1]E + [\beta_2\mu_2 - \beta_1\mu_1]N $$
- **what:** joint (N,E) grand-canonical validation.

<!-- eq:29 -->
$$ \ln \frac{P(\vec{N}, E | \beta_2, \vec{\mu}_2)}{P(\vec{N}, E | \beta_1, \vec{\mu}_1)} = -\big(\beta_2(PV)_2 - \beta_1(PV)_1\big) - [\beta_2 - \beta_1]E + [\beta_2\vec{\mu}_2 - \beta_1\vec{\mu}_1] \cdot \vec{N} $$
- **what:** multi-species generalization; dot product over species.

## Weighted least-squares on histogram ratios (Appendix B)

<!-- eq:30 -->
$$ \text{var}(\ln r_k) = \frac{\text{var}(p_{k,1})}{p_{k,1}^2} + \frac{\text{var}(p_{k,2})}{p_{k,2}^2} = \frac{1-p_{k,1}}{N_1 p_{k,1}} + \frac{1-p_{k,2}}{N_2 p_{k,2}} = \frac{1}{n_{k,1}} - \frac{1}{N_1} + \frac{1}{n_{k,2}} - \frac{1}{N_2} $$
- **what:** first-order variance of the log histogram-ratio in bin k (used to weight linear/ML fits).
- **symbols:** r_k = p_{k,2}/p_{k,1}; p_{k,i} = n_{k,i}/N_i - empirical bin probability; n_{k,i} - counts in bin k of sim i; N_i - total samples of sim i. Bin-probability variance estimator: var(p_k) = p_k(1-p_k)/N.

<!-- eq:31 -->
$$ \text{var}(r_k) = \left(\frac{n_{k,2} N_1}{n_{k,1} N_2^2}\right)^{2}\left(\frac{1}{n_{k,1}} - \frac{1}{N_1} + \frac{1}{n_{k,2}} - \frac{1}{N_2}\right) $$
- **what:** variance of the histogram ratio itself (for nonlinear error estimates).
<!-- CHECK: OCR gives prefactor (n_{k,2} N_1 / (n_{k,1} N^2))^2; interpreted as N_2^2 in the denominator (r_k = p_{k,2}/p_{k,1} = (n_{k,2}/N_2)/(n_{k,1}/N_1)). Verify indices against source. -->

<!-- eq:wls -->
$$ \vec{\alpha} = (X^T W X)^{-1} X^T W Y, \qquad \text{cov}(\vec{\alpha}) = (X^T W X)^{-1} $$
- **what:** weighted linear least squares estimate and its covariance. Y = data vector (log ratios), W diagonal with W_ii = 1/variance_i.

<!-- eq:nls-cov -->
$$ \text{cov}(\vec{\alpha}) = (J^T W J)^{-1}, \qquad J_{ij} = \frac{\partial f(y_i,\vec{\alpha})}{\partial \alpha_j} $$
- **what:** nonlinear least-squares parameter covariance; J is the model Jacobian at the minimum.

## Maximum-likelihood estimation (Appendix C)

<!-- eq:32 -->
$$ \ln \frac{P_2(\vec{X})}{P_1(\vec{X})} = -\vec{\alpha} \cdot \vec{X} $$
- **what:** general Boltzmann-ratio relation; alpha . X is shorthand for alpha_0 + sum_{j=1}^{M} alpha_j X_j (NOT a plain dot product).
- **symbols:** X - vector of M sample variables (E, V, N...); alpha_j - conjugate intensive variables; M+1 parameters total.
<!-- CHECK: source prints "= exp(-alpha.X)" but from eqs 6/32-context the log-ratio equals -alpha.X (the exp belongs to the ratio, not its log). Corrected to -alpha.X. -->

<!-- eq:33 -->
$$ \frac{P(\vec{X}|2)}{P(\vec{X}|1)} = \frac{P(2|\vec{X})}{1 - P(2|\vec{X})}\,\frac{P(1)}{P(2)}, \qquad \frac{P(1)}{P(2)} = \frac{N_1}{N_2} $$
- **what:** Bayes rewrite turning the density ratio into a logistic classification of which simulation a sample came from.

<!-- eq:34 -->
$$ \frac{P(2|\vec{X})}{1 - P(2|\vec{X})} = \exp\left(-M - \alpha_0 - \sum_{j=1}^{M} \alpha_j X_j\right), \qquad M = \ln(N_F/N_R) $$
- **what:** logistic odds in terms of the parameters; M = ln(N_1/N_2) offset from unequal sample counts.

<!-- eq:35 -->
$$ P(1|\vec{X}_i) = \frac{1}{1 + \exp(M + \vec{\alpha} \cdot \vec{X}_i)}, \qquad P(2|\vec{X}_i) = \frac{1}{1 + \exp(-M - \vec{\alpha} \cdot \vec{X}_i)} $$
- **what:** per-sample class probabilities (logistic).

<!-- eq:lnL-general -->
$$ \ln L(\vec{\alpha}|\text{data}) = \sum_{i=1}^{N_1} \ln f(-M - \vec{\alpha} \cdot \vec{X}_i) + \sum_{i=1}^{N_2} \ln f(M + \vec{\alpha} \cdot \vec{X}_i) $$
- **what:** general multivariate log-likelihood; single minimum so the gradient root is unique.
- **symbols:** f(x) = [1 + exp(x)]^{-1} - Fermi function (appendix sign convention; note eq:8 wrote f(x)=[1+exp(-x)]^{-1} with flipped argument, consistent overall).

<!-- eq:36 -->
$$ \text{var}(\alpha_j) = I(\alpha_j)^{-1} = -\frac{1}{\dfrac{\partial^2 \ln L(\alpha)}{\partial \alpha_j^2}} $$
- **what:** ML parameter variance from the Fisher information (negative inverse Hessian of log-likelihood). No extra sample-count correction is needed here (unlike Shirts 2003, since two parameters are fit from two distributions with no implicit beta constraint).
