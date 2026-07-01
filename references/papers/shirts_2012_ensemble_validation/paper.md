# Simple quantitative tests to validate sampling from thermodynamic ensembles

**Michael R. Shirts**, Department of Chemical Engineering, University of Virginia (arXiv:1208.0910).

## Abstract

It is often difficult to quantitatively determine if a new molecular simulation
algorithm or software properly implements sampling of the desired thermodynamic
ensemble. This paper presents simple statistical procedures to sensitively
determine whether a desired thermodynamic ensemble is properly sampled, and
demonstrates them for model systems and MD simulations across constant-volume
and constant-pressure conditions, with an implementation for end users
(`checkensemble`, simtk.org/home/checkensemble).

## 1. Introduction

Molecular simulations (MD and MC) generate statistical samples from a target
ensemble. Subtle errors arise from numerical errors in energy functions,
theoretical errors in algorithms, over-aggressive approximations, and ordinary
programming bugs. Validation is hard because analytical results are unavailable
for complex systems, and observables need long simulations to resolve small
violations.

Some aspects can always be checked directly: NVE total energy conservation
(RMS error scales with step size squared for symplectic integrators); NVT
kinetic energy following the Maxwell-Boltzmann distribution; NPT average
instantaneous pressure from virial + kinetic energy. But there is **no standard
test for the potential-energy distribution** (critical for MC) or for the total
energy of an arbitrary system, and many distributions share the correct average
temperature/pressure without satisfying the correct Boltzmann distribution.
The goal is a physically rigorous test that a method samples the target
distribution in its entirety.

## 2. Theory

Thermodynamic ensembles share the exponential form `exp(-u(x))` where `x` is the
microstate and `u(x)` a reduced energy term (see eq:1-3 for canonical, NPT,
grand canonical). A general test checks directly that collected samples are
consistent with these distributions, plus the requirement that all microstates
with the same generalized energy have equal probability.

The canonical energy density can be written via the density of states Omega(E)
(eq:4). Without knowing Omega(E), no single-state sample set can confirm the
distribution. But the **ratio of two distributions at different beta** (same
other parameters) cancels the unknown Omega(E) (eq:5), and its logarithm is
**linear in E** (eq:6) with slope `alpha_1 = -(beta_2 - beta_1)` independent of
the free energies. This is the basis of the *ensemble consistency* test.
Analogous relations hold for any ensemble whose reduced energy is linear in its
conjugate parameters.

This is a **necessary but not sufficient** test: a deviating slope proves
non-canonical sampling, but agreement does not prove ergodicity or that all
equal-energy states are equally sampled, nor that all phase space is reached.
The paper uses systems known to sample ergodically so no separate
convergence/ergodicity tests are needed.

### 2.1 Visual inspection

Divide the common energy range into ~20-40 bins (need not be equally spaced;
aligning bins between the two data sets simplifies analysis; exclude extreme
bins to avoid small-sample error and zero densities). Estimate P_1(E), P_2(E)
per bin, plot log of the ratio. If linear, the system qualitatively obeys the
equilibrium distribution. A slope below expected means the high-temperature
simulation undersamples that energy (distribution too wide); above means
oversampling (too narrow).

**Fig 1 (semantics):** ensemble validation of 900 TIP3P water in NVT with
Nose-Hoover; linear and nonlinear fits agree with prediction.

### 2.2 Quantitative fitting

Visual inspection of paired replica-exchange probability ratios has been used
before, but not as a general quantitative goodness-of-fit test.

**Linear fitting (2.2.1):** estimate the bin-occupancy error as
`delta p_k = sqrt(p_k(1-p_k)/n)` (standard result) and propagate into the ratio
P_2/P_1 (Appendix B). Enables weighted linear/nonlinear least squares including
small-sample tail effects.

**Nonlinear fitting (2.2.2):** linearization biases parameter estimates, so
minimize the direct residual sum (eq:7) instead, then propagate bin errors into
alpha.

**Maximum likelihood (2.2.3):** eliminates histogram dependence entirely and
includes the tails. Two fit parameters, `alpha_0 = ln(Q_1/Q_2)` and
`alpha_1 = -(beta_2 - beta_1)`, via a logistic log-likelihood (eq:8) with the
Fermi function `f(x) = [1+exp(-x)]^-1`; the first sum is over T_1 energies, the
second over T_2 energies. Everywhere concave with a single maximum -> always
converges. There are nominally four parameters (A_1, A_2, beta_1, beta_2) but
only two are free: fix `beta_ave = (beta_{1,user}+beta_{2,user})/2` as an
energy-scale choice and set A_1+A_0=0 so DeltaA = A_2-A_1 is the second free
variable, giving `alpha_0 + alpha_1 E = beta_ave*DeltaA - Delta_beta*E`. ML gives
no plottable histogram, so run a linear fit alongside for visualization.

### 2.3 Error estimates

Asymptotic covariance estimators for the fit parameters are standard (eq:cov):
`(X^T W X)^-1` (linear), `(J^T W J)^-1` (nonlinear), and inverse Hessian of the
log-likelihood (ML). Bootstrap sampling of the original data is also reliable
(200 bootstrap samples took ~20 min on one core of a 2.7 GHz i7 for 600k
energies per simulation). With enough samples the fit deviates normally, so
compute how many standard deviations the fitted slope beta_2-beta_1 lies from the
user-specified value; consistently > 2-3 sigma across repeats indicates simulation
error. **The number of standard deviations measures certainty of an error, not
its magnitude** (high-precision detection of a tiny error, or low-precision of a
large one).

### 2.4 Choosing the parameter gap

Eq:6 holds for any beta_1, beta_2, but if they are too far apart the two
distributions never overlap with well-defined probabilities; if equal, no
information. There is an ideal intermediate gap: distributions whose mean energy
separation approximately equals the sum of their standard deviations (~2 sigma)
maximize error discrimination. Gap rules (eq:gap-T): `DeltaT/T = sqrt(2 k_B/C_V)`
or `= 2 k_B T/sigma_E` with `sigma_E = T sqrt(C_V k_B)`. Err toward a smaller
gap when data is limited (to keep overlap near E=0). The 1-D harmonic oscillator
(constant Omega(E)) is an exception where sensitivity always increases with gap.

### 2.5 Harmonic-oscillator toy (canonical)

D-dimensional harmonic oscillator, equal spring K, equilibrium x_{i,0}
(eq:ho-energy, eq:ho-Q, eq:ho-P). Use `D=20` (not D=1, whose constant density of
states is atypical: E=0 has nonzero probability at all temperatures). Average
energy `<E> = D/(2 beta)`.

**2.5.1 Noise-injection test:** sample with K=1, beta=1.3 and 0.7, add noise
`dE = nu*|N(0,1)|`, 200 repetitions of 500,000 samples each. True beta_2-beta_1
= 0.6; average energy D/(2 beta) = 16.667. Compare linear, nonlinear, ML fits
with analytic, replicate, and bootstrap errors (Table 1). Bootstrap matches
independent-replicate error; analytic errors match sample std for linear and ML;
nonlinear analytic error underestimates. Deviations exceed 3 sigma consistently
for nu >= 0.0075 (< 1% of k_B T), errors invisible to the eye (Fig 2). ML is the
best overall (no discretization error); linear is also robust.

**Gap optimization (Table 2):** fix nu=0.01, (beta_1+beta_2)/2=1, vary the gap.
Error discrimination peaks at intermediate separation. At the peak
(beta_1=0.7, beta_2=1.3) the center gap 6.6 k_B T ~ sum of std devs 6.9 k_B T,
supporting the "separate by ~sum of std devs" rule. Estimate sigma via
`sigma_E = T sqrt(C_V k_B)` or a short center-temperature run.

### 2.6 Isobaric-isothermal (NPT)

Same principles apply to constant-T,P and constant-T,mu simulations. Three
useful NPT tests: (a) two sims same P, different T -> enthalpy H=E+PV validation
(eq:9, eq:10) with Gibbs free energy G replacing A plus a ln(beta_1/beta_2)
correction; (b) same T, different P -> volume validation (eq:11), slope
-beta(P_2-P_1) on V; (c) both T and P differ -> joint (E,V) multilinear fit
(eq:14). Multidimensional histograms populate poorly, so the joint case uses ML
only (histogram-free). Free-variable choices: DeltaG=G_2-G_1 (set G_1+G_2=0),
Delta_beta=beta_2-beta_1 (fix beta_ave), Delta_P=P_1-P_2 (fix P_ave); the
V-coefficient decomposes as eq:15. ML objectives: enthalpy eq:16, volume eq:17,
joint eq:18.

### 2.7 NPT harmonic-oscillator toy

Modified HO with volume-dependent spring `K=(a/V)^2`, x_0=0, adding a PV work
term (eq:npt-toy-Delta, eq:npt-toy-P). Sample the joint P(E,V) with a Gibbs
sampler alternating P(E|V) and P(V|E). At fixed V, x is Gaussian with
`sigma = V/(a sqrt(beta))`. Sample V by rejection from an exponential envelope
`M exp(-beta P V)`; average efficiency factor `exp(-beta/2)`, independent of P
and a.

**2.7.1 Results (250,000 samples each, ML):** enthalpy (beta_1=2/3, beta_2=2,
P=1) -> beta_2-beta_1 = 1.3341 +/- 0.0040 vs true 4/3; volume (beta=1,
P_1=1.3, P_2=0.7) -> beta(P_2-P_1) = -0.6013 +/- 0.0025 vs true -0.6; joint
(beta_1=0.6, beta_2=0.8, P_1=0.8, P_2=1.2) -> slopes 0.20035 +/- 0.00318 (true
0.2) and -0.48129 +/- 0.00185 (true -0.48). See checks.md.

**2.7.2 NPT gap selection:** enthalpy gap uses C_P: `DeltaT/T = sqrt(2 k_B/C_P)`;
volume gap `|DeltaP| = 2 k_B T/sigma_V = sqrt(2 k_B T/(V kappa_T))` with
`(dV/dP)_T = -sigma_V^2/(k_B T)` and isothermal compressibility
`kappa_T = -(1/V)(dV/dP)_T` (eq:gap-P).

## 3. Molecular systems

### 3.1 Kinetic and potential energy independently obey the test

When potential energy is independent of momenta (typical), the total-energy
distribution factorizes into kinetic and potential parts (eq:sep), so the test
applies to E_kin and E_pot separately (eq:19, eq:20). E_pot validation works for
MC (only potential energies defined). The explicit kinetic ratio (eq:21, eq:22)
depends only on Delta_beta and masses cancel. The kinetic energy is a sum of 3N
squared normal variables -> chi-squared with 3N DOF (minus removed COM DOF),
essentially Gaussian above ~60 DOF (~20 particles) with mean 3N/(2 beta) and
variance 3N/(2 beta^2) (eq:kin-dist). It can be checked directly (Q-Q plots,
Anderson-Darling) and via equipartition `<E_kin> = (k_B T/2) * #DOF`. NPT
separation is analogous (Delta = Q_kin * Delta_pot), enabling NPT MC validation
by removing kinetic energy.

### 3.2 Lennard-Jones argon MD

300 LJ particles, Gromacs 4.6 double precision, Rowley-Nicholson-Parsonage argon
(sigma=0.3405 nm, epsilon=119.8 K, k_B=0.996072 kJ/mol), rho=0.85 rho_c (box
3.5328256 nm), T=0.85 T_c=135.0226 K, velocity Verlet, COM momentum removed each
step, dispersion correction, LJ switch 0.8-0.9 nm, neighborlist 1.0 nm. Default
T_low=132.915, T_high=137.138 K (~0.7x ideal gap, C_V~8.5 kJ/(mol K)). See
checks.md for reduced-time / stability details.

### 3.3 Temperature-control algorithms (Table 3)

All tested thermostats (Bussi-Parrinello, Andersen, Andersen-massive,
Nose-Hoover, stochastic dynamics) are canonical-consistent within ~1 sigma
**except Berendsen**, which gives an overly narrow kinetic-energy distribution
(68 sigma off). NVE kinetic energy deviates but NVE potential energy does not.
DeltaT is derived from Delta_beta via `T_{2,1} = k_B^{-1}(beta_ave +/-
Delta_beta/2)^-1`. Correlation times tau computed with pymbar timeseries;
subsample at 2 tau + 1 using potential-energy correlation times (kinetic
correlation times are artificially short).

**Fig 4 (semantics):** Berendsen kinetic-energy log-ratio slope is 7x too high
(68 sigma); Nose-Hoover is statistically indistinguishable from correct.

### 3.3.1 Large step size (Table 4)

Long steps heat NVE runs; a thermostat hides this as a steady state that need
not be Boltzmann. Bussi-Parrinello, 8-40 fs (48 fs segfaults). Total and kinetic
energy deviate as dt grows, worst near the instability point; potential energy
deviates much less. The full-step velocity-Verlet KE estimator is less accurate
than the half-step-averaged leapfrog KE estimator; at 40 fs the leapfrog KE
gives a markedly better kinetic distribution without changing potential.

### 3.3.2 Abrupt cutoff (Table 5)

Abrupt radial cutoffs create force discontinuities and NVE heating. Cutoffs
r_c=2.0-4.0 sigma with Bussi-Parrinello. Distributions are surprisingly
ensemble-consistent for most cutoffs; only r_c < 3 sigma show clear violations.
Here the average-KE temperature deviation (23 sigma at r_c=2) is *more* sensitive
than the ensemble test (opposite of the step-size case) - use multiple methods.

### 3.3.3 Gap selection for molecular systems (Table 6)

LJ argon, dt=32 fs, C_V=8.5 kJ/(mol K) at 135 K -> sigma_E=36 kJ/mol, ideal gap
~6 K. Maximum error sensitivity for total energy is between 1x and 2x the
estimated gap; a slightly larger gap can help. Kinetic-energy sensitivity peaks
in the 1-2x range using kinetic std. Recommend ~1.5-2x the estimated gap for
molecular systems.

### 3.4 Pressure-control algorithms (Table 7)

Gromacs barostats: Berendsen, Parrinello-Rahman (PR), MTTK. Argon, P=90 bar,
T=125 K, tau_p=5 ps, dt=8 fs. Volume/joint tests DeltaP=120 bar (Berendsen 4
bar); enthalpy tests DeltaT=7.138 K (C_P~10.2 kJ/mol). PR and MTTK give accurate
enthalpy but volume DeltaP off by ~5 bar (~5%, 7-9 sigma). Berendsen fails all
three (needs a much narrower parameter range to overlap; unphysically narrow
volume distributions; anomalously long volume autocorrelation ~110-130 ps).

### 3.5 Water (Tables 8-9)

900 TIP3P, velocity Verlet, PME (cutoff 1.0 nm, order 6, tol 1e-6), LJ switch
0.8-0.9 nm, SETTLE constraints, dt=2 fs, 20 ns. NVT: T=298/301 K, DeltaT=3 K
(slope 0.004023 (k_B T)^-1). All thermostats consistent except Berendsen (wildly
off). NPT: PR and MTTK reasonable for enthalpy; MTTK volume distribution better
than PR (PR pressure lags one time step with leapfrog). Both far better than
Berendsen. See checks.md.

## 4. Tools

`checkensemble` (simtk.org/home/checkensemble): automatic linear/nonlinear
plotting, and linear, nonlinear, and maximum-likelihood parameter analysis, with
example parsers for Gromacs, CHARMM, and Desmond output for NVT and NPT
(enthalpy, volume, joint E-V). Requires only lists of the relevant extensive
variables and the intensive applied variables.

## 5. Conclusions

For Boltzmann-distributed systems (essentially all NVT/NPT molecular
simulations), ensemble consistency can be checked from **pairs of simulations**
differing in one external parameter (T, P, or mu): system-dependent densities of
states cancel, leaving a **linear** relationship between the log probability
ratio and an extensive variable (energy, volume, enthalpy, particle number),
whose slope is fixed entirely by the user's intensive variables. With proper
error analysis this gives quantitative validation. All tested thermostats except
Berendsen pass; among barostats MTTK is best, Parrinello-Rahman acceptable,
Berendsen wrong when volume fluctuations matter. Optimal parameter gaps separate
distribution means by 2-4x the sum of the standard deviations (err smaller with
less data). Underestimating autocorrelation time is the main source of spurious
deviations. These are sensitive **necessary but not sufficient** tests (they do
not guarantee equal sampling of equal-energy states or full phase-space
coverage). Future work: grand canonical support and full hypothesis testing.

---

## Derivation (not implemented): NPT enthalpy distribution (eq:9-10)

Starting from eq:2 with microstate probability
`P(x,V|beta,P) = Delta(beta,P)^-1 exp(-beta E(x) - beta P V)`, integrate out
configurations at fixed instantaneous enthalpy H=E+PV using a delta function:
`P(H|beta,P) = (beta P/h^{3N}) integral_V integral_x delta[E(x)+PV-H]
Delta^-1 exp(-beta(E+PV)) dx dV = (beta P/h^{3N}) Omega'(H,P) Delta^-1 exp(-beta H)`
(eq:npt-H), where Omega'(H,P) counts states with E+PV=H and depends on P but not
beta. The beta*P prefactor cancels units. Taking the ratio at two temperatures,
same P, gives eq:9; its log gives eq:10. Volume-only (eq:11) and joint-(E,V)
(eq:14) relations follow analogously by integrating out E at fixed V, or keeping
both.

## Derivation (not implemented): kinetic/potential separability (eq:sep, 19-22)

Because momenta sample independently of coordinates,
`P(E_pot+E_kin|beta) = Q_kin^-1 Q_pot^-1 Omega(E_pot) Omega(E_kin) exp(-beta E_pot)
exp(-beta E_kin) = P(E_pot|beta) P(E_kin|beta)` (eq:sep). With
`Q_kin = prod_i integral exp(-beta p_i^2/m_i) dp_i = prod_i (m_i/(pi beta))^{3/2}`
the kinetic ratio (eq:21) has all mass terms cancel, yielding eq:22 in the single
parameter Delta_beta. The Gaussian kinetic distribution (eq:kin-dist) follows
from the chi-squared sum with mean 3N/(2 beta) and variance sigma^2 = k_B T^2 C_V
= 3N/(2 beta^2) using the ideal-gas kinetic heat capacity C_V = 3 N k_B/2.

## Derivation (not implemented): grand canonical relations (Appendix A, eq:23-29)

All GC equations mirror the NPT case with -mu replacing P and N replacing V.
From eq:23, examining P(N) at fixed beta, two chemical potentials, gives eq:25;
joint (N,E) gives eq:27; multi-species gives eq:29. Since N is discrete,
histogramming introduces no extra approximation if bins are integer-fine.

## Derivation (not implemented): weighted least squares (Appendix B, eq:30-31)

For a histogram with K bins, N samples, counts n_k, empirical probability
p_k = n_k/N has variance p_k(1-p_k)/N (standard result). First-order propagation
of the log-ratio ln r_k = ln(p_{k,2}/p_{k,1}) gives eq:30; the ratio variance
gives eq:31. For a linear model F=AY, cov(F)=A cov(Y) A^T; with cov(Y)=W
(inverse-variance weights), weighted least squares gives
alpha=(X^T W X)^-1 X^T W Y and cov(alpha)=(X^T W X)^-1 (eq:wls); the nonlinear
case replaces X by the Jacobian J (eq:nls-cov).

## Derivation (not implemented): maximum likelihood (Appendix C, eq:32-36)

The Boltzmann ratio (eq:32) is reinterpreted via Bayes as a logistic
classification of which simulation a sample came from (eq:33), with
P(1)/P(2)=N_1/N_2. Defining M=ln(N_1/N_2) gives the logistic odds (eq:34) and
per-sample class probabilities (eq:35), whose product yields the log-likelihood
(eq:lnL-general). The Fisher information gives the parameter variance (eq:36);
no extra sample-count correction is needed (unlike Shirts 2003) because two
parameters are fit from two distributions with no implicit beta constraint.

## Exercises / limitations noted
- These tests are necessary but not sufficient: they cannot verify equal
  sampling of equal-energy states, nor complete phase-space coverage. Combine
  with ergodicity/convergence checks and simple average-temperature/pressure
  diagnostics.
