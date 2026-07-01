# Identifying the Optimal Integration Time in Hamiltonian Monte Carlo

Michael Betancourt (2016). arXiv:1601.00225v1.

## Abstract

By leveraging the natural geometry of a smooth probabilistic system, Hamiltonian
Monte Carlo (HMC) yields computationally efficient Markov Chain Monte Carlo
estimation, provided the algorithm is sufficiently well-tuned. This paper shows
how the geometric foundations of HMC implicitly identify the optimal choice of
these parameters, especially the integration time. The practical consequences
are considered in existing algorithms and in a new implementation called
Exhaustive Hamiltonian Monte Carlo (XHMC), with illustrative examples.

Keywords: Markov Chain Monte Carlo, Hamiltonian Monte Carlo, Microcanonical Systems.

---

## 1. Hamiltonian Monte Carlo in Theory

The goal is to compute expectations of a function $f: Q \to \mathbb{R}$ with
respect to a target distribution $\pi$ on an $N$-dimensional sample space $Q$,
using MCMC estimators. HMC builds a Markov kernel by mapping the probabilistic
system into a Hamiltonian system whose canonical measure-preserving flow
generates a powerful Markov transition. Performance depends crucially on **how
long** the flow is integrated: too short and the chain devolves into diffusive
exploration; too long and integration is wasteful.

### 1.1 Constructing a Generic Hamiltonian Kernel

Consider the smooth probabilistic system $(Q, \mathcal{B}(Q), \pi)$, with $Q$ a
positively-oriented smooth $N$-dimensional manifold, $\mathcal{B}(Q)$ the
canonical Borel $\sigma$-algebra, and $\pi$ a smooth probability distribution.

Any choice of a disintegration $\xi$ on the cotangent bundle
$\varpi: T^*Q \to Q$ lifts the target onto the cotangent bundle:

$$\pi_H = \varpi^* \pi \wedge \xi.$$

Denoting $\theta$ the tautological one-form and $\Omega = \wedge_{n=1}^N d\theta$
the symplectic volume form, the Hamiltonian is

$$H = -\log \frac{\mathrm{d}(\varpi^* \pi \wedge \xi)}{\mathrm{d}\Omega},$$

giving the Hamiltonian system $(T^*Q, \mathrm{d}\theta, H)$. The lifted target is
the canonical measure on the cotangent bundle,

$$\pi_H = e^{-H}\Omega,$$

which is preserved by the canonical Hamiltonian flow.

Locally, the target decomposes as
$\pi = e^{-V}\, dq^1 \wedge \ldots \wedge dq^n$, where $V$ is the **potential
energy**; the disintegration decomposes as
$\xi = e^{-K}\, dp_1 \wedge \ldots \wedge dp_n + \text{horizontal } n\text{-forms}$,
where $K$ is the **kinetic energy**. The lift then gives the classical form

$$H = -\log \frac{\mathrm{d}\pi_H}{\mathrm{d}\Omega} = K + V.$$

A Hamiltonian Markov transition (a) lifts an initial point from the sample space
to the cotangent bundle by sampling momentum from the fiber, (b) applies the
Hamiltonian flow $\phi_t^H$ for a random time $t \sim \pi_{T(q,p)}$, and (c)
projects back to the sample space. Composing $g = \varpi \circ \phi_t^H \circ l$
yields measure-preserving diffeomorphisms and defines the Hamiltonian kernel

$$\mathcal{T}_{\mathrm{HMC}}(q,A) \equiv \int_G \gamma_q(\mathrm{d}g)\, \mathbb{I}_A(g(q)),$$

with $\mathbb{I}_A$ the indicator function.

### 1.2 Optimal Kernel from the Geometry of Microcanonical Systems

This construction is too general: every choice of cotangent disintegration
$\xi$ and integration-time distribution $\pi_{T(q,p)}$ yields a different
kernel with substantially varying performance.

**1.2.1 The Microcanonical Disintegration.** Hamiltonian systems foliate into
level sets of constant energy,

$$H^{-1}(E) = \{ z \in T^*Q \mid H(z) = E \}.$$

The canonical distribution restricted to each foliated component disintegrates
into microcanonical distributions uniform on each level set, with a marginal
energy distribution $\pi_E = H_* \pi$. Any expectation decouples into a
microcanonical expectation nested in an expectation over energies (Eq. 1).

The microcanonical disintegration is **compatible with the Hamiltonian flow**:
every trajectory is confined to a single level set, and, because the flow
restricted to a level set preserves the microcanonical distribution, long
trajectories explore that distribution. A Hamiltonian Markov chain thus
decouples into a deterministic flow along level sets plus a momentum resampling
that induces a random walk **between** level sets.

The autocorrelation depends on (i) how effectively the flow explores each
microcanonical distribution (controlled by integration time), and (ii) how
effectively momentum resampling explores the marginal energy distribution
(controlled by $\Delta H$ relative to the width of $\pi_E$). Note on Fig 4: when
the expected energy variation $\Delta H$ per resampling matches the width of
$\pi_E$, exploration is rapid; when $\Delta H$ is small, autocorrelations are
large. This provides a diagnostic for a poorly chosen cotangent disintegration.

**1.2.2 Ergodicity of Hamiltonian Flow.** Under dynamical ergodicity, the
Birkhoff ergodic theorem gives temporal average = spatial (microcanonical)
average as $T \to \infty$, motivating $\pi_{T(z)} = U(0, T(z))$. Hamiltonian
systems are not always ergodic; the only generic guarantee is that the time
average converges to the spatial expectation over the **orbit**
$\phi^H(z) = \{\phi_t^H(z), \forall t \in \mathbb{R}\} \subset H^{-1}(H(z))$.

Convergence is initially superlinear (justifying the linear cost of longer
integration) but for long times converges only as $\sqrt{T}$ (diminishing
returns). Optimal performance requires a **maximal integration time** $T(z)$
with $\pi_{T(z)} = U(0, T(z))$ that identifies the transition between the two
regimes uniformly across all level sets — intuitively, just after the
trajectory first traverses the extent of its orbit.

**1.2.3 Poincaré Recurrence and Autocorrelation Functions.** If the Hamiltonian
is proper with compact level sets, all orbits are bounded and Poincaré
recurrence guarantees trajectories return to any neighborhood within a finite
recurrence time. A practical alternative is an auxiliary autocorrelation
function $\kappa(T,z)$ that converges monotonically to zero, relaxed to a
uniform-bound termination criterion (Eq. 2), defining integration times
$T_\kappa(z) = \min\{ t \mid |\kappa(t,z)| \le \delta \}$.

---

## 2. Hamiltonian Monte Carlo in Practice

Implementing the flow requires solving $2n$ first-order ODEs, done numerically
with **symplectic integrators** that preserve the symplectic volume with only
small variations in $H$. By backwards error analysis, a $k$-th order symmetric
symplectic integrator with step size $\epsilon$ exactly follows the flow of a
**modified Hamiltonian**

$$\widetilde{H} = H + \sum_{n=k/2}^{N} \epsilon^{2n} H_{(n)} + \mathcal{O}\!\left(e^{-c/\epsilon}\right),$$

with leading-order behavior $\widetilde{H} = H + \epsilon^k G + \mathcal{O}(\epsilon^{k+2})$.
The discretized flow generates states
$z_L \equiv \Phi^{\widetilde{H}}_{\epsilon, L\cdot\epsilon}(z_0)$, $L \in \mathbb{Z}$.
The integrator introduces error that biases the chain if uncorrected.

### 2.1 Static Implementations

The simplest scheme uses a static integration time $T(z) = T$, i.e. a fixed
number of leapfrog steps $L = T/\epsilon$. Using only the final point (a Dirac
measure $\pi_{T(z)} = \delta_{L\cdot\epsilon}$) as a Metropolis proposal is
**invalid** because the flow is non-reversible; it becomes valid only when
composed with an involution operator $R$ (e.g. momentum flip) satisfying
$\Phi^{\tilde{H}}_{\epsilon,L\cdot\epsilon}\circ R\circ\Phi^{\tilde{H}}_{\epsilon,L\cdot\epsilon} = \mathrm{Id}$.
Acceptance probability then uses $H\circ R(z_L) - H(z_0)$.

Sampling **uniformly from the entire trajectory** $\pi_T = U(0,T)$ (as motivated
by the microcanonical geometry) breaks detailed balance if states are generated
by forward-only integration. To restore detailed balance, one considers the set
$\mathfrak{T}_{z,L}$ of all length-$L$ trajectories containing $z$: first sample
a trajectory $\mathfrak{t} \in \mathfrak{T}_{z_0,L}$ with probability
$\mathbb{P}[\mathfrak{t}|z_0]$, then sample a state $z$ from it with the
Metropolis probabilities (Eq. 3). Detailed balance is guaranteed if the
trajectory sampling probability is equal across all initial states in the
trajectory (Eq. 4).

One immediate solution: sample trajectories in $\mathfrak{T}_{z_0,L}$ uniformly,
$\mathbb{P}[\mathfrak{t}|z_0] = 1/L$. This is implemented by sampling
$L' \sim U[0,L]$ and integrating backwards $L'$ steps and forwards $L-L'$ steps.
Combined with direct Metropolis sampling of the final state, this equals Neal's
Windowed State Algorithm with $W=L$.

### Derivation (not implemented): detailed balance for trajectory sampling

Given the Metropolis weighting (Eq. 3) and the equal-trajectory-probability
condition (Eq. 4), detailed balance follows:

$$\mathbb{P}[z_{1}|z_{2}] \frac{\mathrm{d}\pi_{H}}{\mathrm{d}\Omega}(z_{2})
= \sum_{\mathfrak{t} \in \mathfrak{T}_{(z_{1},z_{2}),L}} \mathbb{P}[z_{1}|\mathfrak{t}]\, \mathbb{P}[\mathfrak{t}|z_{2}] \frac{\mathrm{d}\pi_{H}}{\mathrm{d}\Omega}(z_{2})
= \ldots = \mathbb{P}[z_{2}|z_{1}] \frac{\mathrm{d}\pi_{H}}{\mathrm{d}\Omega}(z_{1}),$$

using $\mathbb{P}[\mathfrak{t}|z_1] = \mathbb{P}[\mathfrak{t}|z_2]$ on the shared
support $\mathfrak{T}_{(z_1,z_2),L}$.

Alternative state-sampling from Eq. 3: sample from the multinomial defined by the
Metropolis probabilities directly, or use a slice sampler ($u \sim U(0,1)$, then
uniformly sample points with $e^{-H(z)}/\sum_{z'\in\mathfrak{t}} e^{-H(z')} > u$).

### 2.2 Dynamic Implementations

To maintain a uniform distribution over trajectories of **dynamic** length $L$,
build the trajectory incrementally, checking a termination criterion after each
expansion (Algorithm 1). Two schemes:

- **Additive:** iteratively expand one step in a random direction. Requires $L$
  termination checks after a length-$L$ proposal and storage of every state.
- **Multiplicative:** only lengths $2^D$; expand a length-$L$ trajectory by
  integrating $L$ additional steps in a random direction. Represented as a
  balanced binary tree; requires only $\log(L)$ checks and $\log(L)$ states in
  memory.

Uniform sampling alone is **not** sufficient for detailed balance when length is
dynamic, because different initial states may reach different terminal lengths.
Each increment must be treated as a proposal, rejecting any extension containing
states from which $\mathbb{P}[\mathfrak{t}_{new}|z] = 0$ (Algorithm 2). For
multiplicative expansion, validation requires that **no internal subtree** of the
proposal satisfies the termination criterion.

Algorithm 1 (naive build):
```
function EXPAND_TRAJECTORY(t)
function CHECK_TERMINATION(t)
function NAIVE_BUILD_TRAJECTORY(t)
   t_new <- EXPAND_TRAJECTORY(t)
   if CHECK_TERMINATION(t_new) then
      return t_new
   else
      NAIVE_BUILD_TRAJECTORY(t_new)
```

Algorithm 2 (build with validation for detailed balance):
```
function EXPAND_TRAJECTORY(t)
function VALIDATE_TRAJECTORY(t)
function CHECK_TERMINATION(t)
function BUILD_TRAJECTORY(t)
   t_new <- EXPAND_TRAJECTORY(t)
   if VALIDATE_TRAJECTORY(t_new) then
      if CHECK_TERMINATION(t_new) then
         return t_new
      else
         NAIVE_BUILD_TRAJECTORY(t_new)
   else
      return t
```

### 2.3 Alternative Schemes

Poor performance stems from the momentum resampling induced by projecting to and
from the cotangent bundle. Constructing a chain on the cotangent bundle directly
can avoid this.

- **Horowitz scheme (1991):** Hamiltonian flow with only **partial** momentum
  resampling. After the prescribed integration time a Metropolis correction is
  applied; on acceptance the final momentum is mixed with new momenta
  (maintaining coherence); on rejection the momentum must be **completely
  negated**, returning to an already-explored neighborhood.
- **Extra-chance schemes (Sohl-Dickstein et al. 2014; Campos & Sanz-Serna
  2015):** apply a fixed number of proposals that do not modify momentum after an
  acceptance, still negating after rejections.

Both maintain coherence only while the symplectic integrator is near the true
flow, devolving to diffusion when the integrator strays and proposals are
rejected. Because symplectic integrators are volume-preserving they do not drift;
excursions are only temporary. To fully exploit the flow one must integrate past
these temporary excursions.

---

## 3. Explicit Termination Criteria

### 3.1 Theoretical Exhaustions

A natural autocorrelation function is the temporal expectation of the temporal
derivative of a scalar function $u$:

$$\kappa_u(T, z) \equiv \frac{1}{T} \int_0^T \mathrm{d}t\, \frac{\mathrm{d}u}{\mathrm{d}t} \circ \phi_t^H(z) = \frac{u \circ \phi_T^H(z) - u(z)}{T}.$$

If $u$ is bounded along the flow, $\kappa_u$ vanishes asymptotically, making it a
potential termination criterion (care: $u$ may recur, breaking monotonicity).

The Hamiltonian itself is conserved (rate of change vanishes trivially, so it is
unsuitable). The other canonical scalar is the **virial** $G = q^i p_i$. When $H$
is proper and trajectories bounded, $G$ is bounded and is a candidate.

**Definition 1 (Exhaustion).** An exhaustion $T_\delta(z)$ is the family of
integration times such that the temporal average of the rate of change of the
virial is uniformly bounded (Eq. exhaustion). A valid exhaustion always exists
for a proper Hamiltonian. Exhaustions reduce tuning to the single threshold
$\delta$; how to choose an optimal $\delta$ remains open.

### 3.2 Numerical Exhaustions

Replace the exact temporal expectation with a Metropolis-corrected expectation
over the numerical trajectory (Eq. numexp). This converges to the continuous
expectation as $|\mathfrak{t}| \to \infty$.

**Definition 2 (Numerical Exhaustion).** A numerical exhaustion $\mathfrak{T}_\delta$
is the set of numerical trajectories whose Metropolis-corrected expectation of
the virial rate is uniformly bounded by $\delta$. The trajectory termination
criterion checks $\mathfrak{t} \in \mathfrak{T}_\delta$; the resulting algorithm
is **Exhaustive Hamiltonian Monte Carlo (XHMC)**.

Divergences of the numerical trajectory serve as diagnostics of an ill-posed
numerical exhaustion.

### 3.3 Riemannian Termination Criteria

A Riemannian metric $g$ defines a family of disintegrations via the kinetic
energy $K(q,p) = A \cdot f(g_q^{-1}(p,p)) + \tfrac12 \log|g_q| + \text{const}$,
and two new scalars — the **effective potential energy**
$\widecheck{V}(q) = V(q) + \tfrac12 \log|g_q| + \text{const}$ and the **effective
kinetic energy** $\check{K}(q,p) = A \cdot f(g_q^{-1}(p,p))$. Because $H$ is
conserved, the autocorrelations induced by these are negations of each other and
give identical times; but they **recur quickly** (at turning points for Gaussian
disintegrations), making them poor criteria.

A Riemannian metric also admits the **generalized No-U-Turn criterion**
(Betancourt 2013), which terminates when $\kappa_{\mathrm{NUTS}}(T) < 0$ using
the running momentum integral $\rho_T$. When the metric is Euclidean this reduces
to the usual No-U-Turn criterion (Hoffman & Gelman 2014); combined with
multiplicative trajectory expansion and a slice sampler this is exactly the
No-U-Turn Sampler (NUTS).

For simple level-set geometries the generalized No-U-Turn criterion is satisfied
when a trajectory has traveled from one side of a level set to the other,
matching the optimal-integration-time intuition. Weakness: because the criterion
is small near the initial point, small oscillations can cause premature vanishing;
and it is more expensive to evaluate than a numerical exhaustion.

---

## 4. Experiments

### 4.1 Graphical Experiments

Two-dimensional Gaussian target with a Euclidean-Gaussian disintegration:
effective potential $\widecheck{V}(q) = \tfrac12 q^i q^j \frac{\delta_{ij} - (1-\delta_{ij})\rho}{1-\rho^2} + \text{const}$
and effective kinetic $\check{K}(q,p) = \tfrac12 p_i p_j \delta^{ij}$ (here
$\delta_{ij}$ is the Kronecker/Dirac delta, not the threshold $\delta$).

- $\rho = 0.99$ (highly correlated): NUTS terminates prematurely; the exhaustive
  criterion yields substantially longer integration times for any $\delta$. The
  temporal expectations of $\check{K}$ and $\widecheck{V}$ vanish long before the
  exhaustion — poor criteria.
- $\rho = 0.7$ (weaker correlation): NUTS no longer suffers premature
  termination and gives more optimal times than the exhaustion.

### 4.2 Performance Experiments

Both NUTS and XHMC implemented with a second-order symplectic **leapfrog**
integrator and multiplicative trajectory expansion. NUTS uses a slice sampler
over the final trajectory; XHMC uses multinomial sampling. An integrator error
cutoff rejects any trajectory $\mathfrak{t}$ around $z_0$ satisfying
$H(z_0) - H(z) > 1000$, $\forall z \in \mathfrak{t}$. XHMC run with nominal
thresholds $\delta = 0.1$ and $\delta = 0.01$. Implemented in Stan/CmdStan
(exhaustions branch, commit c04d34ee77d831a2817cf3c7671aebc50a3bf825).

**4.2.1 IID Gaussian target** ($N=100$, $\rho=0$). Marginal energy distribution
is $\chi^2$ with 100 dof; momentum-resampling variation is $\chi^2$ with 50 dof.
Every trajectory oscillates with period $2\pi$ independent of level set, so the
optimal maximal integration time is $T(z) = 2\pi$, i.e. $L = 2\pi/\epsilon \sim
64$ leapfrog steps. NUTS integrates to about **half** of this; XHMC integrates
much longer, giving worse effective samples per transition and per leapfrog step.
The virial rate decomposes per dimension,
$\mathrm{d}G/\mathrm{d}t = 2\sum_{n=1}^N (T_n - V_n)$; these add incoherently and
the criterion is not satisfied until every dimension has oscillated a full
period, biasing samples toward the initial state.

**4.2.2 Correlated Gaussian target** with covariance
$\Sigma^{ij} = \rho^{|i-j|}$, $\rho = 0.95$. Trajectories are non-periodic but
dynamically ergodic. Up to $\sim 2^7 = 128$ leapfrog steps convergence is
superlinear; afterwards it slows to the $\sqrt{t}$ asymptotic rate. NUTS
identifies near-optimal times; nominal XHMC tunes select times in the inefficient
asymptotic regime (Table 1).

**4.2.3 Nonlinear target:** 1-PL item response theory model for 50 students,
$y_i \sim \text{Bernoulli}(\text{logistic}(\theta - b_i))$,
$b_i \sim \mathcal{N}(0, 10)$, $\theta \sim \mathcal{N}(0, 10)$ (normals given as
mean and standard deviation). The likelihood is non-identified; the posterior has
strong nonlinear correlations requiring dynamic integration times. NUTS
terminates prematurely; XHMC gives much larger effective sample sizes and higher
efficiency (Table 2). This superiority is coincidental: the nominal tunes happen
to identify long times that avoid the asymptotic regime.

---

## 5. Conclusions and Future Work

HMC efficiently explores a target when the flow is integrated long enough to
avoid diffusion but not so long as to waste computation. The analysis both
motivates NUTS and inspires the complementary XHMC algorithm. Neither criterion
robustly identifies optimal times in all cases. XHMC's stronger theoretical
foundation makes it ripe for formal analysis (step-size optimality, ergodicity).
Future directions include analyzing the marginal energy autocorrelation to
identify/optimize cotangent disintegrations, and Rao-Blackwellization by
averaging over the entire trajectory rather than sampling a single point.
