# Equations - Duane et al. 1987, Hybrid Monte Carlo

Notation: $\phi$ = field/configuration (position analog), $\pi$ = conjugate momentum, $S(\phi)$ = action (potential energy analog), $H$ = Hamiltonian, $\tau$ = fictitious "computer" time, $\delta\tau$ = integration step. Detailed-balance / Metropolis structure maps directly onto molecular HMC with $S \to \beta U$ (potential), $\phi \to$ coordinates, $\pi \to$ momenta.

<!-- eq:1 -->
$$\langle \Omega \rangle = \frac{1}{Z} \int [d\phi]\, \exp[-S(\phi)]\, \Omega(\phi)$$
- **what:** Boltzmann-weighted expectation value of an observable.
- **symbols:** $\langle\Omega\rangle$ - expectation of observable; $Z$ - partition function; $\phi$ - field configuration; $S(\phi)$ - action (dimensionless energy); $\Omega(\phi)$ - observable.

<!-- eq:2 -->
$$Z = \int [d\phi]\, \exp[-S(\phi)]$$
- **what:** partition function (normalizing constant).
- **symbols:** as eq:1.

<!-- eq:3 -->
$$\bar{\Omega} = \frac{1}{T} \sum_{t=1}^{T} \Omega(\phi_t)$$
- **what:** Monte Carlo estimator of $\langle\Omega\rangle$ from a chain of $T$ samples.
- **symbols:** $\bar\Omega$ - sample-mean estimator; $\phi_t$ - configuration at chain step $t$; $T$ - number of samples.

<!-- eq:4 -->
$$\bar{\Omega} = \langle \Omega \rangle + \mathcal{O}(1/\sqrt{T})$$
- **what:** estimator error scales as $1/\sqrt{T}$.
- **symbols:** as above.

<!-- eq:5 -->
$$P_S(\phi)\, P_M(\phi \mapsto \phi') = P_S(\phi')\, P_M(\phi' \mapsto \phi)$$
- **what:** detailed-balance condition guaranteeing $P_S$ is the stationary distribution.
- **symbols:** $P_S(\phi) = (1/Z)\exp[-S(\phi)]$ - target distribution; $P_M(\phi\mapsto\phi')$ - Markov transition probability.

<!-- eq:6 -->
$$P_{A}(\phi \mapsto \phi') = \min\!\left(1,\ \frac{P_{S}(\phi')\, P_{C}(\phi' \mapsto \phi)}{P_{S}(\phi)\, P_{C}(\phi \mapsto \phi')}\right)$$
- **what:** generalized Metropolis acceptance probability for arbitrary proposal $P_C$.
- **symbols:** $P_A$ - accept probability; $P_C(\phi\mapsto\phi')$ - candidate/proposal probability; $P_S$ - target.

<!-- eq:7 -->
$$H'(\phi,\pi) \equiv \tfrac{1}{2}\pi^2 + S'(\phi)$$
- **what:** guidance Hamiltonian used to propose moves (kinetic + arbitrary guidance action).
- **symbols:** $H'$ - guidance Hamiltonian; $\pi$ - conjugate momentum ($\tfrac12\pi^2$ = kinetic energy, unit mass); $S'(\phi)$ - guidance action (may differ from $S$).

<!-- eq:8 -->
$$\dot{\phi} = \delta H'/\delta \pi = \pi, \qquad \dot{\pi} = -\delta H'/\delta \phi = -\delta S'/\delta \phi$$
- **what:** Hamilton's equations of motion driving the deterministic proposal trajectory.
- **symbols:** $\dot\phi,\dot\pi$ - $\tau$-derivatives; $\delta S'/\delta\phi$ - functional gradient (force).

<!-- eq:9 -->
$$H(\phi, \pi) \equiv \tfrac{1}{2}\pi^2 + S(\phi)$$
- **what:** acceptance (true) Hamiltonian; defines the sampled equilibrium distribution.
- **symbols:** $H$ - acceptance Hamiltonian; $S(\phi)$ - target action of eq:1.

<!-- eq:10 -->
$$P_{G}(\pi) \propto \exp(-\pi^2/2)$$
- **what:** momentum-refresh (heatbath): draw fresh momenta from a unit-variance Gaussian each iteration.
- **symbols:** $P_G(\pi)$ - Gaussian momentum distribution, mean 0, unit variance (mass = 1).

<!-- eq:11 -->
$$P_{H}\big((\phi, \pi) \mapsto (\phi', \pi')\big) = \delta\big[(\phi', \pi') - (\phi(\tau_{0}), \pi(\tau_{0}))\big]$$
- **what:** deterministic trajectory proposal is a delta function at the time-$\tau_0$ endpoint.
- **symbols:** $P_H$ - proposal kernel in phase space; $\tau_0$ - total trajectory time; $(\phi(\tau_0),\pi(\tau_0))$ - endpoint of integrating eq:8.

<!-- eq:12 -->
$$P_{A}\big((\phi,\pi)\mapsto(\phi',\pi')\big) = \min\!\big(1,\ \exp(\delta H)\big), \qquad \delta H = H(\phi,\pi) - H(\phi',\pi')$$
- **what:** HMC Metropolis accept step on the energy change of the trajectory. Equivalent to $\min(1,\exp(-\Delta H))$ with $\Delta H = H_{\text{new}} - H_{\text{old}}$.
- **symbols:** $\delta H$ - negative energy change; $H$ - acceptance Hamiltonian eq:9. <!-- CHECK: paper OCR "5H=-H(O',re')-H(O,re)" garbled; sign chosen so identity eq:15 (min(1,exp(-δH)) with δH=H_new-H_old) is consistent with accept=min(1,exp(δH_here)). -->

<!-- eq:13 -->
$$P_{M}(\phi \mapsto \phi') = \int [d\pi]\, [d\pi']\, P_{G}(\pi)\, P_{H}\big((\phi, \pi) \mapsto (\phi', \pi')\big)\, P_{A}\big((\phi, \pi) \mapsto (\phi', \pi')\big)$$
- **what:** full HMC transition on the field: momentum draw $\times$ deterministic trajectory $\times$ accept, marginalized over momenta.
- **symbols:** $P_M$ - marginal field transition; $P_G,P_H,P_A$ - eqs 10, 11, 12.

<!-- eq:18 -->
$$\pi(\delta \tau/2) = \pi(0) - [\delta S(0)/\delta \phi]\, \delta \tau/2$$
- **what:** leapfrog initial half-kick of momentum.
- **symbols:** $\delta\tau$ - integrator step; $\delta S/\delta\phi$ - force (gradient of action).

<!-- eq:19 -->
$$\phi(\tau + \delta \tau) = \phi(\tau) + \pi(\tau + \delta \tau/2)\cdot \delta \tau$$
- **what:** leapfrog drift: full position update using half-step momentum.
- **symbols:** as above.

<!-- eq:20 -->
$$\pi(\tau + \delta \tau/2) = \pi(\tau - \delta \tau/2) - [\delta S(\tau)/\delta \phi]\cdot \delta \tau$$
- **what:** leapfrog full-kick of momentum at intermediate steps ($n-1$ of them).
- **symbols:** as above; $\delta S(\tau)/\delta\phi$ - force evaluated at current $\phi(\tau)$.

<!-- eq:21 -->
$$\pi(\tau_0) = \pi(\tau_0 - \delta \tau/2) - [\delta S(\tau_0)/\delta \phi]\cdot \delta \tau/2$$
- **what:** leapfrog final half-kick of momentum; completes reversible, area-preserving map.
- **symbols:** $\tau_0 = n\,\delta\tau$ - trajectory endpoint.

<!-- eq:22 -->
$$P_S(\phi) = \frac{1}{Z} \int [d\bar{\psi}][d\psi]\, \exp[-S(\phi) - \bar{\psi}\mathcal{M}\psi] = \frac{1}{Z'} \det(\mathcal{M})\, \exp[-S(\phi)] = \frac{1}{Z''} \int [d\chi^*][d\chi]\, \exp[-S(\phi) - \chi^{*}(\mathcal{M}^{\dagger}\mathcal{M})^{-1}\chi]$$
- **what:** pseudofermion representation - integrate out Grassmann fermions to a bosonic $\chi$ field with the fermion-matrix inverse in the action.
- **symbols:** $\bar\psi,\psi$ - Grassmann fermion fields; $\mathcal{M}$ - fermion matrix; $\chi,\chi^*$ - bosonic pseudofermion fields; $\det\mathcal{M}$ - fermion determinant.

<!-- eq:23 -->
$$H(\phi, \pi) = \tfrac{1}{2}\pi^2 + S(\phi) + \chi^{*}(\mathcal{M}^{\dagger}\mathcal{M})^{-1}\chi$$
- **what:** Hamiltonian including pseudofermion contribution.
- **symbols:** as eqs 9, 22.

<!-- eq:24 -->
$$\dot{\phi} = \pi, \qquad \dot{\pi} = -\delta S/\delta \phi + \chi^{*}(\mathcal{M}^{\dagger}\mathcal{M})^{-1}\big[\mathcal{M}^{\dagger}\, \delta\mathcal{M}/\delta \phi + (\delta \mathcal{M}^{\dagger}/\delta \phi)\, \mathcal{M}\big](\mathcal{M}^{\dagger}\mathcal{M})^{-1}\chi$$
- **what:** equations of motion including the pseudofermion force; requires one conjugate-gradient solve of $(\mathcal{M}^\dagger\mathcal{M})^{-1}\chi$ per MD step.
- **symbols:** $\delta\mathcal{M}/\delta\phi$ - derivative of fermion matrix wrt field; $\chi$ held fixed during MD, refreshed by heatbath $\chi=\mathcal{M}^\dagger\eta$, $\eta$ Gaussian.

<!-- eq:25 -->
$$H = \frac{1}{2} \sum_{x\mu} \pi_{\mu}^{2}(x) + S(U, \bar{\psi}, \psi)$$
- **what:** lattice-QED test Hamiltonian.
- **symbols:** $x$ - lattice site; $\mu$ - direction index; $\pi_\mu(x)$ - momentum conjugate to link; $U$ - gauge link variables.

<!-- eq:26 -->
$$S(U, \bar{\psi}, \psi) = \beta \sum_{\substack{x\nu\mu \\ \nu>\mu}} \left[ 1 - \operatorname{Re} U_{\mu\nu}(x) \right] + \sum_{x\nu} \bar{\psi}(x)\, (-D^2 + m^2)^{-1}\, \psi(y)$$
- **what:** Wilson gauge action plus staggered-fermion term for compact QED.
- **symbols:** $\beta$ - coupling constant; $U_{\mu\nu}(x)$ - plaquette (product of links); $\operatorname{Re}$ - real part; $-D^2$ - lattice Laplacian/Dirac operator; $m$ - fermion mass.
