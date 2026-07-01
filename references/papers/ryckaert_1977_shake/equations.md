# Equations - Ryckaert, Ciccotti, Berendsen (1977): SHAKE constraint dynamics

Reduced/CGS units as in the original (positions in Å, time in s). No unit
convention beyond consistency of $m$, $h$, and forces.

<!-- eq:1 -->
$$\sigma_k(\{\mathbf{r}\}) = 0 \qquad (k = 1, ..., l).$$
- **what:** general holonomic constraint k; l constraints total.
- **symbols:** $\sigma_k$ - kth constraint function; $\{\mathbf{r}\}$ - all positions ($\mathbb{R}^{3N}$); $l$ - number of constraints; $k$ - constraint index.

<!-- eq:2 -->
$$\mathbf{G}_i = -\sum_{k=1}^{l} \lambda_k(t) \, \nabla_i \sigma_k.$$
- **what:** constraint force on particle i as a sum of Lagrange multipliers times constraint gradients.
- **symbols:** $\mathbf{G}_i$ - constraint force on particle i ($\mathbb{R}^3$); $\lambda_k(t)$ - kth Lagrange multiplier (time-dependent scalar); $\nabla_i\sigma_k$ - gradient of $\sigma_k$ w.r.t. $\mathbf{r}_i$.

<!-- eq:2.1 -->
$$\sigma_k(\{\mathbf{r}(t)\}) \equiv (\mathbf{r}_i(t) - \mathbf{r}_j(t))^2 - d_{ij}^2 = 0 \qquad (k = 1,...,l).$$
- **what:** rigid-bond distance constraint: squared separation of the pair (i,j) equals fixed $d_{ij}^2$.
- **symbols:** $\mathbf{r}_i,\mathbf{r}_j$ - positions of the bonded pair; $d_{ij}$ - fixed bond length; $(\cdot)^2$ - squared Euclidean norm.

<!-- eq:2.2 -->
$$m_i \ddot{\mathbf{r}}_i = \mathbf{F}_i + \mathbf{G}_i = -\nabla_i V - \sum_{k=1}^{l} \lambda_k \nabla_i \sigma_k \qquad (i = 1, ..., N).$$
- **what:** Lagrange equation of motion of the first kind; total force = physical force + constraint force.
- **symbols:** $m_i$ - mass of particle i; $\ddot{\mathbf{r}}_i$ - acceleration; $\mathbf{F}_i=-\nabla_i V$ - physical force; $V$ - potential energy; $N$ - number of particles.

<!-- eq:3.1 -->
$$u(h) = -u(-h) + 2u(0) + h^2\ddot{u}(0) + O(h^4).$$
- **what:** Verlet position update (any coordinate u), used as the base integrator.
- **symbols:** $u(h)$ - coordinate at $t_0+h$; $u(-h)$ - previous step; $u(0)$ - current; $\ddot{u}(0)$ - acceleration = force/mass; $h$ - time step.

<!-- eq:3.2 -->
$$\dot{u}(0) = (u(h) - u(-h))/2h + O(h^2).$$
- **what:** central-difference velocity estimate at the current step.
- **symbols:** $\dot{u}(0)$ - velocity at $t_0$; $h$ - time step.

<!-- eq:3.3 -->
$$\mathbf{r}_{i}'(h) = -\mathbf{r}_{i}(-h) + 2\mathbf{r}_{i}(0) + (h^{2}/m_{i}) \mathbf{F}_{i}(0), \qquad \delta \mathbf{r}_{i}(h) = (h^{2}/m_{i}) \sum_{k=1}^{l} \gamma_{k} [\nabla_{i} \sigma_{k}]_{t_{0}}.$$
- **what:** unconstrained Verlet position $\mathbf{r}_i'$, plus the constraint correction $\delta\mathbf{r}_i$; constrained position is $\mathbf{r}_i(h)=\mathbf{r}_i'(h)+\delta\mathbf{r}_i(h)$.
- **symbols:** $\mathbf{r}_i'(h)$ - unconstrained new position; $\delta\mathbf{r}_i$ - constraint displacement; $\gamma_k$ - undetermined parameter replacing $\lambda_k^{(n-2)}$; $[\nabla_i\sigma_k]_{t_0}$ - constraint gradient at $t_0$ ($=2\mathbf{r}_{ij}$ for bond constraints).

<!-- eq:3.4 -->
$$\sigma_k(\{\mathbf{r}_i(h, \{\gamma_k\})\}) = 0 \qquad (k = 1, ..., l).$$
- **what:** the l parameters $\gamma_k$ are fixed by requiring all constraints hold exactly at the new step.
- **symbols:** as above; $\mathbf{r}_i(h,\{\gamma_k\})=\mathbf{r}_i'(h)+\delta\mathbf{r}_i(h)$.

<!-- eq:3.5 -->
$$\gamma_k(0) = \lambda_k(0) + O(h^2) \qquad (k = 1, ..., l).$$
- **what:** the parameter equals the true Lagrange multiplier to $O(h^2)$; error consistent with Verlet.
- **symbols:** $\gamma_k(0)$ - parameter; $\lambda_k(0)$ - true multiplier at $t_0$.

<!-- eq:3.6 -->
$$(\mathbf{r}_i(t) - \mathbf{r}_j(t))^2 - d_{ij}^2 = 0.$$
- **what:** rigid bond constraint (same as 2.1), specialized form used in Section 3.
- **symbols:** as 2.1.

<!-- eq:3.7 -->
$$2(\mathbf{r}_{j}'(h) - \mathbf{r}_{i}'(h)) \cdot \left(-h^{2} \sum_{k=1}^{l} \gamma_{k} \left[ \left(\tfrac{\nabla_{j}}{m_{j}} - \tfrac{\nabla_{i}}{m_{i}}\right) \sigma_{k} \right]_{0}\right) + h^{4} \sum_{k=1}^{l} \sum_{k'=1}^{l} \gamma_{k} \gamma_{k'} \left[ \left(\tfrac{\nabla_{j}}{m_{j}} - \tfrac{\nabla_{i}}{m_{i}}\right) \sigma_{k} \right]_{0} \cdot \left[ \left(\tfrac{\nabla_{j}}{m_{j}} - \tfrac{\nabla_{i}}{m_{i}}\right) \sigma_{k'} \right]_{0} = d_{ij}^{2} - (\mathbf{r}_{j}'(h) - \mathbf{r}_{i}'(h))^{2}.$$
- **what:** quadratic system for the $\gamma_k$ (matrix method); linear term $\propto h^2$, quadratic term $\propto h^4$. Solve by iteration starting from $\gamma_k=0$.
- **symbols:** $\mathbf{r}_i'(h)$ - unconstrained positions; $\gamma_k,\gamma_{k'}$ - parameters; $m_i,m_j$ - masses of the bonded pair; the bracketed gradient difference $(\nabla_j/m_j - \nabla_i/m_i)\sigma_k$ is the mass-weighted constraint gradient.

<!-- eq:4.1 -->
$$\langle |\delta \mathbf{r}(t)| \rangle = \left[ \sum_{i=1}^{256} (\mathbf{r}_{i}^{(1)}(t) - \mathbf{r}_{i}^{(2)}(t))^{2}/256 \right]^{1/2}.$$
- **what:** RMS deviation between the generalized-coordinate (method 1) and Cartesian (method 2) trajectories over 256 particles.
- **symbols:** $\mathbf{r}_i^{(1)}$ - method-1 position; $\mathbf{r}_i^{(2)}$ - method-2 position; 256 = number of particles (64 butanes x 4).

<!-- eq:4.2 -->
$$\langle r(t) \rangle = \left[ \sum_{i=1}^{256} (\mathbf{r}_i^{(1)}(t) - \mathbf{r}_i^{(1)}(0))^2 / 256 \right]^{1/2}.$$
- **what:** RMS displacement of method-1 trajectory from its own initial configuration (reference scale).
- **symbols:** as 4.1.

<!-- eq:5.2 -->
$$\delta \mathbf{r}_{i} = -\frac{1}{m_{i}} (\Delta t)^{2} \sum_{k=1}^{l} \gamma_{k} [\nabla_{i} \sigma_{k}]_{t_{0}} = -\frac{2(\Delta t)^{2}}{m_{i}} \sum_{k=1}^{l} \gamma_{k} \mathbf{r}_{ij}(t_{0}).$$
- **what:** Verlet constraint displacement; uses $\nabla_i\sigma_k = 2\mathbf{r}_{ij}$ for bond constraints.
- **symbols:** $\Delta t$ - time step (=h); $\mathbf{r}_{ij}=\mathbf{r}_i-\mathbf{r}_j$ - kth bond vector; $\gamma_k$ - parameter.

<!-- eq:5.3 -->
$$\delta \mathbf{r}_i = \sum_{j} g_{ij} \mathbf{r}_{ij}(t_0) / m_i.$$
- **what:** total constraint correction on particle i as a sum over its constrained partners; $g_{ij}=g_{ji}$, $g_{ij}=0$ for non-constrained pairs.
- **symbols:** $g_{ij} = -2(\Delta t)^2 \gamma_k$ - scalar constraint magnitude for pair (i,j); $\mathbf{r}_{ij}(t_0)$ - bond vector at $t_0$ (fixed reference direction).

<!-- eq:5.4a -->
$$\delta^k \mathbf{r}_i = g_{ij} \mathbf{r}_{ij}(t_0) / m_i.$$
- **what:** SHAKE per-constraint correction to particle i (partner of the sweep).
- **symbols:** $\delta^k\mathbf{r}_i$ - correction from constraint k only; $g_{ij}$ - solved per constraint.

<!-- eq:5.4b -->
$$\delta^k \mathbf{r}_j = -g_{ij} \mathbf{r}_{ij}(t_0) / m_j.$$
- **what:** equal-and-opposite SHAKE correction to partner j (conserves momentum along the bond).
- **symbols:** as 5.4a.

<!-- eq:5.5 -->
$$(\mathbf{r}' + \delta \mathbf{r})^2 - d^2 = 0, \qquad \delta \mathbf{r} = \left(\frac{1}{m_i} + \frac{1}{m_j}\right) g\,\mathbf{r}.$$
- **what:** single-constraint condition after applying the correction; $\mathbf{r}'$ is the current (partially corrected) inter-particle vector.
- **symbols:** $\mathbf{r}=\mathbf{r}_i(t_0)-\mathbf{r}_j(t_0)$ - reference bond vector at $t_0$; $\mathbf{r}'$ - current pair vector after unconstrained step and earlier corrections; $g=g_{ij}$; $d=d_{ij}$ - target length.

<!-- eq:5.6 -->
$$2\left(\frac{1}{m_i} + \frac{1}{m_j}\right) g\, (\mathbf{r} \cdot \mathbf{r}') + \left(\frac{1}{m_i} + \frac{1}{m_j}\right)^2 g^2 \mathbf{r}^2 = d^2 - \mathbf{r}'^2.$$
- **what:** the per-constraint quadratic in g solved by SHAKE. For efficiency drop the $g^2$ term (first order): $g \approx (d^2-\mathbf{r}'^2)/[2(1/m_i+1/m_j)(\mathbf{r}\cdot\mathbf{r}')]$. Iterate over all constraints until residuals < tolerance.
- **symbols:** $m_i,m_j$ - masses of the pair; $\mathbf{r}\cdot\mathbf{r}'$ - dot product of reference and current bond vectors; $d^2-\mathbf{r}'^2$ - current constraint violation.

## Appendix - n-alkane semirigid chain

<!-- eq:A.1 -->
$$(\mathbf{r}_{i+1} - \mathbf{r}_i)^2 - a^2 = 0 \qquad (i = 1,..., n-1).$$
- **what:** adjacent-group bond constraint (C-C bond).
- **symbols:** $a$ - C-C bond length (= 1.53 Å); $n$ - chain length (groups).

<!-- eq:A.2 -->
$$(\mathbf{r}_{i+2} - \mathbf{r}_i)^2 - b^2 = 0 \qquad (i = 1,..., n-2), \quad b = 2a \sin(\theta/2).$$
- **what:** second-neighbor constraint encoding the fixed C-C-C bond angle via a fixed 1-3 distance b.
- **symbols:** $b$ - 1-3 distance; $\theta$ - C-C-C angle (= 109.28deg); $a$ - bond length.

<!-- eq:A.4 -->
$$\mathbf{G}_{i} = -2\lambda_{2i-4}(\mathbf{r}_{i} - \mathbf{r}_{i-2}) - 2\lambda_{2i-3}(\mathbf{r}_{i} - \mathbf{r}_{i-1}) + 2\lambda_{2i-1}(\mathbf{r}_{i+1} - \mathbf{r}_{i}) + 2\lambda_{2i}(\mathbf{r}_{i+2} - \mathbf{r}_{i}).$$
- **what:** total constraint force on chain group i from its four incident constraints (two bond, two angle); drop terms for nonexistent neighbors at chain ends.
- **symbols:** odd-indexed $\lambda$ = bond-constraint multipliers, even-indexed = angle-constraint multipliers; group indexing $i$ from 1 to n.

<!-- eq:A.6 -->
$$\mathbf{G}_{i} = 2\big[-\gamma_{2i-4}(\mathbf{r}_{i} - \mathbf{r}_{i-2}) - \gamma_{2i-3}(\mathbf{r}_{i} - \mathbf{r}_{i-1}) + \gamma_{2i-1}(\mathbf{r}_{i+1} - \mathbf{r}_{i}) + \gamma_{2i}(\mathbf{r}_{i+2} - \mathbf{r}_{i})\big].$$
- **what:** same as A.4 but with the parameters $\gamma_k$ replacing the multipliers $\lambda_k$; used inside the Verlet correction $\delta\mathbf{r}_i(h)=(h^2/m_i)\mathbf{G}_i$.
- **symbols:** $\gamma_k$ - $(2n-3)$ undetermined parameters; indices as A.4.

## Derivations note

The high-order Taylor construction (eqs 2.3, 2.4, 2.4b, 2.4c, 2.5b) and the
error-order analysis (2.6-2.11, 2.8) are derivation/justification steps kept in
`paper.md` under the `Derivation` marker; they establish that replacing
$\lambda_k^{(n-2)}$ by $\gamma_k$ preserves $O(\Delta t^{m+1})$ accuracy. The
implementable core is: unconstrained Verlet step (3.3) + solve constraints for the
correction, either by the matrix quadratic (3.7 / A.7) or by the SHAKE sweep
(5.4-5.6). The Appendix matrix elements A.8-odd / A.8-even give the explicit
banded A matrix for the n-alkane chain.
