# RATTLE - Implementable Equations

Reduced/consistent units throughout; masses `m_i` explicit (no `mass=1` assumption).
`h` is the time step. Prime on sums ($\sum'$) means "over atoms $j$ constrained to $i$".

<!-- eq:2.1 -->
$$ r(t+h) = 2r(t) - r(t-h) + h^2 f[r(t)] $$
- **what:** standard (position) Verlet integrator step for $\ddot{r}=f[r]$.
- **symbols:** $r$ - Cartesian coordinates ($\mathbb{R}^{3N}$); $h$ - time step; $f[r]=\ddot{r}$ - acceleration field (force/mass); $t$ - time.

<!-- eq:2.2 -->
$$ \dot{r}(t) = [r(t+h) - r(t-h)]/2h $$
- **what:** central-difference velocity estimate from Verlet positions (error $O(h^2)$).
- **symbols:** $\dot{r}$ - velocity ($\mathbb{R}^{3N}$); other symbols as above.

<!-- eq:2.3 -->
$$ r(t+h) = r(t) + h\dot{r}(t) + h^2 f[r(t)]/2 $$
- **what:** velocity-Verlet position update (first half).
- **symbols:** $\dot{r}(t)$ - velocity at $t$; $f[r(t)]$ - acceleration at $t$.

<!-- eq:2.4 -->
$$ \dot{r}(t+h) = \dot{r}(t) + h[f[r(t)] + f[r(t+h)]]/2 $$
- **what:** velocity-Verlet velocity update; averages accelerations at $t$ and $t+h$.
- **symbols:** $f[r(t+h)]$ - acceleration evaluated at the new positions.

<!-- eq:2.7 -->
$$ r(t+h) = r(t) + h\dot{r}(t) + h^{2}[f[r(t)] + g_{RR}(t)]/2 $$
- **what:** RATTLE position step; $g_{RR}$ chosen so new positions satisfy the distance constraints (A1).
- **symbols:** $g_{RR}(t)$ - position-constraint force approximation (per unit mass) at $t$.

<!-- eq:2.8 -->
$$ \dot{r}(t+h) = \dot{r}(t) + h[f[r(t)] + g_{RR}(t) + f[r(t+h)] + g_{RV}(t)]/2 $$
- **what:** RATTLE velocity step; $g_{RV}$ chosen so new velocities satisfy the velocity constraints (A2).
- **symbols:** $g_{RV}(t)$ - velocity-constraint force approximation (per unit mass) associated with $t+h$.

<!-- eq:A0 -->
$$ \sigma_{ij}(\{\mathbf{r}(t)\}) \equiv [\mathbf{r}_i(t) - \mathbf{r}_j(t)]^2 - d_{ij}^2 $$
- **what:** distance-constraint residual for the bonded pair $(i,j)$.
- **symbols:** $\mathbf{r}_i$ - position of atom $i$ ($\mathbb{R}^3$); $d_{ij}$ - fixed target distance; $\sigma_{ij}$ - scalar constraint residual (=0 when satisfied).

<!-- eq:A2 -->
$$ \left[\dot{\mathbf{r}}_{i}(t) - \dot{\mathbf{r}}_{j}(t)\right] \cdot \left[\mathbf{r}_{i}(t) - \mathbf{r}_{j}(t)\right] = 0 $$
- **what:** velocity (time-derivative) form of the constraint: relative velocity is perpendicular to the bond.
- **symbols:** $\dot{\mathbf{r}}_i$ - velocity of atom $i$ ($\mathbb{R}^3$).

<!-- eq:A2b -->
$$ m_i \ddot{\mathbf{r}}_i = \mathbf{F}_i + \mathbf{G}_i $$
- **what:** Newton's equation with physical force plus constraint force.
- **symbols:** $m_i$ - mass of atom $i$; $\mathbf{F}_i$ - non-constraint force; $\mathbf{G}_i$ - constraint force.

<!-- eq:A2c -->
$$ \mathbf{G}_i = -\sum_j{}' \lambda_{ij}(t)\, \nabla_i \sigma_{ij} = -2\sum_j{}' \lambda_{ij}(t)\, \mathbf{r}_{ij}(t) $$
- **what:** constraint force as a sum of Lagrange multipliers times bond vectors ($\nabla_i\sigma_{ij}=2\mathbf{r}_{ij}$).
- **symbols:** $\lambda_{ij}$ - time-dependent Lagrange multiplier ($\lambda_{ij}=\lambda_{ji}$); $\mathbf{r}_{ij}\equiv\mathbf{r}_i-\mathbf{r}_j$; prime - sum over constrained neighbors of $i$.

<!-- eq:A3 -->
$$ \mathbf{r}_{i}(t+h) = \mathbf{r}_{i}(t) + h\dot{\mathbf{r}}_{i}(t) + (h^{2}/2m_{i})\left[\mathbf{F}_{i}(t) - 2\sum_{j}{}' \lambda_{RRij}(t)\, \mathbf{r}_{ij}(t)\right] $$
- **what:** per-atom RATTLE position update; solve for $\lambda_{RRij}(t)$ so (A1) holds at $t+h$.
- **symbols:** $\lambda_{RRij}(t)$ - position-stage multiplier; $\mathbf{r}_{ij}(t)=\mathbf{r}_i(t)-\mathbf{r}_j(t)$; $\mathbf{F}_i(t)$ - force at $t$.

<!-- eq:A4 -->
$$ \dot{\mathbf{r}}_{i}(t+h) = \dot{\mathbf{r}}_{i}(t) + (h/2m_{i}) \left[ \mathbf{F}_{i}(t) - 2 \sum_{j}{}' \lambda_{RRij}(t)\, \mathbf{r}_{ij}(t) + \mathbf{F}_{i}(t+h) - 2 \sum_{j}{}' \lambda_{RVij}(t+h)\, \mathbf{r}_{ij}(t+h) \right] $$
- **what:** per-atom RATTLE velocity update; solve for $\lambda_{RVij}(t+h)$ so (A2) holds at $t+h$.
- **symbols:** $\lambda_{RVij}(t+h)$ - velocity-stage multiplier; $\mathbf{F}_i(t+h)$ - force at new positions; $\mathbf{r}_{ij}(t+h)$ - new bond vector.

<!-- eq:C3 -->
$$ \mathbf{q}_i = \dot{\mathbf{r}}_i(t) + (h/2m_i)\, \mathbf{F}_i(t) - (1/m_i) \sum_j g_{ij}\, \mathbf{r}_{ij}(t) $$
- **what:** auxiliary "half-step velocity" collecting the position-stage impulses; $g_{ij}=h\lambda_{RRij}(t)$.
- **symbols:** $\mathbf{q}_i$ - auxiliary velocity ($\mathbb{R}^3$); $g_{ij}$ - scaled position multiplier.

<!-- eq:C4 -->
$$ \mathbf{r}_i(t+h) = \mathbf{r}_i(t) + h\mathbf{q}_i $$
- **what:** compact RATTLE position update in terms of $\mathbf{q}_i$.
- **symbols:** as above.

<!-- eq:C5 -->
$$ \dot{\mathbf{r}}_i(t+h) = \mathbf{q}_i + (h/2m_i)\mathbf{F}_i(t+h) - (1/m_i)\sum_j k_{ij}\mathbf{r}_{ij}(t+h) $$
- **what:** compact RATTLE velocity update; $k_{ij}=h\lambda_{RVij}(t+h)$.
- **symbols:** $k_{ij}$ - scaled velocity multiplier.

<!-- eq:C6 -->
$$ \mathbf{q}_i = \dot{\mathbf{r}}_i(t) + (h/2m_i)\, \mathbf{F}_i(t) \qquad i=1,\dots,N $$
- **what:** initial guess for $\mathbf{q}_i$ before the position-constraint iteration (all $g_{ij}=0$).
- **symbols:** as above.

<!-- eq:C7 -->
$$ \mathbf{s} = \mathbf{r}_i(t) + h\mathbf{q}_i - \mathbf{r}_j(t) - h\mathbf{q}_j $$
- **what:** current trial bond displacement for constraint $(i,j)$ during position iteration.
- **symbols:** $\mathbf{s}$ - trial displacement $\mathbf{r}_i(t+h)-\mathbf{r}_j(t+h)$ ($\mathbb{R}^3$).
<!-- CHECK: OCR wrote both terms with subscript i; corrected the second pair to atom j so s is the i-j displacement. -->

<!-- eq:C11 -->
$$ g = \frac{s^2 - d_{ij}^2}{2h\,[\mathbf{s} \cdot \mathbf{r}_{ij}(t)]\,(m_i^{-1} + m_j^{-1})} $$
- **what:** single-constraint position correction magnitude (SHAKE-style, neglecting $O(g^2)$).
- **symbols:** $s^2=|\mathbf{s}|^2$; $\mathbf{r}_{ij}(t)$ - reference bond vector at $t$; then update $\mathbf{q}_i \mathrel{-}= g\,\mathbf{r}_{ij}(t)/m_i$, $\mathbf{q}_j \mathrel{+}= g\,\mathbf{r}_{ij}(t)/m_j$.

<!-- eq:C12 -->
$$ \dot{\mathbf{r}}_i(t+h) = \mathbf{q}_i + h\mathbf{F}_i(t+h)/2m_i \qquad i=1,\dots,N $$
- **what:** initial guess for velocities before the velocity-constraint iteration (all $k_{ij}=0$).
- **symbols:** as above.

<!-- eq:C15 -->
$$ k = \frac{\mathbf{r}_{ij}(t+h) \cdot [\dot{\mathbf{r}}_i(t+h) - \dot{\mathbf{r}}_j(t+h)]}{d_{ij}^2\,(m_i^{-1} + m_j^{-1})} $$
- **what:** single-constraint velocity correction magnitude (exact, no approximation).
- **symbols:** $\mathbf{r}_{ij}(t+h)$ - new bond vector; then update $\dot{\mathbf{r}}_i \mathrel{-}= k\,\mathbf{r}_{ij}(t+h)/m_i$, $\dot{\mathbf{r}}_j \mathrel{+}= k\,\mathbf{r}_{ij}(t+h)/m_j$.

## Derivations note

The order-of-error analysis in Appendix B (Eqs. B1-B4) is proof algebra establishing
that RATTLE's local error is $O(h^3)$ and global error is $O(h^2)$; it is not
implemented. See `paper.md` Appendix B.
