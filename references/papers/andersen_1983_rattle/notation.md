# RATTLE - Notation and Conventions

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $t$ | time | time | current step |
| $h$ | integration time step | time | fixed; can change between steps |
| $N$ | number of degrees of freedom (atoms) | int | storage is $3N$ reals |
| $r$, $\mathbf{r}_i$ | Cartesian positions (all / atom $i$) | length, $\mathbb{R}^{3N}$ / $\mathbb{R}^3$ | Cartesian per-atom |
| $\dot{\mathbf{r}}_i$ | velocity of atom $i$ | length/time, $\mathbb{R}^3$ | |
| $\ddot{\mathbf{r}}_i$ | acceleration of atom $i$ | length/time$^2$ | |
| $\mathbf{r}_{ij}$ | bond vector $\mathbf{r}_i - \mathbf{r}_j$ | $\mathbb{R}^3$ | $\mathbf{r}_{ij}=-\mathbf{r}_{ji}$ |
| $m_i$ | mass of atom $i$ | mass | appears explicitly (NOT reduced to 1) |
| $f[r]$ | acceleration field $= \ddot r$ | length/time$^2$ | force divided by mass |
| $\mathbf{F}_i$ | physical (non-constraint) force on atom $i$ | force, $\mathbb{R}^3$ | inter + intramolecular non-constraint |
| $\mathbf{G}_i$ | constraint force on atom $i$ | force, $\mathbb{R}^3$ | $=-2\sum'_j \lambda_{ij}\mathbf{r}_{ij}$ |
| $g[r,\dot r]$ | constraint acceleration (exact) | length/time$^2$ | contains Lagrange multipliers |
| $g_{RR}(t)$ | RATTLE position-stage constraint accel. | length/time$^2$ | makes positions satisfy (A1) |
| $g_{RV}(t)$ | RATTLE velocity-stage constraint accel. | length/time$^2$ | makes velocities satisfy (A2) |
| $g_s$ | SHAKE constraint approximation | length/time$^2$ | |
| $\sigma_{ij}$ | distance-constraint residual | length$^2$ | $=[\mathbf{r}_i-\mathbf{r}_j]^2-d_{ij}^2$; $\sigma_{ij}=\sigma_{ji}$; $=0$ satisfied |
| $d_{ij}$ | fixed target bond distance for pair $(i,j)$ | length | constant |
| $\lambda_{ij}$ | exact Lagrange multiplier | | $\lambda_{ij}=\lambda_{ji}$, time-dependent |
| $\lambda_{RRij}(t)$ | RATTLE position-stage multiplier | | solved so (A1) holds at $t+h$ |
| $\lambda_{RVij}(t+h)$ | RATTLE velocity-stage multiplier | | solved so (A2) holds at $t+h$ |
| $g_{ij}$ | scaled position multiplier $=h\lambda_{RRij}(t)$ | | Appendix C |
| $k_{ij}$ | scaled velocity multiplier $=h\lambda_{RVij}(t+h)$ | | Appendix C |
| $\mathbf{q}_i$ | auxiliary half-step velocity | $\mathbb{R}^3$ | Appendix C working variable |
| $\mathbf{s}$ | trial bond displacement during position iter. | $\mathbb{R}^3$ | $\mathbf{r}_i(t+h)-\mathbf{r}_j(t+h)$ |
| $g$ (scalar) | per-constraint position correction magnitude | | Eq. C11 |
| $k$ (scalar) | per-constraint velocity correction magnitude | | Eq. C15 |

## Conventions

- Primed sum $\sum_j{}'$ = sum only over atoms $j$ bonded to $i$ by a constraint.
- $\nabla_i \sigma_{ij} = 2(\mathbf{r}_i - \mathbf{r}_j) = 2\mathbf{r}_{ij}$.
- Constraints treated: fixed pairwise distances only (bond-length constraints). Bond
  angles are imposed as extra fixed distances between the outer atoms.
- Two-stage structure: (1) position constraints via SHAKE-like iteration on $\sigma_{ij}=0$;
  (2) velocity constraints via a separate iteration enforcing
  $\mathbf{r}_{ij}\cdot(\dot{\mathbf{r}}_i-\dot{\mathbf{r}}_j)=0$.
- Correction sign convention: atom $i$ gets $-g\,\mathbf{r}_{ij}/m_i$ (or $-k\,\mathbf{r}_{ij}/m_i$),
  atom $j$ gets $+g\,\mathbf{r}_{ij}/m_j$ (or $+k\,\mathbf{r}_{ij}/m_j$) -> conserves momentum.
- Error orders: local error $O(h^3)$, global error $O(h^2)$; energy conserved to $O(h^2)$.
- Memory: $3N$ locations; positions at $t+h$ overwrite positions at $t$, and $\mathbf{q}_i$
  overwrite velocities at $t$, in place.
