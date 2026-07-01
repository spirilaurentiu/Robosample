# Checks / fixtures - Ryckaert, Ciccotti, Berendsen (1977)

## Model parameters (n-alkane / n-butane)

- C-C bond length: $a = 1.53\,\text{Å}$.
- C-C-C bond angle: $\theta = 109.28^\circ$ (109deg28'), i.e. tetrahedral-ish.
- Fixed 1-3 distance: $b = 2a\sin(\theta/2)$. Given: with $a=1.53$ Å,
  $\theta=109.47^\circ$ tetrahedral, $b = 2(1.53)\sin(54.735^\circ) \approx 2.499\,\text{Å}$.
  With $\theta=109.28^\circ$: $b \approx 2(1.53)\sin(54.64^\circ) \approx 2.497\,\text{Å}$.
  (expect: 1-3 distance $\approx 2.50$ Å for a C-C-C group)
- n-butane: 4 point groups, 3 rigid bonds + 2 rigid angles = **5 constraints/molecule**.
- n-alkane (n groups): $(n-1)$ bond + $(n-2)$ angle = $(2n-3)$ constraints;
  $(n-3)$ internal rotational DOF.

## n-Butane liquid MD benchmark

- System: 64 n-butane molecules = 256 particles, cubic box, periodic BC.
- Density: $\rho = 0.675\,\text{g/cm}^3$; temperature $T \approx 200\,\text{K}$.
- Time step: $h = 1.95\times 10^{-15}\,\text{s}$.
- Comparison period: $T = 1.56\times 10^{-13}\,\text{s}$ (~80 steps).
- At $t=t_0+T$, deviation between generalized-coord (method 1) and Cartesian-SHAKE
  (method 2) trajectories:
  - given: `identical initial config, h=1.95e-15 s, run T=1.56e-13 s`
  - expect: $\langle|\delta\mathbf{r}(t_0+T)|\rangle = 1.1\times10^{-4}\,\text{Å}$
  - expect: $\langle r(t_0+T)\rangle = 0.6\,\text{Å}$ (self-displacement scale)
  - expect: relative average discrepancy $\approx 2\times10^{-4}$.
- Total energy: no drift over long times; oscillates with amplitude
  $\sim 10^{-3}$ of the kinetic energy. (regression: energy conservation to ~0.1% KE)
- CPU time per integration step (IBM 370/168):
  - method 1 (generalized coords + Gear): 3.25 s/step
  - method 2 (Cartesian + undetermined params): 1.30 s/step
  - Cartesian is ~2.5x faster per step.

## Constraint-satisfaction convergence (matrix method, eq 3.7)

- On liquid n-butane: **3-4 iterations** of the quadratic solve suffice to satisfy
  constraints to relative displacement discrepancy of order $10^{-10}$.
  - given: `Verlet MD step on butane, iterate eq 3.7 from gamma_k=0`
  - expect: relative constraint residual $\sim 10^{-10}$ after 3-4 iterations.

## SHAKE vs matrix method benchmark (single decane)

- System: single decane molecule (10 groups, 17 constraints), time step
  $4\times10^{-16}\,\text{s}$.
- Both give identical numerical results within the specified tolerance.
- CPU (CDC Cyber 74-16):
  - At relative tolerance $10^{-7}$: SHAKE and matrix method use the **same**
    CPU time ($\approx 90\,\text{ms}$ for the resetting step).
  - Higher accuracy (tighter tolerance): matrix method faster.
  - Lower accuracy: SHAKE preferred.
  - SHAKE iteration count and CPU grow roughly $\propto -\log(\text{tolerance})$.
  - At tolerance $10^{-10}$: SHAKE ~2x slower than at $10^{-7}$.

## Algorithmic invariants (unit-test targets)

- After a constrained step, every bond/angle constraint residual
  $|\sigma_k| < \text{tol}$ (SHAKE guarantees this to tolerance; matrix guarantees
  it to iteration convergence).
- SHAKE correction is momentum-conserving along each bond: the corrections on the
  two partners are equal and opposite in mass-weighted form
  ($m_i\delta^k\mathbf{r}_i = -m_j\delta^k\mathbf{r}_j$ per constraint), so total
  linear momentum is preserved.
- Constraint force acts purely along the current bond direction
  $\mathbf{r}_{ij}(t_0)$ (central force), so it does no work on the rigid DOF.
