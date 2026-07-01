# RATTLE - Check Fixtures

The paper reports qualitative/analytical results rather than result tables. The
implementable fixtures are the convergence orders and the exact per-step constraint
invariants.

## Error scaling (2D pendulum in a gravitational potential)

- given: RATTLE (or SHAKE) integration of a 2D pendulum (rigid rod = one distance
  constraint) in a gravitational potential, varied time step $h$.
  expect: global error in energy, coordinates, and velocities is **quadratic in $h$**
  (i.e. $\propto h^2$). Halving $h$ reduces global error by ~4x.
- given: same, compare RATTLE vs SHAKE error coefficients.
  expect: the coefficients of the quadratic errors are **similar in magnitude** for the
  two algorithms.

## Analytical order guarantees (regression assertions)

- given: one RATTLE step from exact initial data.
  expect: local (single-step) error $= O(h^3)$; global error over fixed interval $=
  O(h^2)$ — same as velocity Verlet for unconstrained dynamics and same as SHAKE.
- given: any RATTLE step output.
  expect: energy conserved along trajectory to within $O(h^2)$.

## Per-step exact constraint invariants (unit tests for the solver)

For every constrained pair $(i,j)$ after a RATTLE step, to within the chosen tolerance:

- position invariant: $|\mathbf{r}_i(t+h) - \mathbf{r}_j(t+h)|^2 - d_{ij}^2 = 0$
  (Eq. A1). expect residual $< \text{tol}$.
- velocity invariant: $[\dot{\mathbf{r}}_i(t+h) - \dot{\mathbf{r}}_j(t+h)] \cdot
  [\mathbf{r}_i(t+h) - \mathbf{r}_j(t+h)] = 0$ (Eq. A2). expect $< \text{tol}$.
  (This velocity invariant is what RATTLE guarantees and SHAKE does NOT.)
- momentum conservation: the paired corrections $\mp g\,\mathbf{r}_{ij}/m_{i,j}$ and
  $\mp k\,\mathbf{r}_{ij}/m_{i,j}$ leave total linear momentum unchanged. expect
  $\sum_i m_i \dot{\mathbf{r}}_i$ unchanged by the constraint-correction loops.

## Solver behavior

- given: SHAKE-style position iteration (Eq. C11) and velocity iteration (Eq. C15),
  iterating over constraints.
  expect: both loops converge; each terminates when all constraints are satisfied to
  within the desired tolerance. Velocity correction $k$ (Eq. C15) is exact (no $O(g^2)$
  neglect), so the velocity loop typically converges faster than the position loop.

## Cost / storage

- given: N degrees of freedom.
  expect: storage = $3N$ floating point numbers (same as SHAKE and Verlet). RATTLE
  computes two constraint-force sets vs SHAKE's one, but overall cost is dominated by
  the intermolecular force evaluation, so wall-time is comparable to SHAKE.
