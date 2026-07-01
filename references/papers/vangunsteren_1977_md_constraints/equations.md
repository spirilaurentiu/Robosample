# Equations - van Gunsteren & Berendsen 1977 (MD & constraint dynamics)

<!-- eq:2.1 -->
$$\frac{d^2 \mathbf{r}_i}{dt^2} = \mathbf{F}_i(\mathbf{r}_1, \mathbf{r}_2, \dots, \mathbf{r}_N)/m_i, \quad i = 1, 2, \dots N,$$
- **what:** Newton's equations of motion for N particles (the MD system to integrate).
- **symbols:** $\mathbf{r}_i$ - position of particle i (R^3); $\mathbf{F}_i$ - force on i (R^3); $m_i$ - mass of i; $N$ - number of particles; $t$ - time.

<!-- eq:2.2 -->
$$\mathbf{F}_i(\mathbf{r}_1, \mathbf{r}_2, \dots \mathbf{r}_N) = -\boldsymbol{\nabla}_i V(\mathbf{r}_1, \mathbf{r}_2, \dots \mathbf{r}_N), \quad i = 1, 2, \dots N.$$
- **what:** Conservative force is minus the gradient of the potential energy.
- **symbols:** $V$ - potential energy (scalar); $\boldsymbol{\nabla}_i$ - gradient w.r.t. $\mathbf{r}_i$.

<!-- eq:2.3a -->
$$y'' = f(y),$$
- **what:** Second-order initial value problem form of the EOM; $f$ is the (mass-scaled) force function.
- **symbols:** $y, y', y''$ - 3N-dimensional state vector and its time derivatives; $f$ - acceleration function (double sum over particles).

<!-- eq:2.4 -->
$$t_n = nh, \quad n = 0, 1, 2, \dots,$$
- **what:** Equally-spaced mesh (time grid).
- **symbols:** $h$ - time step; $n$ - step index; $t_n$ - time at step n.

<!-- eq:2.5 -->
$$y_n, \; hy'_n, \; h^2 y''_n/2, \; \dots, \; h^{k-1} y_n^{(k-1)}/(k-1)!$$
- **what:** Stored state in the Nordsieck (N) representation: scaled derivatives up to order k-1.
- **symbols:** $y_n^{(j)}$ - jth time derivative at step n; $k$ - number of values (order = k-1); $h$ - time step.

<!-- eq:2.6 -->
$$y_n, \; hy'_n, \; hy'_{n-1}, \; \dots, \; hy'_{n-k+2}$$
- **what:** Stored state in the Adams representation (value plus past scaled first derivatives).
- **symbols:** as above; $y'_{n-j}$ - first derivative at previous mesh points.

<!-- eq:2.7 -->
$$\mathbf{y}_n(N) \equiv [y_n, \; hy'_n, \; h^2 y''_n/2, \; \dots, \; h^{k-1} y_n^{(k-1)}/(k-1)!]^{T}.$$
- **what:** N-representation column vector (state carried between steps).
- **symbols:** $\mathbf{y}_n(N)$ - Nordsieck state column vector at step n; T - transpose.

<!-- eq:2.8 -->
$$\mathbf{y}_{n+1,\,(p)} = \mathbf{A}\mathbf{y}_n,$$
- **what:** Predictor step: extrapolate the state with matrix A.
- **symbols:** $\mathbf{A}$ - predictor matrix (Pascal triangle in N-rep, eq A.1); $(p)$ - predicted value.

<!-- eq:2.9 -->
$$\mathbf{y}_{n+1} = \mathbf{y}_{n+1,\,(p)} + \mathbf{a}\,\frac{h^2}{2!}\left[f(y_{n+1,\,(p)}) - y''_{n+1,\,(p)}\right].$$
- **what:** Corrector step for a 2nd-order ODE: add correction proportional to local residual.
- **symbols:** $\mathbf{a}$ - corrector column vector (Gear coefficients, table in checks.md); bracket - residual = true accel minus predicted 2nd derivative.

<!-- eq:2.10 -->
$$\mathbf{y}_n(F) \equiv [y_n, \; hy'_n, \; h^2 y''_n/2, \; \dots, \; h^2 y''_{n-k+3}/2]^{T}.$$
- **what:** Force (F) representation column vector: position, velocity, then present and past accelerations. Required for constraint dynamics.
- **symbols:** $\mathbf{y}_n(F)$ - F-representation state; $y''_{n-j}$ - past accelerations (forces).

<!-- eq:2.11a -->
$$y_{n+1} = 2y_n - y_{n-1} + h^2 y''_n,$$
- **what:** Verlet position update.
- **symbols:** $y_{n\pm1}$ - positions at neighboring steps; $y''_n$ - acceleration at step n.

<!-- eq:2.11b -->
$$y'_n = (y_{n+1} - y_{n-1})/2h.$$
- **what:** Verlet central-difference velocity.
- **symbols:** $y'_n$ - velocity at step n.

<!-- eq:2.12 -->
$$\mathbf{y}_n(V) \equiv [y_n, \; h^2 y''_n/2, \; y_{n-1}]^{T}.$$
- **what:** Verlet (V) representation column vector (3-value: current position, scaled accel, previous position).
- **symbols:** as above.

<!-- eq:3.1 -->
$$\text{SHAKE}(y_1, y_2, y_3).$$
- **what:** SHAKE operator: reset non-constraint positions $y_2$ to constrained positions $y_3$, with displacement direction $(y_3-y_2)$ set by reference positions $y_1$. Iterative to relative tolerance *tol*.
- **symbols:** $y_1$ - reference positions; $y_2$ - input (non-constraint) positions; $y_3$ - output constrained positions.

<!-- eq:3.2 -->
$$\mathbf{y}_{n+1,\,(p)} = \mathbf{B}\mathbf{y}_n,$$
- **what:** Predictor step in the F-representation (constraint dynamics).
- **symbols:** $\mathbf{B}$ - F-representation predictor matrix (eq A.2a, numeric in checks.md).

<!-- eq:3.3 -->
$$f_{\text{tot}} = f_{\text{free}} + f_{\text{constr}}.$$
- **what:** Total force decomposed into free (interaction, excluding constrained DOF) plus constraint force.
- **symbols:** $f_{\text{tot}}$ - total force; $f_{\text{free}}$ - potential force excluding constrained DOF; $f_{\text{constr}}$ - unknown constraint force.

<!-- eq:3.4 -->
$$\mathbf{y}_{n+1,\,\text{(tot)}} = \mathbf{y}_{n+1,\,\text{(p)}} + \mathbf{b}\,\tfrac{1}{2} h^2 \left[f_{\text{tot}}(y_{n+1,\,\text{(p)}}) - y''_{n+1,\,\text{(p)}}\right].$$
- **what:** Full corrector including constraint force (cannot be applied directly since $f_{\text{tot}}$ unknown; used in final step 6).
- **symbols:** $\mathbf{b}$ - F-rep corrector vector (eq A.7); $b_0$ - its first component.

<!-- eq:3.5 -->
$$\mathbf{y}_{n+1,\,\text{(free)}} = \mathbf{y}_{n+1,\,(p)} + \mathbf{b}\,\tfrac{1}{2} h^2 \left[f_{\text{free}}(y_{n+1,\,(p)}) - y''_{n+1,\,(p)}\right].$$
- **what:** Corrector using only the free (non-constraint) force; gives positions without constraint effect.
- **symbols:** as above; $(free)$ - free (unconstrained) corrected positions.

<!-- eq:3.6 -->
$$\text{SHAKE}(y_{n+1,\,(p)}, \; y_{n+1,\,\text{(free)}}, \; y_{n+1,\,\text{(tot)}}).$$
- **what:** Apply SHAKE to the free positions to obtain the constraint-satisfying total positions.
- **symbols:** arguments = (reference, input free, output total) per eq 3.1 convention.

<!-- eq:3.7 -->
$$f_{\text{constr}}(y_{n+1,\,(p)}) = \left[y_{n+1,\,\text{(tot)}} - y_{n+1,\,\text{(free)}}\right] / (b_0 \tfrac{1}{2} h^2).$$
- **what:** Recover the constraint force from the position displacement produced by SHAKE.
- **symbols:** $b_0$ - first component of corrector vector $\mathbf{b}$; displacement over $b_0 h^2/2$ gives the force.

<!-- eq:3.8 -->
$$y_2 = y_{n+1,\,(p)} - b_0\,\tfrac{1}{2} h^2 \, y''_{n+1,\,(p)}$$
- **what:** Stability modification (a): reference positions to be re-SHAKEd so eq 3.7 gives a correct constraint force.
- **symbols:** as above.

<!-- eq:3.9 -->
$$\text{SHAKE}(y_n, y_2, y_3).$$
- **what:** Extra SHAKE (step 2a) ensuring $y_2$ (eq 3.8) satisfies the constraints.
- **symbols:** $y_n$ - reference; $y_2$ - eq 3.8 input; $y_3$ - reset output used in step 3.

<!-- eq:3.10 -->
$$y'_{n+1} = \left\{ B_{11} y_{n+1} + \sum_{\substack{i=0\\i\neq 1}}^{k-1} (B_{1i} - B_{11} B_{0i}) y_{n,i} + (b_1 - B_{11} b_0)\,\tfrac{1}{2} h^2 \left[f_{\text{tot}}(y_{n+1,\,(p)}) - y''_{n+1,\,(p)}\right] \right\} / h.$$
- **what:** Stability modification (b): recompute velocity from $y_{n+1}$ and $y''_{n+1}$ to avoid propagating errors in past velocities.
- **symbols:** $B_{ij}$ - elements of the F-rep predictor matrix $\mathbf{B}$; $b_0, b_1$ - first two components of corrector $\mathbf{b}$; $y_{n,i}$ - ith component of $\mathbf{y}_n$.

<!-- eq:A.2a -->
$$\mathbf{B} = \mathbf{T}\mathbf{A}\mathbf{T}^{-1},$$
- **what:** Transform the N-representation predictor matrix to the F-representation.
- **symbols:** $\mathbf{T}$ - N->F transformation matrix (eq A.4); $\mathbf{A}$ - Pascal predictor (eq A.1).

<!-- eq:A.2b -->
$$\mathbf{b} = \mathbf{T}\mathbf{a}.$$
- **what:** Transform the N-representation corrector vector to the F-representation.
- **symbols:** $\mathbf{a}$ - N-rep Gear corrector; $\mathbf{b}$ - F-rep corrector.

<!-- eq:A.3 -->
$$\mathbf{y}_n(F) = \mathbf{T}\mathbf{y}_n(N).$$
- **what:** Definition of the transformation matrix T (maps N-state to F-state).
- **symbols:** as above.

<!-- eq:A.7 -->
$$b_i = a_i, \quad i = 0, 1, 2; \qquad b_i = 0, \quad i > 2.$$
- **what:** F-representation corrector vector: equals the Gear corrector in its first three entries, zero beyond. Zeros beyond index 2 are what make SHAKE-compatible constraint incorporation possible.
- **symbols:** $b_i$ - ith F-rep corrector component; $a_i$ - ith N-rep Gear corrector component.
