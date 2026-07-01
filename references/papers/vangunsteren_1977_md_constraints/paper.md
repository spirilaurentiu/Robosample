# Algorithms for macromolecular dynamics and constraint dynamics

W. F. van Gunsteren and H. J. C. Berendsen, *Molecular Physics* **34**(5), 1311-1327 (1977). DOI: 10.1080/00268977700102571

## Abstract

The application of molecular dynamics (MD) to macromolecules is investigated. The protein trypsin inhibitor (BPTI), consisting of 454 united atoms, is used as an example. Different algorithms for integrating the equations of motion are compared, both theoretically and in practice. It is examined to what extent the chain structure of a macromolecule allows a reduction of the computational effort by the introduction of constraints in the dynamics of the chain.

A calculational scheme is proposed by which constraints can be incorporated in predictor-corrector algorithms. The optimum choice of an algorithm depends on the desired accuracy of the solution and on the character of the forces acting on the molecule, i.e. whether these are noisy or not. For nonconstraint dynamics a Gear predictor-corrector algorithm yields the best results, whereas for constraint dynamics the Gear and Verlet algorithms produce comparable results. The application of bond-length constraints reduces the required computer time by a factor of 3. The inclusion of bond-angle constraints is not recommended.

## 1. Introduction

A macromolecule differs from a simple liquid by the more complex nature of the interaction potentials and forces and by the presence of a large number of covalent bonds. Since these bonds represent the highest components in the frequency distribution of the molecular motion, the introduction of constraints for bond lengths and angles is expected to increase computational efficiency.

Two points are investigated:

1. Finding the best presently available algorithm for MD of a macromolecule. The most frequently used integrators are those of Verlet and Gear; Beeman formulated an algorithm claimed to be better than both. The mathematical treatment of numerical solution of ODEs is incomplete, so the answer is partly dependent on the physical system.
2. Finding the best algorithm for dynamics in the presence of constraints. Introducing constraints is only physically correct if (a) the frequency components of the motion along the eliminated degrees of freedom are well separated from the other frequencies, and (b) the coupling between both types of motion is weak.

The evaluation is on the basis of time step, accuracy and computer efficiency.

## 2. Algorithms for MD without constraints

The MD method consists essentially of solving the set of coupled second-order differential equations

<!-- eq:2.1 -->
$$\frac{d^2 \mathbf{r}_i}{dt^2} = \mathbf{F}_i(\mathbf{r}_1, \mathbf{r}_2, \dots, \mathbf{r}_N)/m_i, \quad i = 1, 2, \dots N,$$

governing the classical dynamics of a system of $N$ particles, to find positions $\{\mathbf{r}_i(t)\}$, velocities $\{\mathbf{v}_i(t)\}$, etc. as a function of time. Initial positions $\{\mathbf{r}_i(t_0)\}$ and velocities $\{\mathbf{v}_i(t_0)\}$ must be specified. In systems of interest $N \sim 1000$. The force is derived from a potential $V$:

<!-- eq:2.2 -->
$$\mathbf{F}_i(\mathbf{r}_1, \mathbf{r}_2, \dots \mathbf{r}_N) = -\boldsymbol{\nabla}_i V(\mathbf{r}_1, \mathbf{r}_2, \dots \mathbf{r}_N), \quad i = 1, 2, \dots N.$$

For simplicity the force is taken conservative (position dependent only), but the discussion also holds for non-conservative forces.

Mathematically this is an initial value problem:

<!-- eq:2.3a -->
$$y'' = f(y),$$

<!-- eq:2.3b -->
$$y_0 = y(t_0), \quad y_0' = y'(t_0),$$

where $y, y', y''$ are 3N-dimensional vectors and their derivatives with respect to $t$. The solution is approximated at discrete mesh points, normally equally spaced:

<!-- eq:2.4 -->
$$t_n = nh, \quad n = 0, 1, 2, \dots,$$

where $h$ is the spacing. A $k$-step algorithm uses preceding values of $y$ or its derivatives up to and including $t_{n-k}$.

The most time-consuming part of any solution method is evaluation of the function $f$, because it involves at least a double summation over the $N$ particles. This rules out Runge-Kutta and extrapolation methods (several function evaluations per step). Predictor-corrector methods can be applied directly to a higher-order equation, avoiding the doubling of information from converting to first-order form. The number of corrector iterations in MD is usually taken equal to 1 (one $f$ evaluation per step).

### Representations of multi-value predictor-corrector methods

A $k$-value method uses $k$ previously calculated values of $y$ or its successive derivatives. The representation is determined by which of $y, y', y'', \dots$ are used.

Nordsieck (N-representation) saves and uses

<!-- eq:2.5 -->
$$y_n, \; hy'_n, \; h^2 y''_n/2, \; \dots, \; h^{k-1} y_n^{(k-1)}/(k-1)!$$

The Adams method saves

<!-- eq:2.6 -->
$$y_n, \; hy'_n, \; hy'_{n-1}, \; \dots, \; hy'_{n-k+2}$$

In the N-representation the column vector is

<!-- eq:2.7 -->
$$\mathbf{y}_n(N) \equiv [y_n, \; hy'_n, \; h^2 y''_n/2, \; \dots, \; h^{k-1} y_n^{(k-1)}/(k-1)!]^{T}.$$

The predictor step is

<!-- eq:2.8 -->
$$\mathbf{y}_{n+1,\,(p)} = \mathbf{A}\mathbf{y}_n,$$

where matrix $\mathbf{A}$ is determined by the predictor. The corrector for a second-order ODE is

<!-- eq:2.9 -->
$$\mathbf{y}_{n+1} = \mathbf{y}_{n+1,\,(p)} + \mathbf{a}\,\frac{h^2}{2!}\left[f(y_{n+1,\,(p)}) - y''_{n+1,\,(p)}\right].$$

The second term is the amount by which the ODE is not locally satisfied by $\mathbf{y}_{n+1,\,(p)}$. The column vector $\mathbf{a}$ characterizes the corrector. The corrector may be repeated for a fixed number of iterations or until convergence.

The usual predictor is a Taylor-series extrapolation of polynomial degree $(k-1)$ satisfying the retained values; then $\mathbf{A}$ equals the Pascal triangle (Appendix, eq A.1). Only the last $(k-2)$ coefficients of $\mathbf{a}$ ($a_2, a_3, \dots$) are determined by stability; the first two ($a_0, a_1$) are chosen to optimize accuracy (per Gear).

### Beeman and Verlet as special cases

The Beeman algorithm is written in the force (F) representation

<!-- eq:2.10 -->
$$\mathbf{y}_n(F) \equiv [y_n, \; hy'_n, \; h^2 y''_n/2, \; \dots, \; h^2 y''_{n-k+3}/2]^{T}.$$

Transforming to the N-representation shows it differs from Gear only in $(a_0, a_1)$. The Gear coefficients produce results at least twice as accurate as Beeman's.

The Verlet algorithm consists of

<!-- eq:2.11a -->
$$y_{n+1} = 2y_n - y_{n-1} + h^2 y''_n,$$

<!-- eq:2.11b -->
$$y'_n = (y_{n+1} - y_{n-1})/2h.$$

In the representation

<!-- eq:2.12 -->
$$\mathbf{y}_n(V) \equiv [y_n, \; h^2 y''_n/2, \; y_{n-1}]^{T},$$

transforming eq (2.11a) to the N-representation shows Verlet is equivalent to a 3-value predictor-corrector algorithm *without* evaluating the corrector.

Changing the time-step size or the degree ($k$) of the algorithm is most easily done in the N-representation (only last-step information; only $\mathbf{a}$ changes, $\mathbf{A}$ fixed). Conclusion: the best MD algorithm without constraints is a $k$-value predictor-corrector in the N-representation, predictor from Taylor expansion, corrector from Gear stability/accuracy.

## 3. Algorithms for MD with constraints

A MD calculation may be sped up by reducing degrees of freedom, e.g. eliminating bond-stretching vibrations via constraint dynamics (bond lengths kept fixed). Generalized coordinates / Lagrangian equations are impractical for macromolecules; cartesian coordinates must be used, also in the presence of constraints.

Two methods (Ryckaert, Ciccotti, Berendsen 1977) integrate cartesian equations of motion under holonomic constraints. The procedure suitable for large systems is **SHAKE**: at each step the constraints are satisfied by adding displacement vectors to the position vectors that result from a non-constraint time step. Because constraints are interdependent, SHAKE iterates over all constraints in succession to relative accuracy *tol*. Notation:

<!-- eq:3.1 -->
$$\text{SHAKE}(y_1, y_2, y_3).$$

Positions $y_2$ from the non-constraint step are reset to constrained positions $y_3$; the direction of the displacement vectors $(y_3 - y_2)$ is determined by reference positions $y_1$.

Until now SHAKE was used only with Verlet (which uses only positions $(y_n, y_{n-1})$). The scheme below incorporates the effect of the constraint forces into the values of $y''$, enabling arbitrary predictor-corrector algorithms. The algorithm must be written in a representation containing (besides one value of $y'$) only values of $y$ and/or $y''$ at the mesh points — so constraint dynamics cannot be done in the N-representation. The natural choice is the F-representation (2.10). The predictor matrix $\mathbf{A}$ and corrector $\mathbf{a}$ are transformed from N to F giving $\mathbf{B}$ and $\mathbf{b}$ (Appendix A.2, eqs A.6, A.7).

### Computational scheme for a MD step with constraints

1. Predictor step:

<!-- eq:3.2 -->
$$\mathbf{y}_{n+1,\,(p)} = \mathbf{B}\mathbf{y}_n,$$

giving $y_{n+1,(p)}$ and $y''_{n+1,(p)}$.

2. The total force splits into constraint and free parts:

<!-- eq:3.3 -->
$$f_{\text{tot}} = f_{\text{free}} + f_{\text{constr}}.$$

$f_{\text{free}}$ is derived from the interaction potential excluding interaction along the constrained degrees of freedom; $f_{\text{constr}}$ are the unknown constraint forces. The full corrector cannot be applied directly because $f_{\text{tot}}$ is unknown:

<!-- eq:3.4 -->
$$\mathbf{y}_{n+1,\,\text{(tot)}} = \mathbf{y}_{n+1,\,\text{(p)}} + \mathbf{b}\,\tfrac{1}{2} h^2 \left[f_{\text{tot}}(y_{n+1,\,\text{(p)}}) - y''_{n+1,\,\text{(p)}}\right].$$

2a. (Stability modification a) Ensure that the reference positions satisfy the constraints. Compute

<!-- eq:3.8 -->
$$y_2 = y_{n+1,\,(p)} - b_0\,\tfrac{1}{2} h^2 \, y''_{n+1,\,(p)}$$

and reset via

<!-- eq:3.9 -->
$$\text{SHAKE}(y_n, y_2, y_3).$$

The resulting positions $y_3$ are used in step 3.

3. Apply the corrector without constraints:

<!-- eq:3.5 -->
$$\mathbf{y}_{n+1,\,\text{(free)}} = \mathbf{y}_{n+1,\,(p)} + \mathbf{b}\,\tfrac{1}{2} h^2 \left[f_{\text{free}}(y_{n+1,\,(p)}) - y''_{n+1,\,(p)}\right].$$

This yields positions at $t_{n+1}$ without the effect of constraint forces.

4. Incorporate the constraint effect by resetting with SHAKE:

<!-- eq:3.6 -->
$$\text{SHAKE}(y_{n+1,\,(p)}, \; y_{n+1,\,\text{(free)}}, \; y_{n+1,\,\text{(tot)}}).$$

5. The constraint forces are then obtained from

<!-- eq:3.7 -->
$$f_{\text{constr}}(y_{n+1,\,(p)}) = \left[y_{n+1,\,\text{(tot)}} - y_{n+1,\,\text{(free)}}\right] / (b_0 \tfrac{1}{2} h^2).$$

(Obtained by subtracting the first row of eq 3.5 from that of eq 3.4.) The total force $f_{\text{tot}}$ follows from eq (3.3).

6. Use $f_{\text{tot}}$ in the corrector (3.4).

### Stability modification (b): velocity update

When velocities $y'_{n+1}$ are calculated from eq (3.4), errors in preceding $y'_n, y'_{n-1}$ may propagate. Instead compute $y'_{n+1}$ from $y_{n+1}$ and $y''_{n+1} = f_{\text{tot}}(y_{n+1,(p)})$ using

<!-- eq:3.10 -->
$$y'_{n+1} = \left\{ B_{11} y_{n+1} + \sum_{\substack{i=0\\i\neq 1}}^{k-1} (B_{1i} - B_{11} B_{0i}) y_{n,i} + (b_1 - B_{11} b_0)\,\tfrac{1}{2} h^2 \left[f_{\text{tot}}(y_{n+1,\,(p)}) - y''_{n+1,\,(p)}\right] \right\} / h.$$

### Derivation (not implemented): elimination of $y'_n$

Equation (3.10) is obtained by eliminating $y'_n$ from the two corrector rows

<!-- eq:3.11a -->
$$y_{n+1} = \sum_{i=0}^{k-1} B_{0i} y_{n,i} + b_0\,\tfrac{1}{2} h^2 \left[f_{\text{tot}}(y_{n+1,\,(p)}) - y''_{n+1,\,(p)}\right]$$

<!-- eq:3.11b -->
$$h\,y'_{n+1} = \sum_{i=0}^{k-1} B_{1i} y_{n,i} + b_1\,\tfrac{1}{2} h^2 \left[f_{\text{tot}}(y_{n+1,\,(p)}) - y''_{n+1,\,(p)}\right],$$

where $y_{n,i}$ denotes the $i$th component of the vector $\mathbf{y}_n$, obtained by inserting (3.2) into (3.4).

The extra computer time versus non-constraint MD is almost entirely determined by SHAKE and depends directly on *tol*. Conclusion: the best algorithm for MD with constraints is a $k$-value predictor-corrector in the F-representation (Taylor predictor, Gear corrector) with the SHAKE scheme of this section.

## 4. Macromolecular dynamics with and without constraints

The MD method was applied to BPTI (bovine pancreatic trypsin inhibitor, 58 residues), a single molecule in vacuum (solvent neglected).

### 4.1. Parameters of the calculation

The empirical interaction function of Gelin and Karplus was used, containing bond stretching, bond-angle bending, dihedral-angle twisting, hydrogen bonds, non-bonded (Van der Waals) and electrostatic (Coulomb) contributions. The extended-atom concept is used (hydrogens incorporated into the heavy atoms they bind). Number of extended atoms = 458 (454 in BPTI + 4 hydrogen-bonded water molecules), giving 1374 cartesian degrees of freedom. Number of possible bond-length constraints = 468; possible angle constraints = 626.

Starting configuration: from the X-ray structure with all velocities zero, a 100-step run with the 6-value Gear predictor-corrector, time step $9.778 \times 10^{-16}$ s, reached ~150 K. All velocities were then multiplied by 1.74 and 100 more steps taken, approaching equilibrium at ~300 K. Analysis in comparison runs started after 10 steps.

### 4.2. Results

Algorithms compared over 100-step runs. Diagnostics: root mean square fluctuation of total energy relative to RMS fluctuation of kinetic energy, plus the drift of total energy per step (least-squares linear fit vs time). These depend on the algorithm, $k$, $h$, and (for constraints) *tol*.

#### 4.2.1. Non-constraint dynamics

Verlet requires storage of 3 vectors, 1.2 s CPU/step on a CDC Cyber 74-16. Gear requires $k+1$ vectors, 1.3 s/step ($k=4$) to 1.5 s/step ($k=8$).

For $h < 0.02$ (time-step units $4.889 \times 10^{-14}$ s) the total-energy fluctuation is more than a factor 100 smaller than the kinetic-energy fluctuation. Beyond $h = 0.02$ the algorithms become successively unstable. Below $h \approx 0.03$ Gear is more accurate than Verlet; above it Verlet is better (Verlet uses a degree-2 polynomial, Gear uses degrees 3-7; higher degree blows up at smaller time step). Maximum small-step accuracy is reached at $k = 7$, implying the force is rather predictable (weakly damped harmonic bond-stretching vibrations).

Drift per step behaves like the fluctuation. At $h = 0.01$: drift $\approx 10^{-6}$ kcal/mole (Gear) vs $\approx 10^{-3}$ kcal/mole (Verlet).

For simple liquids (more heavily damped, less predictable motion) a large $k$-value may be pointless; Verlet may match Gear. Same applies when forces are noisy (truncation, tabulation, stochastic components).

#### 4.2.2. Constraint dynamics

Constrained bond lengths are set to the bond-stretch potential minimum. SHAKE requires one extra stored vector. *tol* should cause an error smaller than the algorithm's truncation error; for relevant $h > 0.01$, $tol = 10^{-6}$ is small enough (energy fluctuation changes < 1% on decreasing *tol*).

With Verlet + constraints ($tol = 10^{-6}$): for the same accuracy the time step can be taken 4 times as large as in non-constraint dynamics, while CPU time per step rises from 1.2 s to ~1.6 s (constraint dynamics ~1.3x slower per step). Net speedup ≈ factor 3.

For Gear with constraints, larger $k$ yields no improvement over $k = 4$; for $h > 0.02$ even Verlet beats Gear — the force has a strong random character once the high-frequency bond-stretching forces are removed.

Bond-angle constraints: SHAKE becomes very time-consuming for the enlarged constraint set, and no time-step increase is allowed because bond-angle vibration frequencies are not well separated from other vibrations. Bond-angle constraints are not worthwhile.

## 5. Summary and conclusions

- The best presently available integrator for MD is a multi-value ($k$-value) predictor-corrector with Gear stability/accuracy parameters. Non-constraint: use the N-representation. Constraint: use the F-representation. A scheme was proposed to incorporate constraints into a multi-value predictor-corrector.
- Non-constraint BPTI: most accurate results from Gear at large $k$ (~7). Lower $k$ allows larger time step if less accuracy is needed.
- Constraint BPTI (fixed bond lengths): Verlet and Gear give comparable accuracy; Gear slightly better at small $h$, Verlet at large $h$; Gear with larger $k$ gives no improvement over $k=4$. Bond-length constraints speed up MD by ~factor 3.
- Bond-angle constraints do not pay, physically or computationally.

Practical procedure: (1) decide on bond-length constraints; (2) choose Gear (no constraints) or Gear/Verlet (constraints, depending on step); (3) choose an accuracy (upper limit on RMS total-energy fluctuation); (4) for constraints, find the largest *tol* compatible with that accuracy; (5) determine optimum $h$ and $k$ from test calculations; (6) run production.

## Appendix

### A.1. The $k$-value predictor-corrector in the N-representation

In the N-representation the predictor matrix $\mathbf{A}$ is the Pascal triangle ($k \le 8$):

<!-- eq:A.1 -->
$$\mathbf{A} = \begin{bmatrix} 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 \\ 0 & 1 & 2 & 3 & 4 & 5 & 6 & 7 \\ 0 & 0 & 1 & 3 & 6 & 10 & 15 & 21 \\ 0 & 0 & 0 & 1 & 4 & 10 & 20 & 35 \\ 0 & 0 & 0 & 0 & 1 & 5 & 15 & 35 \\ 0 & 0 & 0 & 0 & 0 & 1 & 6 & 21 \\ 0 & 0 & 0 & 0 & 0 & 0 & 1 & 7 \\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \end{bmatrix}$$

<!-- CHECK: the printed A shows only nonzero upper-triangular Pascal entries; the lower-triangle zeros and last rows/columns are reconstructed from "A = Pascal triangle (k<8)". Binomial C(col,row) with 0-based indexing. -->

The corrector column vector $\mathbf{a}$ ensuring optimum stability/accuracy for a second-order ODE (2.3) is tabulated below (from Gear, p. 154). See `checks.md` for the numeric table.

### A.2. The $k$-value predictor-corrector in the F-representation

$\mathbf{B}$ and $\mathbf{b}$ are obtained from $\mathbf{A}$ and $\mathbf{a}$ via the transformation $\mathbf{T}$:

<!-- eq:A.2a -->
$$\mathbf{B} = \mathbf{T}\mathbf{A}\mathbf{T}^{-1},$$

<!-- eq:A.2b -->
$$\mathbf{b} = \mathbf{T}\mathbf{a}.$$

$\mathbf{T}$ follows from

<!-- eq:A.3 -->
$$\mathbf{y}_n(F) = \mathbf{T}\mathbf{y}_n(N).$$

The transformation matrix, its inverses per $k$, and the resulting $\mathbf{B}$ matrices are concrete numeric fixtures — see `checks.md`.

For $\mathbf{b}$:

<!-- eq:A.7 -->
$$b_i = a_i, \quad i = 0, 1, 2; \qquad b_i = 0, \quad i > 2.$$

The last $(k-3)$ rows of $\mathbf{B}$ are zero except for a single 1, and the last $(k-3)$ coefficients of $\mathbf{b}$ are zero. This is fortunate: otherwise $y''$ calculated at previous steps would change at later steps, making incorporation of SHAKE impossible.
