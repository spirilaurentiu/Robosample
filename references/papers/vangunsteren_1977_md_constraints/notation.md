# Notation - van Gunsteren & Berendsen 1977

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| $\mathbf{r}_i$ | position of particle $i$ | length (cartesian), R^3 | cartesian coords used throughout (not generalized) |
| $\mathbf{v}_i = y'$ | velocity of particle $i$ | length/time, R^3 | |
| $\mathbf{F}_i$ | force on particle $i$ | force | $=-\nabla_i V$ (conservative) |
| $m_i$ | mass of particle $i$ | mass | |
| $V$ | potential energy | energy (kcal/mole in results) | Gelin-Karplus empirical function |
| $N$ | number of particles | int | ~1000 typical; 458 extended atoms for BPTI |
| $y, y', y''$ | 3N-dim state, velocity, acceleration vectors | R^{3N} | $y''=f(y)$ is the ODE |
| $f$ | acceleration function (mass-scaled force) | R^{3N} | double sum over particles; 1 evaluation/step |
| $h$ | time step | dimensionless in figures; unit $= 4.889\times10^{-14}$ s | mesh spacing $t_n=nh$; unstable beyond $h\approx0.02$ |
| $n$ | step index | int | |
| $k$ | number of stored values (k-value method) | int, $k\le 8$ | polynomial order $=k-1$; degree = k-1 |
| $\mathbf{A}$ | N-rep predictor matrix | $k\times k$, Pascal triangle | binomial entries $A_{rc}=\binom{c}{r}$ (0-based) |
| $\mathbf{a}$ | N-rep corrector column vector (Gear) | R^k | $a_2=1$ fixed; $a_0,a_1$ optimize accuracy; $a_{\ge2}$ from stability |
| $\mathbf{B}$ | F-rep predictor matrix | $k\times k$ | $\mathbf{B}=\mathbf{T}\mathbf{A}\mathbf{T}^{-1}$ |
| $\mathbf{b}$ | F-rep corrector vector | R^k | $b_i=a_i$ for $i\le2$, else 0 |
| $\mathbf{T}$ | N->F transformation matrix | $k\times k$ | $\mathbf{y}_n(F)=\mathbf{T}\mathbf{y}_n(N)$ |
| $y_{n,i}$ | ith component of state vector $\mathbf{y}_n$ | scalar | |
| $(p)$ | predicted quantity | subscript | after eq 2.8 / 3.2 |
| $(free)$ | corrected using only free force | subscript | |
| $(tot)$ | corrected using total (free+constraint) force | subscript | |
| $f_{\text{free}}$ | interaction force excluding constrained DOF | force | |
| $f_{\text{constr}}$ | constraint force (unknown, solved via SHAKE) | force | |
| *tol* | SHAKE relative constraint accuracy | dimensionless | $10^{-6}$ adequate for $h>0.01$ |
| $E_{\text{tot}}, E_{\text{kin}}$ | total / kinetic energy | kcal/mole | RMS fluctuation is the accuracy diagnostic |

## Representation conventions

- **N (Nordsieck) representation** (eq 2.5, 2.7): stores scaled successive derivatives $h^j y^{(j)}/j!$. Predictor matrix = Pascal triangle. Best for non-constraint dynamics (easy to change $h$ and $k$).
- **F (force) representation** (eq 2.10): stores position, velocity, and present + past accelerations $h^2 y''_{n-j}/2$. Required for constraint dynamics because SHAKE resets positions and the constraint effect must enter through $y''$. The trailing zeros in $\mathbf{b}$ keep past $y''$ values frozen.
- **V (Verlet) representation** (eq 2.12): 3-value, $[y_n, h^2 y''_n/2, y_{n-1}]$. Verlet = 3-value predictor-corrector with the corrector step omitted.

## Units note

- Time-step unit in the results figures is $4.889\times10^{-14}$ s (so $h=0.01 \Rightarrow 4.889\times10^{-16}$ s per step, matching the equilibration step of $9.778\times10^{-16}$ s at $h=0.02$).
- Energies reported in kcal/mole.
- Forces are per-atom cartesian; extended-atom (united-atom) model — hydrogens folded into heavy atoms.
