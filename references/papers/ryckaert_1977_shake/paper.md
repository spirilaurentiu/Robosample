# Numerical Integration of the Cartesian Equations of Motion of a System with Constraints: Molecular Dynamics of n-Alkanes

Ryckaert, Ciccotti, Berendsen (1977), *Journal of Computational Physics* 23, 327-341.

## Abstract

A numerical algorithm integrating the 3N Cartesian equations of motion of a
system of N points subject to holonomic constraints is formulated. The relations
of constraint remain perfectly fulfilled at each step of the trajectory despite
the approximate character of numerical integration. The method is applied to a
molecular dynamics simulation of a liquid of 64 n-butane molecules and compared
to a simulation using generalized coordinates. The method should be useful for
molecular dynamics calculations on large molecules with internal degrees of
freedom.

## 1. Introduction

Molecular dynamics (MD) has been applied to molecular systems with internal
degrees of freedom (N2, H2O, C4H10). In applying MD three problems arise:
(a) the choice of a suitable mechanical model, (b) the derivation of the equations
of motion, and (c) the choice of an efficient numerical integration algorithm.

In polyatomic molecules the fast internal vibrations are usually decoupled from
rotational/translational motion and can be frozen by introducing a certain number
of rigid bonds and angles into the molecular skeleton (N2 becomes a rod, H2O a
rigid triangle, C4H10 a nonrigid solid with one internal rotation). The classical
way to treat such systems is in generalized coordinates (Lagrange-Hamilton
formalism), but as the number of internal degrees of freedom grows it rapidly
becomes hard to write the equations of motion explicitly.

Orban and Ryckaert suggested using Cartesian equations of motion for n-alkane
chains built of n CH2/CH3 point units connected by n-1 rigid bonds and n-2 rigid
angles. These are the Lagrange equations of motion of the first kind, in which the
forces of constraint appear explicitly. The difficulty: constraints are satisfied
exactly at some initial time but drift away in numerical integration of a large
coupled ODE system, because the computed trajectory deviates from the true one as
time proceeds. Reducing the time step to keep constraint discrepancy small is very
inefficient. This paper develops a method still based on Cartesian coordinates but
producing a trajectory in which all constraints are fulfilled exactly at each step,
without any loss of precision.

Consider a system of N interacting points subject to l holonomic constraints:

<!-- eq:1 -->
$$\sigma_k(\{\mathbf{r}\}) = 0 \qquad (k = 1, ..., l).$$

The force on each point splits into a potential force $\mathbf{F}_i$ and a
constraint force $\mathbf{G}_i$ due to all constraints $\sigma_k$ involving the
ith particle:

<!-- eq:2 -->
$$\mathbf{G}_i = -\sum_{k=1}^{l} \lambda_k(t) \, \nabla_i \sigma_k,$$

where the $\{\lambda_k(t)\}$ are a set of l Lagrangian multipliers depending only
on time.

The key idea: instead of solving directly for the top-order derivative
$\{\lambda_k^{(n-2)}(t_0)\}$, replace it by a set of parameters $\{\gamma_k\}$ that
enforce the constraints exactly at time $t_0+\Delta t$. The substitution is shown
to leave the coordinate propagation exact to the algorithm's order $(\Delta t)^m$
while automatically satisfying the constraints. In Section 3 the procedure is
combined with the Verlet algorithm, for which the constraint forces at time $t$
need not be evaluated explicitly; their effect is obtained from the l parameters
$\{\gamma_k\}$, solved by an iterative quadratic process. Section 5 introduces
SHAKE, an alternative iteration that satisfies each constraint in succession and
is suited to large molecules with complicated constraint sets.

## 2. Integration of the Cartesian Equations of Motion (Method of Lagrangian Multipliers)

Let $\mathbf{r}_i$ and $\dot{\mathbf{r}}_i$ be position and velocity of particle i
($i=1,...,N$) and $V(\{\mathbf{r}\})$ the potential energy. The system is subject
to l holonomic constraints (rigid bonds):

<!-- eq:2.1 -->
$$\sigma_k(\{\mathbf{r}(t)\}) \equiv (\mathbf{r}_i(t) - \mathbf{r}_j(t))^2 - d_{ij}^2 = 0 \qquad (k = 1,...,l),$$

where k labels the rigid bond (ij) of length $d_{ij}$. The 3N Lagrangian equations
of motion of the first kind are:

<!-- eq:2.2 -->
$$m_i \ddot{\mathbf{r}}_i = \mathbf{F}_i + \mathbf{G}_i = -\nabla_i V - \sum_{k=1}^{l} \lambda_k \nabla_i \sigma_k \qquad (i = 1, ..., N).$$

The 3N equations (2.2) plus the l constraints (2.1) give 3N+l equations for 3N+l
unknowns ($\{\mathbf{r}(t)\}$, $\{\lambda(t)\}$). Given an initial configuration
satisfying all constraints, the trajectory is uniquely determined.

### Derivation (not implemented) - analytic Taylor solution for the multipliers

Seeking the analytic solution of (2.2), expand the multipliers in a Taylor series:

<!-- eq:2.3 -->
$$\lambda_k(t) = \sum_{n=0}^{\infty} \lambda_k^{(n)}(t_0) \frac{(t - t_0)^n}{n!}.$$

The coordinates then follow as an explicit function of the set
$\{\lambda^{(n)}(t_0)\}$:

<!-- eq:2.4 -->
$$\mathbf{r}_{i}(t, \{\lambda^{(n)}(t_{0})\}) = \mathbf{r}_{i}(t_{0}) + \dot{\mathbf{r}}_{i}(t_{0})(t - t_{0}) + \frac{1}{m_{i}} \sum_{n=2}^{\infty} \frac{(t - t_{0})^{n}}{n!} \left\{ \mathbf{F}_{i}^{(n-2)}(t_{0}) - \sum_{k=1}^{l} \sum_{p=0}^{n-2} \binom{n-2}{p} \lambda_{k}^{(p)}(t_{0}) [(\nabla_{i}\sigma_{k})^{(n-2-p)}]_{t_{0}} \right\},$$

with the time derivatives of the physical force and constraint gradient given by

<!-- eq:2.4b -->
$$\mathbf{F}_i^{(s)}(t_0) \equiv \left[\frac{d^s}{dt^s} \mathbf{F}_i(t)\right]_{t_0} = \sum_{j=1}^{N} \frac{\partial \mathbf{F}_i^{(s-1)}(t_0)}{\partial \mathbf{r}_j} \cdot \dot{\mathbf{r}}_j(t_0),$$

<!-- eq:2.4c -->
$$[(\nabla_i \sigma_k)^{(s)}]_{t_0} \equiv \left[\frac{d^s}{dt^s} (\nabla_i \sigma_k)\right]_{t_0} = \sum_{j=1}^N \nabla_j [(\nabla_i \sigma_k)^{(s-1)}]_{t_0} \cdot \dot{\mathbf{r}}_j(t_0).$$

Because the constraints hold at all times, all time derivatives of (2.1) vanish:
$\sigma_k^{(s+2)}(\{\mathbf{r}(t_0)\}) = 0$. Expanding these derivatives and
substituting (2.4) gives, after isolating the highest-order terms ($\alpha=0$,
$\beta=s$), a linear system for $\lambda_p^{(s)}$:

<!-- eq:2.5b -->
$$\sum_{j=1}^{N} \frac{1}{m_{j}} \left\{ \mathbf{F}_{j}^{(s)}(t_{0}) - \sum_{p=1}^{l} \lambda_{p}^{(s)}(t_{0}) [\nabla_{j}\sigma_{p}]_{t_{0}} \right\} \cdot [\nabla_{j}\sigma_{k}]_{t_{0}} + \mathscr{F}\big(\{\lambda_{p}^{(s-1)}(t_{0}), ..., \lambda_{p}^{(0)}(t_{0})\}, \{\mathbf{r}(t_{0}), \dot{\mathbf{r}}(t_{0})\}\big) = 0,$$

where $\mathscr{F}$ is a known function of lower-order multiplier derivatives and of
positions/momenta at $t_0$. All $\{\lambda_p^{(s)}\}$ ($s=0,1,...$) follow by
successively inverting (2.5b), reducing (2.2) to second-order ODEs.

### The undetermined-parameter substitution

The constraints obtained this way are only fulfilled to the order in the time step
implicit in the algorithm, and the discrepancy grows faster than linearly. If the
algorithm has coordinate error of order $(\Delta t)^{m+1}$ and uses force
derivatives up to order $n-2$, then the constraint residual is at worst

<!-- eq:2.6 -->
$$\sigma(\{\mathbf{r}^{A}(t)\}) = O[(\Delta t)^{(m+1)}],$$

where $\mathbf{r}^A$ denotes the algorithm's values. Instead, compute from (2.5)
only the first $(n-3)$ derivatives of the $\lambda_k$ at $t_0$, and replace the
top derivative $\{\lambda_k^{(n-2)}(t_0)\}$ by parameters $\{\gamma_k\}$ chosen so
that

<!-- eq:2.7 -->
$$\sigma_k(\{\mathbf{r}(t,\gamma_{k})\}) = 0.$$

Because $\lambda_k^{(n-2)}$ enters expansion (2.4) always multiplied by
$(\Delta t)^n$, it follows that

<!-- eq:2.8 -->
$$\lambda_k^{(n-2)}(t_0) - \gamma_k = O[(\Delta t)^{m+1-n}].$$

The difference between the trajectory computed with $\gamma_k$ and with the true
$\lambda_k^{(n-2)}$ is of $O[(\Delta t)^{m+1}]$ - the same order as the algorithm
error - but the constraints are now perfectly fulfilled.

Write the algorithm value of $\mathbf{r}_i(t_0+\Delta t)$ as the sum of a part
independent of $\gamma_k$ and a correction linear in $\gamma_k$:

<!-- eq:2.9 -->
$$\mathbf{r}_{i}(t_0 + \Delta t) = \mathbf{r}_{i}'(t_{0} + \Delta t, \{\lambda_{k}^{(0)}, ..., \lambda_{k}^{(n-3)}\}) + \delta \mathbf{r}_{i}(t_{0} + \Delta t, \{\gamma_{k}\}),$$

where $\mathbf{r}_i'$ is independent of $\gamma_k$ and $\delta\mathbf{r}_i$ is
linear in $\gamma_k$ and depends only on $\{\mathbf{r}(t_0)\}$. Substituting into
(2.1) gives, for the kth constraint on pair (i,j):

<!-- eq:2.10 -->
$$2(\mathbf{r}_{j}' - \mathbf{r}_{i}') \cdot (\delta \mathbf{r}_{j} - \delta \mathbf{r}_{i}) + (\delta \mathbf{r}_{j} - \delta \mathbf{r}_{i})^{2} = d_{ij}^{2} - (\mathbf{r}_{j}' - \mathbf{r}_{i}')^{2}.$$

The displacements from the Taylor expansion (2.4) are

<!-- eq:2.11 -->
$$\delta \mathbf{r}_{i}(t_{0} + \Delta t, \{\gamma_{k}\}) = \frac{1}{m_{i}} \frac{(\Delta t)^{n}}{n!} \sum_{k=1}^{l} \gamma_{k} (\nabla_{i} \sigma_{k})_{t_{0}}.$$

The l equations (2.10), one per constraint, form a matrix equation for the vector
$\{\gamma_k\}$. Remarks: (i) they are nonlinear but solvable efficiently by
iteration (justified as $\Delta t\to 0$ since the nonlinear terms are small);
initiate with $\{\gamma_k=0\}$ in the nonlinear terms and iterate to the physical
solution. (ii) For the Verlet algorithm $n=2$ and no $\lambda$'s need to be
computed at all.

## 3. Numerical Integration Using the Method of Undetermined Parameters

The Verlet algorithm ($h$ is the time step) is

<!-- eq:3.1 -->
$$u(h) = -u(-h) + 2u(0) + h^2\ddot{u}(0) + O(h^4),$$

<!-- eq:3.2 -->
$$\dot{u}(0) = (u(h) - u(-h))/2h + O(h^2).$$

Combining with (2.9), (2.10), (2.11), the update for each point becomes
$\mathbf{r}_i(h) = \mathbf{r}_i'(h) + \delta\mathbf{r}_i(h)$ with

<!-- eq:3.3 -->
$$\mathbf{r}_{i}'(h) = -\mathbf{r}_{i}(-h) + 2\mathbf{r}_{i}(0) + (h^{2}/m_{i}) \mathbf{F}_{i}(0), \qquad \delta \mathbf{r}_{i}(h) = (h^{2}/m_{i}) \sum_{k=1}^{l} \gamma_{k} [\nabla_{i} \sigma_{k}]_{t_{0}},$$

with the $\{\gamma_k\}$ obtained from

<!-- eq:3.4 -->
$$\sigma_k(\{\mathbf{r}_i(h, \{\gamma_k\})\}) = 0 \qquad (k = 1, ..., l).$$

From (2.8),

<!-- eq:3.5 -->
$$\gamma_k(0) = \lambda_k(0) + O(h^2) \qquad (k = 1, ..., l).$$

The trajectory is exact to third order in h (as usual for Verlet), the constraints
are perfectly obeyed, and the $\{\lambda_k(0)\}$ need never be evaluated: the
constraint-force computation at $t=0$ is converted into evaluating the l parameters
$\{\gamma_k\}$.

For rigid bond constraints (2.1),

<!-- eq:3.6 -->
$$(\mathbf{r}_i(t) - \mathbf{r}_j(t))^2 - d_{ij}^2 = 0,$$

substituting the displacement (2.11) into (2.10) yields the quadratic system

<!-- eq:3.7 -->
$$2(\mathbf{r}_{j}'(h) - \mathbf{r}_{i}'(h)) \cdot \left(-h^{2} \sum_{k=1}^{l} \gamma_{k} \left[ \left(\tfrac{\nabla_{j}}{m_{j}} - \tfrac{\nabla_{i}}{m_{i}}\right) \sigma_{k} \right]_{0}\right) + h^{4} \sum_{k=1}^{l} \sum_{k'=1}^{l} \gamma_{k} \gamma_{k'} \left[ \left(\tfrac{\nabla_{j}}{m_{j}} - \tfrac{\nabla_{i}}{m_{i}}\right) \sigma_{k} \right]_{0} \cdot \left[ \left(\tfrac{\nabla_{j}}{m_{j}} - \tfrac{\nabla_{i}}{m_{i}}\right) \sigma_{k'} \right]_{0} = d_{ij}^{2} - (\mathbf{r}_{j}'(h) - \mathbf{r}_{i}'(h))^{2}.$$

For l rigid bonds this is a system of l quadratic equations in $\gamma_k$, solved
by iteration: the nth iterate $\{\gamma_k^{(n)}\}$ is obtained by solving the
linearized equations that substitute $\{\gamma_k^{(n-1)}\}$ into the quadratic
terms of (3.7). The process starts from $\{\gamma_k^{(0)}=0\}$ and converges in a
few steps because the quadratic terms are $\propto h^4$ and small for usual MD time
steps. In the n-butane test, three or four iterations satisfy the constraints to a
relative displacement discrepancy of order $10^{-10}$.

Algorithm (3.3) is not self-starting: two successive configurations
$\{\mathbf{r}(0)\}$, $\{\mathbf{r}(h)\}$ are needed to initiate. Given the initial
state with all constraints satisfied, $\{\mathbf{r}(h)\}$ can be obtained to
$O(h^3)$ by applying the Section-2 method to the second-order Taylor expansion, by
the same procedure as (3.3), (3.4).

## 4. The Case of n-Butane: Comparison with Integration in Generalized Coordinates

Both methods were applied to 64 n-butane molecules. The n-butane molecule is a
four-point system with three rigid C-C bonds between adjacent points and two rigid
C-C-C angles of 109deg28' between adjacent bonds. The system is in a cubic box with
periodic boundary conditions; density $\rho=0.675\,\mathrm{g/cm^3}$ matches liquid
butane and the kinetic energy corresponds to $T\approx 200\,\mathrm{K}$.

**Method 1 - generalized coordinates.** Lagrange equations of the second kind in
generalized coordinates: three center-of-mass coordinates per molecule, three
Eulerian angles for orientation, and one internal-rotation angle about the central
C-C bond. Integrated with the Gear algorithm.

**Method 2 - Cartesian coordinates.** Lagrange equations of the first kind
integrated by the method of undetermined parameters of Section 3.

Let $\{\mathbf{r}^{(1)}(t)\}$ and $\{\mathbf{r}^{(2)}(t)\}$ be the 256-particle
trajectories from methods 1 and 2 from the same initial configuration at $t_0$.
Define the mean deviation between the two trajectories and the mean displacement:

<!-- eq:4.1 -->
$$\langle |\delta \mathbf{r}(t)| \rangle = \left[ \sum_{i=1}^{256} (\mathbf{r}_{i}^{(1)}(t) - \mathbf{r}_{i}^{(2)}(t))^{2}/256 \right]^{1/2},$$

<!-- eq:4.2 -->
$$\langle r(t) \rangle = \left[ \sum_{i=1}^{256} (\mathbf{r}_i^{(1)}(t) - \mathbf{r}_i^{(1)}(0))^2 / 256 \right]^{1/2}.$$

Both calculations used time step $h = 1.95\times 10^{-15}$ s over a period
$T = 1.56\times 10^{-13}$ s (long enough that the normalized center-of-mass
velocity autocorrelation function falls from 1 to 0.3). At $t=t_0+T$:
$\langle|\delta\mathbf{r}(t_0+T)|\rangle = 1.1\times10^{-4}\,\text{Å}$ and
$\langle r(t_0+T)\rangle = 0.6\,\text{Å}$, a relative average discrepancy of
$\approx 2\times10^{-4}$. No total-energy drift was observed on longer times; the
energy oscillates around a stable value with amplitude $\sim 10^{-3}$ of the
kinetic energy. Computer times per integration step: 3.25 s (method 1) and 1.30 s
(method 2) on an IBM 370/168. The Cartesian method is 2.5x cheaper per step, mainly
due to its greater simplicity (a single reference frame and variable set).

## 5. An Alternative Procedure for Coordinate Resetting (SHAKE)

The Section-3 matrix method requires a matrix inversion per step to compute the
$\gamma_k$. An alternative based on a physical picture is valid for algorithms
requiring no force derivatives (e.g. Verlet). Let the kth constraint be

<!-- eq:5.1 -->
$$\sigma_k = (\mathbf{r}_i - \mathbf{r}_j)^2 - d_{ij}^2 = 0.$$

Constraints are satisfied by adding displacement vectors $\delta\mathbf{r}_i$ to
the unconstrained-step positions $\mathbf{r}_i'$ so that $\sigma_k=0$ for
$\mathbf{r}_i = \mathbf{r}_i' + \delta\mathbf{r}_i$. For Verlet the displacement is

<!-- eq:5.2 -->
$$\delta \mathbf{r}_{i} = -\frac{1}{m_{i}} (\Delta t)^{2} \sum_{k=1}^{l} \gamma_{k} [\nabla_{i} \sigma_{k}]_{t_{0}} = -\frac{2(\Delta t)^{2}}{m_{i}} \sum_{k=1}^{l} \gamma_{k} \mathbf{r}_{ij}(t_{0}),$$

where $\mathbf{r}_{ij} = \mathbf{r}_i - \mathbf{r}_j$ is the vector of the kth rigid
bond. Defining $g_{ij} = -2(\Delta t)^2 \gamma_k$, the term
$g_{ij}\mathbf{r}_{ij}(t_0)/m_i$ is the kth constraint's contribution to
$\delta\mathbf{r}_i$, and $g_{ji}\mathbf{r}_{ji}(t_0)/m_j$ its contribution to
$\delta\mathbf{r}_j$, with $g_{ji}=g_{ij}$ (equal and opposite constraint forces
along the bond). Any convergent procedure satisfying all constraints by
displacements of the form

<!-- eq:5.3 -->
$$\delta \mathbf{r}_i = \sum_{j} g_{ij} \mathbf{r}_{ij}(t_0) / m_i,$$

with $g_{ji}=g_{ij}$ for constrained pairs and $g_{ij}=0$ otherwise, gives results
equivalent to the matrix method.

**SHAKE** is an iterative method that considers all constraints in succession. For
particle i and its partner j, correct positions for the kth constraint by

<!-- eq:5.4a -->
$$\delta^k \mathbf{r}_i = g_{ij} \mathbf{r}_{ij}(t_0) / m_i,$$

<!-- eq:5.4b -->
$$\delta^k \mathbf{r}_j = -g_{ij} \mathbf{r}_{ij}(t_0) / m_j.$$

The position $\mathbf{r}_i'$ from the unconstrained step is corrected with
$\sum_k \delta^k\mathbf{r}_i$ and the next particle is considered. For each
constraint a quadratic equation in $g_{ij}$ results. Define the current
inter-particle vector (with partial corrections from earlier constraints applied)
$\mathbf{r}'$, and the incremental correction $\delta\mathbf{r} = \delta^k\mathbf{r}_i - \delta^k\mathbf{r}_j$. Writing $\mathbf{r} = \mathbf{r}_i(t_0)-\mathbf{r}_j(t_0)$,
$g=g_{ij}$, $d=d_{ij}$, the constraint reads

<!-- eq:5.5 -->
$$(\mathbf{r}' + \delta \mathbf{r})^2 - d^2 = 0, \qquad \delta \mathbf{r} = \left(\frac{1}{m_i} + \frac{1}{m_j}\right) g\,\mathbf{r},$$

which yields the quadratic

<!-- eq:5.6 -->
$$2\left(\frac{1}{m_i} + \frac{1}{m_j}\right) g\, (\mathbf{r} \cdot \mathbf{r}') + \left(\frac{1}{m_i} + \frac{1}{m_j}\right)^2 g^2 \mathbf{r}^2 = d^2 - \mathbf{r}'^2.$$

Satisfying the kth constraint partially destroys the previous ones ($k'<k$), so the
sweep is iterated until every constraint residual is below tolerance. The total
correction $\{\delta\mathbf{r}_i\}$ is of the form (5.3) with $g_{ij}$ the sum of
all $g_{ij}$ from successive iterations. For efficiency (5.6) is solved only to
first order per constraint (dropping the $g^2$ term); the iterative nature still
drives each quadratic within tolerance.

Benchmark on a single decane molecule, time step $4\times10^{-16}$ s: both the
initial and unconstrained-step configurations taken from a dynamics run. The
relative performance depends strongly on required accuracy - SHAKE iteration count
and time grow roughly linearly in $-\log(\text{tolerance})$, while the matrix
method spends nearly all time on inversion and reaches high accuracy cheaply. For a
relative tolerance $10^{-7}$ both used the same CPU time (90 ms on a CDC Cyber
74-16); for higher accuracy the matrix method is faster, for lower accuracy SHAKE
is preferred. At tolerance $10^{-10}$ SHAKE was roughly twice as slow as at
$10^{-7}$. SHAKE is generally applicable to very large molecules with hundreds of
atoms; for small molecules in accurate runs the matrix method is preferred.

## 6. Conclusions

An algorithm was described for integrating the Cartesian equations of motion of an
N-point system subject to holonomic constraints, with constraints satisfied
exactly at each step. For liquid n-butane the method has efficiency similar to
integration in generalized coordinates (and 2.5x cheaper CPU per step). The
suitable time step is of the same order as those used in other MD studies. Adding
holonomic constraints to a system of N free interacting points introduces no new
technical problems: the time step is about the same, the CPU cost does not increase
much, and apart from computing the constraint forces the program is equivalent to a
simple-liquid MD program. The method is adequate for molecules with many internal
degrees of freedom, where generalized coordinates become impractical.

## Appendix: Cartesian Forces of Constraint for an n-Point Semirigid Chain (n-Alkane)

The n-alkane model is a semirigid linear chain of n CH2/CH3 groups treated as
interacting points: $(n-1)$ fixed C-C bonds of length $a = 1.53\,\text{Å}$ and
$(n-2)$ fixed C-C-C angles $\theta = 109.28^\circ$. The internal motion is $(n-3)$
rotations about the C-C bonds. Physical forces come from pair interactions
(e.g. Lennard-Jones) between groups on different molecules or on the same molecule
(from the third neighbor onward).

Number the groups $i=1,...,n$ from one terminal group to the other, with positions
$\{\mathbf{r}_i\}$, velocities $\{\dot{\mathbf{r}}_i\}$, physical forces
$\{\mathbf{F}_i\}$, constraint forces $\{\mathbf{G}_i\}$. The constraints are:

<!-- eq:A.1 -->
$$(\mathbf{r}_{i+1} - \mathbf{r}_i)^2 - a^2 = 0 \qquad (i = 1,..., n-1; \text{ bond constraint}),$$

<!-- eq:A.2 -->
$$(\mathbf{r}_{i+2} - \mathbf{r}_i)^2 - b^2 = 0 \qquad (i = 1,..., n-2; \text{ angle constraint}), \quad b = 2a \sin(\theta/2).$$

Label constraints and multipliers: odd labels $c=1,3,5,...,(2n-3)$ are the
successive bond constraints (A.1) between adjacent pairs; even labels
$c=2,4,...,(2n-4)$ are the angle constraints (A.2) between second-neighbor groups.
The equations of motion for each group are

<!-- eq:A.3 -->
$$m_i \ddot{\mathbf{r}}_i = \mathbf{F}_i + \mathbf{G}_i,$$

where

<!-- eq:A.4 -->
$$\mathbf{G}_{i} = -2\lambda_{2i-4}(\mathbf{r}_{i} - \mathbf{r}_{i-2}) - 2\lambda_{2i-3}(\mathbf{r}_{i} - \mathbf{r}_{i-1}) + 2\lambda_{2i-1}(\mathbf{r}_{i+1} - \mathbf{r}_{i}) + 2\lambda_{2i}(\mathbf{r}_{i+2} - \mathbf{r}_{i}).$$

(For the terminal groups $i=1,2,...,n-1,n$, terms involving nonexistent groups are
dropped.) Integrating (A.3) with the Section-3 algorithm gives

<!-- eq:A.5 -->
$$\mathbf{r}_{i}(h) = \mathbf{r}_{i}'(h) + \delta \mathbf{r}_{i}(h),$$

where

<!-- eq:A.5b -->
$$\mathbf{r}_{i}'(h) = -\mathbf{r}_{i}(-h) + 2\mathbf{r}_{i}(0) + (h^{2}/m_{i}) \mathbf{F}_{i}(0), \qquad \delta \mathbf{r}_i(h) = (h^2/m_i) \mathbf{G}_i,$$

with the constraint force in terms of the parameters $\gamma_k$:

<!-- eq:A.6 -->
$$\mathbf{G}_{i} = 2\big[-\gamma_{2i-4}(\mathbf{r}_{i} - \mathbf{r}_{i-2}) - \gamma_{2i-3}(\mathbf{r}_{i} - \mathbf{r}_{i-1}) + \gamma_{2i-1}(\mathbf{r}_{i+1} - \mathbf{r}_{i}) + \gamma_{2i}(\mathbf{r}_{i+2} - \mathbf{r}_{i})\big].$$

The $\{\gamma_k\}$ are $(2n-3)$ parameters corresponding to $\{\lambda_k(0)\}$
through (3.5), solutions of the $(2n-3)$ equations obtained by substituting (A.5)
into (A.1) and (A.2). Written in matrix form (kth line = kth constraint):

<!-- eq:A.7 -->
$$A(\{\mathbf{r}(0)\}, \{\mathbf{r}'(h)\})\,\mathbf{\gamma} = \mathbf{B}(\{\mathbf{r}(0)\}, \{\mathbf{r}'(h)\}, \{\mathbf{\gamma}\}),$$

where A is a $(2n-3)\times(2n-3)$ matrix (elements below), $\gamma$ the column
vector $(\gamma_1,...,\gamma_{2n-3})$, and B the column vector collecting all terms
nonlinear in $\gamma$.

**If k is odd** (constraint between groups $(k+1)/2$ and $(k+1)/2+1$; using $m$ for
the uniform group mass):

<!-- eq:A.8-odd -->
$$B_{k} = a^{2} - (\mathbf{r}'_{(k+3)/2} - \mathbf{r}'_{(k+1)/2})^{2} - \frac{h^{2}}{m^{2}} (\mathbf{G}_{(k+3)/2} - \mathbf{G}_{(k+1)/2})^2,$$

$$A_{k,k-3} = \frac{4h^{2}}{m} (\mathbf{r}'_{(k+3)/2} - \mathbf{r}'_{(k+1)/2}) \cdot (\mathbf{r}_{(k+1)/2} - \mathbf{r}_{(k-3)/2}),$$

$$A_{k,k-2} = \frac{4h^{2}}{m} (\mathbf{r}'_{(k+3)/2} - \mathbf{r}'_{(k+1)/2}) \cdot (\mathbf{r}_{(k+1)/2} - \mathbf{r}_{(k-1)/2}),$$

$$A_{k,k-1} = -\frac{4h^{2}}{m} (\mathbf{r}'_{(k+3)/2} - \mathbf{r}'_{(k+1)/2}) \cdot (\mathbf{r}_{(k+3)/2} - \mathbf{r}_{(k-1)/2}),$$

$$A_{k,k} = -\frac{8h^{2}}{m} (\mathbf{r}'_{(k+3)/2} - \mathbf{r}'_{(k+1)/2}) \cdot (\mathbf{r}_{(k+3)/2} - \mathbf{r}_{(k+1)/2}),$$

$$A_{k,k+1} = -\frac{4h^{2}}{m} (\mathbf{r}'_{(k+3)/2} - \mathbf{r}'_{(k+1)/2}) \cdot (\mathbf{r}_{(k+5)/2} - \mathbf{r}_{(k+1)/2}),$$

$$A_{k,k+2} = \frac{4h^{2}}{m} (\mathbf{r}'_{(k+3)/2} - \mathbf{r}'_{(k+1)/2}) \cdot (\mathbf{r}_{(k+5)/2} - \mathbf{r}_{(k+3)/2}),$$

$$A_{k,k+3} = \frac{4h^{2}}{m} (\mathbf{r}'_{(k+3)/2} - \mathbf{r}'_{(k+1)/2}) \cdot (\mathbf{r}_{(k+7)/2} - \mathbf{r}_{(k+3)/2}).$$

**If k is even** (constraint between groups $k/2$ and $k/2+2$):

<!-- eq:A.8-even -->
$$B_{k} = b^{2} - (\mathbf{r}'_{(k+4)/2} - \mathbf{r}'_{k/2})^{2} - \frac{h^{4}}{m^{2}} (\mathbf{G}_{(k+4)/2} - \mathbf{G}_{k/2})^{2},$$

$$A_{k,k-4} = \frac{4h^{2}}{m} (\mathbf{r}'_{(k+4)/2} - \mathbf{r}'_{k/2}) \cdot (\mathbf{r}_{k/2} - \mathbf{r}_{(k-4)/2}),$$

$$A_{k,k-3} = \frac{4h^{2}}{m} (\mathbf{r}'_{(k+4)/2} - \mathbf{r}'_{k/2}) \cdot (\mathbf{r}_{k/2} - \mathbf{r}_{(k-2)/2}),$$

$$A_{k,k-1} = -\frac{4h^{2}}{m} (\mathbf{r}'_{(k+4)/2} - \mathbf{r}'_{k/2}) \cdot (\mathbf{r}_{(k+2)/2} - \mathbf{r}_{k/2}),$$

$$A_{k,k} = -\frac{8h^{2}}{m} (\mathbf{r}'_{(k+4)/2} - \mathbf{r}'_{k/2}) \cdot (\mathbf{r}_{(k+4)/2} - \mathbf{r}_{k/2}),$$

$$A_{k,k+1} = -\frac{4h^{2}}{m} (\mathbf{r}'_{(k+4)/2} - \mathbf{r}'_{k/2}) \cdot (\mathbf{r}_{(k+4)/2} - \mathbf{r}_{(k+2)/2}),$$

$$A_{k,k+3} = \frac{4h^{2}}{m} (\mathbf{r}'_{(k+4)/2} - \mathbf{r}'_{k/2}) \cdot (\mathbf{r}_{(k+6)/2} - \mathbf{r}_{(k+4)/2}),$$

$$A_{k,k+4} = \frac{4h^{2}}{m} (\mathbf{r}'_{(k+4)/2} - \mathbf{r}'_{k/2}) \cdot (\mathbf{r}_{(k+8)/2} - \mathbf{r}_{(k+4)/2}).$$

<!-- CHECK: OCR gave A_{k,k+3} (even) index as (r_{(k+4)/2}-r_{(k+4)/2}) which is zero/nonsensical; corrected to (r_{(k+6)/2}-r_{(k+4)/2}) by pattern analogy with the odd case. Verify against original PDF. -->
