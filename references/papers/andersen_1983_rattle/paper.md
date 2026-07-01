# RATTLE: A "Velocity" Version of the SHAKE Algorithm for Molecular Dynamics Calculations

**Author:** Hans C. Andersen (Department of Chemistry, Stanford University), 1983.
**Venue:** Journal of Computational Physics 52, 24-34.

## Abstract

An algorithm, called RATTLE, for integrating the equations of motion in molecular
dynamics calculations for molecular models with internal constraints. It is similar
to SHAKE, one of the standard methods for such calculations. RATTLE calculates the
positions and velocities at the next time step from the positions and velocities at
the present time step, without requiring information about earlier history. Like
SHAKE, it is based on the Verlet algorithm and retains the simplicity of using
Cartesian coordinates for each atom to describe the configuration of a molecule with
internal constraints. RATTLE guarantees that the coordinates AND velocities of the
atoms satisfy the internal constraints at each time step. RATTLE has two advantages
over SHAKE: (1) on computers of fixed precision it is of higher precision than SHAKE,
and (2) since it deals directly with the velocities, it is easier to modify for use
with constant-temperature and constant-pressure MD methods and with nonequilibrium
MD methods that rescale atomic velocities.

## I. Introduction

Various models are used in molecular dynamics for molecular fluids: vibrating models,
rigid models, and flexible models with internal constraints. For rigid molecules,
Euler angles or quaternions can represent rotational degrees of freedom. An
alternative is the SHAKE method of Ryckaert et al., which retains the simplicity of
Cartesian coordinates and avoids the complications of Euler angles and quaternions
while incorporating the effects of constrained geometry. SHAKE can also be used for
flexible molecules with internal constraints.

SHAKE is based on the Verlet algorithm (the explicit central difference method). The
Verlet algorithm as usually implemented has several drawbacks: the velocities are not
among the integration variables and can be obtained only with extra effort/storage;
it is difficult to start with chosen coordinates and velocities; it is difficult to
change the time step and continue; and precision is easily lost because quantities of
very different magnitude are added.

Because velocities are not integration variables, several modern MD methods are hard
to implement on top of Verlet/SHAKE: stochastic collisions for constant-temperature
MD; constant-pressure MD (whose equations of motion contain the velocities, which
Verlet cannot solve); and velocity rescaling in nonequilibrium MD.

For atomic fluids and vibrating models, these drawbacks are overcome by the "velocity"
version of the Verlet algorithm (velocity Verlet), which explicitly computes the
velocities as part of solving the equations of motion. SHAKE cannot be expressed in
an analogous way directly. This paper presents a generalization of SHAKE with the
advantages of velocity Verlet, called RATTLE: the positions and velocities at time
$t$ are used to compute the positions and velocities at $t + h$. At each time, the
positions satisfy the internal constraints to within any desired accuracy, and (unlike
SHAKE) the velocities also satisfy the (time-derivative) constraints.

## II. Derivation of the Algorithm RATTLE

Consider a differential equation of the form

<!-- eq:2.0a -->
$$ \ddot{r}(t) = f[r], $$

where $r$ denotes the set of Cartesian coordinates specifying the configuration of the
system. The Verlet algorithm for solving this equation is

<!-- eq:2.1 -->
$$ r(t+h) = 2r(t) - r(t-h) + h^2 f[r(t)]. $$

One starts with $r(0)$ and $r(h)$ and calculates each succeeding value from the two
preceding values. Eq. (2.1) is correct except for local errors of order $h^3$; over a
finite interval the global error is of order $h^2$.

One problem is that the calculation adds a term of order $h^2$ to terms of order $h^0$,
so only a few significant figures of the force $f[r(t)]$ are utilized and precision can
be lost. The velocity version below avoids this.

A second problem is that velocities do not appear explicitly. To evaluate the velocity
in Verlet, the following approximation can be used:

<!-- eq:2.2 -->
$$ \dot{r}(t) = [r(t+h) - r(t-h)]/2h. $$

This makes an error of order $h^2$. Note the velocity at time $t$ can be obtained only
after $r(t+h)$ has been obtained, which makes stochastic collisions and
velocity-dependent accelerations difficult.

The velocity version of the Verlet algorithm eliminates both problems:

<!-- eq:2.3 -->
$$ r(t+h) = r(t) + h\dot{r}(t) + h^2 f[r(t)]/2, $$

<!-- eq:2.4 -->
$$ \dot{r}(t+h) = \dot{r}(t) + h[f[r(t)] + f[r(t+h)]]/2. $$

These are equivalent to the Verlet algorithm (2.1) plus (2.2). The local error in each
is of order $h^3$; the global error is of order $h^2$, the same as the usual Verlet
form. Velocity Verlet explicitly involves the velocities, so stochastic collisions and
velocity rescaling at time $t$ are done by changing $\dot{r}(t)$ before computing the
$t+h$ quantities. Storage requirement is $3N$ locations for $N$ degrees of freedom.

### Constrained equations of motion

When internal constraints are present (e.g. fixed internuclear distances or bond
angles), the equations of motion become

<!-- eq:2.0b -->
$$ \ddot{r}(t) = f[r(t)] + g[r(t), \dot{r}(t)], $$

where $f$ includes the physical forces and $g$ includes the constraint forces. The
constraint forces $g$ depend on all details of the mechanical state, contain
time-dependent Lagrange multipliers, and their functional form depends on the nature
of the constraints.

The Verlet algorithm for this equation is

<!-- eq:2.0c -->
$$ r(t+h) = 2r(t) - r(t-h) + h^{2}[f[r(t)] + g[r(t), \dot{r}(t)]]. $$

Even if the exact $g$ were known, the intramolecular constraints would eventually be
violated because the algorithm is not exact. Ryckaert et al. noticed this could be
solved by using not the exactly correct $g$ but an approximation for $g$ that requires
$r(t+h)$ to satisfy the constraints exactly (or to desired accuracy); this
approximation makes errors of the same order as the local error of Verlet. The SHAKE
algorithm is thus

<!-- eq:2.0d -->
$$ r(t+h) = 2r(t) - r(t-h) + h^{2}[f[r(t)] + g_{s}(t)], $$

where $g_s$ is the SHAKE approximation for $g$.

### Why SHAKE has no simple velocity form

A naive "velocity version" of SHAKE would give

<!-- eq:2.5 -->
$$ r(t+h) = r(t) + h\dot{r}(t) + h^{2}[f[r(t)] + g[r(t), \dot{r}(t)]]/2, $$

<!-- eq:2.6 -->
$$ \dot{r}(t+h) = \dot{r}(t) + h[f[r(t)] + g[r(t), \dot{r}(t)] + f[r(t+h)] + g[r(t+h), \dot{r}(t+h)]]/2. $$

Using (2.5) one could compute $r(t+h)$ by replacing $g[r(t),\dot{r}(t)]$ with an
approximation that makes $r(t+h)$ satisfy the constraints. But then (2.6) cannot be
used to get $\dot{r}(t+h)$ because there is no way to evaluate the second $g$ appearing
there: by (2.5) we need $g(t)$ before computing anything, but by (2.6) we need
$g(t+h)$ before computing $\dot{r}(t+h)$. This is inconsistent with a simple iterative
scheme.

### The RATTLE resolution

There is no requirement that the same approximation for $g$ be used in the position
equation as in the velocity equation. Using (2.5), compute $r(t+h)$ by choosing the
$g$ so that $r(t+h)$ satisfies the constraints exactly:

<!-- eq:2.7 -->
$$ r(t+h) = r(t) + h\dot{r}(t) + h^{2}[f[r(t)] + g_{RR}(t)]/2. $$

Knowing $r(t+h)$, compute $f[r(t+h)]$. Then using (2.6) compute $\dot{r}(t+h)$ by
choosing the second $g$ so that the resulting $\dot{r}(t+h)$ satisfies the time
derivatives of the constraints exactly:

<!-- eq:2.8 -->
$$ \dot{r}(t+h) = \dot{r}(t) + h[f[r(t)] + g_{RR}(t) + f[r(t+h)] + g_{RV}(t)]/2. $$

RATTLE makes two separate approximations, $g_{RR}$ and $g_{RV}$, for the constraint
forces. As a result, both the positions and the velocities can be required to satisfy
the constraints.

## III. Properties of the Algorithm RATTLE

RATTLE, defined by Eqs. (2.7) and (2.8) and the conditions for choosing $g_{RR}(t)$ and
$g_{RV}(t)$, has these important characteristics:

1. The positions and velocities at one time are used to calculate positions and
   velocities at the next time, without using information from previous times.
2. The coordinates at each time satisfy the intramolecular constraints.
3. The velocities at each time satisfy the (time-derivative) constraints.
4. The precision is comparable to velocity Verlet rather than the original Verlet form.

The algorithm can be initiated by choosing positions and velocities at one time, so it
is easy to implement stochastic collisions and changes in the time step magnitude.

The global error of RATTLE is of order $h^2$ for small $h$, as in Verlet for
unconstrained dynamics and in SHAKE (proof in Appendix B). Consequently, energy along
a trajectory is conserved only to within errors of order $h^2$.

Test calculations were performed on a pendulum swinging in a gravitational potential in
two dimensions, using both RATTLE and SHAKE. The numerical work verifies that the
global errors in energy, coordinates, and velocities are quadratic in the time step.
The coefficients of the quadratic errors are similar in magnitude for the two
algorithms. The method has also been programmed for MD calculations for water and
appears to work satisfactorily.

RATTLE requires two sets of constraint forces ($g_{RR}$ and $g_{RV}$) rather than one,
so this part of the computation takes longer than SHAKE. However, since intermolecular
force calculation typically dominates the compute time, RATTLE is in practice as
efficient as SHAKE. Both require the same storage, approximately $3N$ floating point
numbers for $N$ degrees of freedom.

## Appendix A: Details of RATTLE

We restrict attention to constraints requiring pairs of mass points to remain a fixed
distance apart. If $i$ and $j$ are a constrained pair, define

<!-- eq:A0 -->
$$ \sigma_{ij}(\{\mathbf{r}(t)\}) \equiv [\mathbf{r}_i(t) - \mathbf{r}_j(t)]^2 - d_{ij}^2, $$

where $\mathbf{r}_i$ is the position of atom $i$, $m_i$ its mass, and $d_{ij}$ the fixed
distance between atoms $i$ and $j$. The constraint is

<!-- eq:A1 -->
$$ \sigma_{ij}(\{\mathbf{r}(t)\}) = 0. $$

The time derivatives of the constraints give constraints on the velocities:

<!-- eq:A2 -->
$$ \left[\dot{\mathbf{r}}_{i}(t) - \dot{\mathbf{r}}_{j}(t)\right] \cdot \left[\mathbf{r}_{i}(t) - \mathbf{r}_{j}(t)\right] = 0. $$

The equations for constrained dynamics are

<!-- eq:A2b -->
$$ m_i \ddot{\mathbf{r}}_i = \mathbf{F}_i + \mathbf{G}_i, $$

where $\mathbf{F}_i$ is the force due to intermolecular and non-constraint
intramolecular interactions and $\mathbf{G}_i$ is the constraint force on atom $i$,
given by

<!-- eq:A2c -->
$$ \mathbf{G}_i = -\sum_j{}' \lambda_{ij}(t) \nabla_i \sigma_{ij}, $$

where the prime denotes a sum over only those atoms $j$ connected to atom $i$ by a
constraint, and the $\lambda_{ij}$ are time-dependent Lagrange multipliers. Note
$\sigma_{ij} = \sigma_{ji}$ and $\lambda_{ij} = \lambda_{ji}$. Since
$\nabla_i \sigma_{ij} = 2\mathbf{r}_{ij}$ with
$\mathbf{r}_{ij} \equiv \mathbf{r}_i - \mathbf{r}_j$, the constraint force reduces to
$\mathbf{G}_i = -2\sum_j{}' \lambda_{ij}\mathbf{r}_{ij}$.

The first equation of RATTLE, Eq. (2.7), in this notation is

<!-- eq:A3 -->
$$ \mathbf{r}_{i}(t+h) = \mathbf{r}_{i}(t) + h\dot{\mathbf{r}}_{i}(t) + (h^{2}/2m_{i})\left[\mathbf{F}_{i}(t) - 2\sum_{j}{}' \lambda_{RRij}(t)\, \mathbf{r}_{ij}(t)\right], $$

where $\mathbf{r}_{ij}(t) = \mathbf{r}_i(t) - \mathbf{r}_j(t)$. The
$\lambda_{RRij}(t)$ are chosen so that the constraints (A1) are satisfied at $t+h$. The
second equation of RATTLE, Eq. (2.8), is

<!-- eq:A4 -->
$$ \dot{\mathbf{r}}_{i}(t+h) = \dot{\mathbf{r}}_{i}(t) + (h/2m_{i}) \left[ \mathbf{F}_{i}(t) - 2 \sum_{j}{}' \lambda_{RRij}(t)\, \mathbf{r}_{ij}(t) + \mathbf{F}_{i}(t+h) - 2 \sum_{j}{}' \lambda_{RVij}(t+h)\, \mathbf{r}_{ij}(t+h) \right]. $$

The $\lambda_{RVij}(t+h)$ are chosen to satisfy the time derivatives of the constraints
(A2) at $t+h$. The iterative method of Ryckaert et al. is applicable to solving for
$\lambda_{RRij}$ and $\lambda_{RVij}$ (see Appendix C).

## Appendix B: Global Error of RATTLE

### Derivation (not implemented)

We prove that the local error of RATTLE equals that of velocity Verlet for
unconstrained dynamics; the global errors are then both of order $h^2$.

Consider Eq. (A3). If $\lambda_{RRij}(t)$ were replaced by the exact (unknown)
$\lambda_{ij}(t)$, then (A3) would be a second-order Taylor expansion for
$\mathbf{r}_i(t+h)$; the result would have errors of order $h^3$ and the constraints
would be violated by amounts of order $h^3$. Therefore, if the $\lambda_{RRij}(t)$ are
chosen to satisfy the constraints exactly, they differ from the correct
$\lambda_{ij}(t)$ by amounts of order $h$,

<!-- eq:B1 -->
$$ \lambda_{RRij}(t) = \lambda_{ij}(t) + O(h), $$

since $\lambda_{RRij}$ is multiplied by a factor $h^2$ in (A3). It follows that (A3)
predicts positions at $t+h$ differing from the Taylor prediction by order $h^3$, hence
from the correct positions by order $h^3$. Thus the local error in $\mathbf{r}(t+h)$
from RATTLE is of order $h^3$, same as velocity Verlet.

Next consider Eq. (A4). Replace $2\sum_j{}' \lambda_{RRij}(t)\mathbf{r}_{ij}(t)$ by the
equivalent quantity

<!-- eq:B-aux -->
$$ 2 \sum_{j}{}' \lambda_{ij}(t) \mathbf{r}_{ij}(t) + 2 \sum_{j}{}' \left[ \lambda_{RRij}(t) - \lambda_{ij}(t) \right] \mathbf{r}_{ij}(t+h) + 2 \sum_{j}{}' \left[ \lambda_{RRij}(t) - \lambda_{ij}(t) \right] \left[ \mathbf{r}_{ij}(t) - \mathbf{r}_{ij}(t+h) \right]. $$

Then (A4) becomes

<!-- eq:B2 -->
$$ \dot{\mathbf{r}}_{i}(t+h) = \dot{\mathbf{r}}_{i}(t) + (h/2m_{i}) \left[ \mathbf{F}_{i}(t) - 2 \sum_{j}{}' \lambda_{ij}(t) \mathbf{r}_{ij}(t) + \mathbf{F}_{i}(t+h) - 2 \sum_{j}{}' (\lambda_{RVij}(t+h) + \lambda_{RRij}(t) - \lambda_{ij}(t)) \mathbf{r}_{ij}(t+h) - 2 \sum_{j}{}' [\lambda_{RRij}(t) - \lambda_{ij}(t)] [\mathbf{r}_{ij}(t) - \mathbf{r}_{ij}(t+h)] \right]. $$

The last sum contributes an amount of order $h^3$ to $\dot{\mathbf{r}}_i$ because of
(B1). Neglecting that term, if $\lambda_{RVij}(t+h)$ were chosen so that

<!-- eq:B3 -->
$$ \lambda_{RVij}(t+h) + \lambda_{RRij}(t) - \lambda_{ij}(t) = \lambda_{ij}(t+h), $$

then (B2) would be a Taylor series for $\dot{\mathbf{r}}_i$ correct to order $h^2$. The
resulting velocities would contain errors of order $h^3$ and fail to satisfy (A2) by
order $h^3$. Restoring the last sum introduces additional errors of order $h^3$. By
making changes of order $h^2$ to $\lambda_{RVij}(t+h)$, the constraints can be
satisfied exactly, giving

<!-- eq:B4 -->
$$ \lambda_{RVij}(t+h) = \lambda_{ij}(t+h) + \lambda_{ij}(t) - \lambda_{RRij}(t) + O(h^2). $$

The RATTLE velocities differ from the exact Taylor expression by order $h^3$. Hence the
local error in $\dot{\mathbf{r}}(t+h)$ is of order $h^3$, the same as velocity Verlet.

## Appendix C: Iterative Procedure for RATTLE Calculations

Suppose positions, velocities, and intermolecular forces are known at time $t$, and we
wish to calculate the corresponding quantities at $t+h$. This is a modification of the
Ryckaert et al. iterative procedure for SHAKE. Define

<!-- eq:C1 -->
$$ g_{ij} = h\lambda_{RRij}(t), $$

<!-- eq:C2 -->
$$ k_{ij} = h\lambda_{RVij}(t+h), $$

<!-- eq:C3 -->
$$ \mathbf{q}_i = \dot{\mathbf{r}}_i(t) + (h/2m_i)\, \mathbf{F}_i(t) - (1/m_i) \sum_j g_{ij}\, \mathbf{r}_{ij}(t). $$

Then Eqs. (A3) and (A4) can be expressed as

<!-- eq:C4 -->
$$ \mathbf{r}_i(t+h) = \mathbf{r}_i(t) + h\mathbf{q}_i, $$

<!-- eq:C5 -->
$$ \dot{\mathbf{r}}_i(t+h) = \mathbf{q}_i + (h/2m_i)\mathbf{F}_i(t+h) - (1/m_i)\sum_j k_{ij}\mathbf{r}_{ij}(t+h). $$

### Position iteration (solve for the $\mathbf{q}_i$)

To start, let

<!-- eq:C6 -->
$$ \mathbf{q}_i = \dot{\mathbf{r}}_i(t) + (h/2m_i)\, \mathbf{F}_i(t), \qquad i = 1,\dots, N. $$

The iterative loop begins. Pick a constraint involving atoms $i$ and $j$. Let

<!-- eq:C7 -->
$$ \mathbf{s} = \mathbf{r}_i(t) + h\mathbf{q}_i - \mathbf{r}_j(t) - h\mathbf{q}_j. $$
<!-- CHECK: OCR had "r_i(t) + h q_i(t) - r_i(t) - h q_i(t)"; s must be the i-j displacement, so the second pair is atom j. -->

$\mathbf{s}$ is the current approximation for the vector displacement of atoms $i$ and
$j$. If $|\mathbf{s}|^2 - d_{ij}^2$ differs from zero by less than an acceptable
tolerance, pick a new constraint. Otherwise correct $\mathbf{q}_i$ and $\mathbf{q}_j$.
Let

<!-- eq:C8 -->
$$ \mathbf{r}_i^{T} = \mathbf{r}_i(t) + h[\mathbf{q}_i - g\,\mathbf{r}_{ij}(t)/m_i], $$

<!-- eq:C9 -->
$$ \mathbf{r}_j^{T} = \mathbf{r}_j(t) + h[\mathbf{q}_j + g\,\mathbf{r}_{ij}(t)/m_j]. $$

These are the new values for $\mathbf{r}_i(t+h)$ and $\mathbf{r}_j(t+h)$ when the
corrections proportional to $g$ are made to $\mathbf{q}_i$ and $\mathbf{q}_j$. Choose
$g$ so that

<!-- eq:C10 -->
$$ |\mathbf{r}_i^{T} - \mathbf{r}_j^{T}|^2 = d_{ij}^2. $$

Solving for $g$ (neglecting quantities of order $g^2$):

<!-- eq:C11 -->
$$ g = \frac{s^2 - d_{ij}^2}{2h\,[\mathbf{s} \cdot \mathbf{r}_{ij}(t)]\,(m_i^{-1} + m_j^{-1})}. $$

Then replace $\mathbf{q}_i$ by $\mathbf{q}_i - g\,\mathbf{r}_{ij}(t)/m_i$ and
$\mathbf{q}_j$ by $\mathbf{q}_j + g\,\mathbf{r}_{ij}(t)/m_j$, and pick a new constraint.
Continue until all constraints are satisfied to within tolerance.

Once $\mathbf{q}_i$ and $\mathbf{r}_i(t+h)$ are known for all $i$, the forces at $t+h$
can be calculated. Before doing so, place the positions at $t+h$ in the memory that
held positions at $t$, and place the $\mathbf{q}_i$ in the memory that held velocities
at $t$, so the algorithm uses just $3N$ memory locations.

### Velocity iteration (solve for the $\dot{\mathbf{r}}_i(t+h)$)

To start, let

<!-- eq:C12 -->
$$ \dot{\mathbf{r}}_i(t+h) = \mathbf{q}_i + h\mathbf{F}_i(t+h)/2m_i, \qquad i = 1,\dots,N. $$

The iterative loop begins. Pick a constraint involving atoms $i$ and $j$. Calculate the
dot product $\mathbf{r}_{ij}(t+h)\cdot[\dot{\mathbf{r}}_i(t+h)-\dot{\mathbf{r}}_j(t+h)]$.
If it differs from zero by less than an acceptable tolerance, pick another constraint.
Otherwise correct the two velocities. Let

<!-- eq:C13 -->
$$ \dot{\mathbf{r}}_i^{T} = \dot{\mathbf{r}}_i(t+h) - k\,\mathbf{r}_{ij}(t+h)/m_i, $$

<!-- eq:C14 -->
$$ \dot{\mathbf{r}}_j^{T} = \dot{\mathbf{r}}_j(t+h) + k\,\mathbf{r}_{ij}(t+h)/m_j, $$

with $k$ chosen so that $\dot{\mathbf{r}}_i^{T} - \dot{\mathbf{r}}_j^{T}$ is
perpendicular to $\mathbf{r}_{ij}(t+h)$:

<!-- eq:C15 -->
$$ k = \frac{\mathbf{r}_{ij}(t+h) \cdot [\dot{\mathbf{r}}_i(t+h) - \dot{\mathbf{r}}_j(t+h)]}{d_{ij}^2\,(m_i^{-1} + m_j^{-1})}. $$

Then replace $\dot{\mathbf{r}}_i(t+h)$ by $\dot{\mathbf{r}}_i^{T}$ and
$\dot{\mathbf{r}}_j(t+h)$ by $\dot{\mathbf{r}}_j^{T}$, and pick another constraint.
Continue until all velocity constraints are satisfied to within tolerance.
