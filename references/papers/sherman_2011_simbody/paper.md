# Simbody: multibody dynamics for biomedical research

Michael A. Sherman, Ajay Seth, Scott L. Delp (Stanford University).
Procedia IUTAM 2 (2011) 241-261, 2011 Symposium on Human Body Dynamics.

## TOC

- 1. Introduction
- 2. Simbody overview (scope, top-level architecture System/State/Study, state handling)
- 3. Time stepping (DAE/DEM formulation, integrators, accuracy control, real-time)
- 4. Formulation of the multibody system (equations of motion, solving for accelerations)
- 5. Contact modeling (Hertz/Hunt-Crossley/Stribeck, Elastic Foundation, rigid contact)
- 6. Future directions
- Appendix A. Contact model details

## Abstract

Simbody is an open-source, extensible, high-performance multibody-mechanics library
aimed at biomedical research. It supports neuromuscular, prosthetic, and biomolecular
simulation, plus biologically-inspired robot/avatar control. Simbody is the dynamics
engine behind OpenSim. This article reviews the architecture, theory, and computational
methods Simbody uses. Simbody provides a minimal-coordinate O(n) multibody formulation
with robust constraint handling and compliant contact.

## 1. Introduction

Multibody dynamics tools developed for mechanical/aerospace engineering are hard to
apply directly to biological systems: biomechanical joints do not perform simple fixed-axis
rotations and may have several moving parts; soft-tissue contact involves large deformation;
joint actuation is redundant; parameters are hard to measure and measurements are noisy.
Concepts like "generalized coordinate" or "moment arm" become hard to define precisely.
Game physics engines (e.g. ODE) attain real-time performance via simplified theory that may
not converge to correct results ("ODE should not be used for quantitative engineering").
Methods built on sound theory with selectable accuracy can be fast yet converge to high
fidelity. Simbody is developed under the NIH-funded Simbios center as part of the SimTK
biosimulation toolkit.

## 2. Simbody overview

Simbody is a C++ API providing robust, high-performance, minimal-coordinate O(n) multibody
dynamics. Applications range from biomolecular machines (amino/nucleic-acid components) to
pathological-gait musculoskeletal models to biologically inspired robots. Accuracy is
user-controllable from interactive real-time to high-fidelity.

### 2.1 Scope

Beyond the equations of motion, Simbody includes contact modeling, numerical integration
and differentiation, constraint stabilization and redundancy handling, assembly analysis,
optimization/root finding, linear algebra, event isolation/handling, accuracy control,
threading, visualization, and real-time interaction. User-written extensions are supported
for force/constraint elements and for novel internal-coordinate joints whose kinematics can
be based on empirical measurements. User-specifiable accuracy requires variable-step
integration with error estimation, which in turn requires strict management of system state
to avoid referencing out-of-date computations after a step rejection.

### 2.2 Top-level architecture (System / State / Study)

The three primary objects are the **System**, **State**, and **Study**.

- A **System** encapsulates the model components (bodies, joints, force elements) and the
  computation code. It defines the parameterization but is itself **stateless** and unchanged
  during a study.
- A **State** holds a complete set of values for the System's parameters. The response of a
  System is completely determined by the state values presented to it. "State" here means
  *everything* variable about a System: continuous time/position/velocity, discrete variables,
  event memory, modeling choices, and *instance variables* (e.g. masses, lengths).
- A **Study** couples a System with one or more States and represents a computational
  experiment. Any Study's result can be expressed as a state (or trajectory - a series of
  states) plus quantities the System can compute from those values.

Examples of Studies: evaluation (query a quantity), assembly (satisfy position constraints
like loop closures), inverse kinematics (fit marker observations), dynamics (satisfy Newton's
laws), energy minimization, Monte Carlo (trajectory satisfying a probability distribution such
as Boltzmann), design/parameter-fitting, and modeling (select among algorithmic choices).
Separate State objects allow trajectories to be saved/restored by copying State objects, with
no hidden state missed, and trial states to be generated/discarded cheaply.

### 2.3 Handling of state (realization cache, stages)

State-dependent computations (e.g. reaction forces) can be expensive; Simbody caches them in
a **realization cache** held inside the State. *Realizing* a State means presenting it to a
System to compute the physical consequences of its values. Evaluation proceeds in ordered
**stages**: positions before velocities, velocities before forces, applied forces before
accelerations (topology is the first/construction stage). A change to a state variable at
stage s invalidates all cache entries at stages s and above; access to an invalidated value
triggers recomputation or an error. Cache values are not logically part of the state (they
only affect speed, not behavior). This is an architectural replacement for error-prone ad hoc
"isValid" flags.

## 3. Time stepping

Simbody simulations are **hybrid**: continuous time evolution plus discrete events. A
**time stepper** consists of a numerical integrator (advancing smooth intervals and detecting
pending events) and an event handler.

### 3.1 Formulation as seen by the time stepper (DAE / DEM)

The equations seen by the time stepper are the ODE (Eq. 1), algebraic constraints (Eq. 2),
and event-detection functions (Eq. 3). Together (1) and (2) form a differential-algebraic
equation (DAE). More precisely the continuous system is a **differential equation on a
manifold (DEM)**: when constraint equations (2) hold, their time derivatives hold
automatically, i.e. the DEM condition (Eq. 4) `c = 0 => cdot = 0`. Simbody produces
differential equations that guarantee this, allowing conventional ODE integrators (augmented
with coordinate projection) rather than a general DAE solver.

### 3.2 Advancing the continuous system

Integrators offered: error-controlled variable-step explicit Runge-Kutta (biomechanics /
real-time), a velocity Verlet method (biomolecular systems), and a variable-order backward
difference formula (BDF) implicit integrator **CPODES** for stiff systems (a modification of
CVODE adding coordinate projection). All integrators offer continuous (dense) output so step
size is decoupled from reporting interval, and monitor event functions (Eq. 3) for sign
transitions to isolate event times.

Because of the DEM condition (Eq. 4), a trajectory started on the manifold stays on it under
exact integration; truncation error causes drift. Simbody removes drift via **coordinate
projection** using Eq. 2, which is superior to Baumgarte stabilization (no feedback gains to
tune, guarantees the solution lies on the manifold each step, and - if the projection is
normal to the manifold in a suitable norm - reduces the error estimate, permitting larger
steps). Coordinate projection is cheap when there are few constraints (small projection
matrices); Simbody's internal-coordinate joints eliminate most constraints.

### 3.3 Controlling accuracy

Simbody exposes a single scalar accuracy parameter alpha (≈ "% relative error" / number of
significant digits), related to digits by Eq. 5, `alpha = 10^-n`. alpha resembles the relative
tolerance (rtol) of standard integrators, but there is no absolute-tolerance (atol) equivalent.
Instead Simbody combines alpha with internally-computed scaling to form two step-acceptance
tests (Eqs. 6, 7) on weighted state error and weighted constraint violation, using diagonal
weighting matrices W and T that map heterogeneous errors (lengths, angles, velocities) to unit
errors. This gives non-expert users a single predictable "knob."

### 3.4 Real-time interaction

Real-time interaction requires matching simulated time to clock time and steady frame delivery.
Fixed-step integration wastes work because step size must handle the worst-case interval.
Simbody instead uses variable-step integration with a "local lag" delay buffer (adapted from
networked games): frames are collected regularly in simulated time but arrive at variable real
rate, and extracted at a fixed real rate delayed slightly (up to t_delay, e.g. 100-120 ms,
imperceptible to humans). Interpolation (dense output) lets steps up to t_delay be taken while
delivering regularly-spaced frames; each step has the same cost, interpolated frames are nearly
free. Separation of System from State makes this easy (copy whole State objects into the buffer).

## 4. Formulation of the multibody system

Simbody uses a generalized-coordinate formulation with the fewest possible coordinates. It does
not reduce the system to an ODE; instead it chooses coordinates as a basis and restricts motion
to a constraint manifold, giving a compact representation with robust handling of constraint
stabilization, poor conditioning, and redundant constraints. The unconstrained solver is a
recursive O(n) spatial-operator method (Rodriguez and Jain), extending Schwieters's templatized
C++ implementation. The m constraint equations are adjoined, introducing an O(m^3) term; but
biologically-realistic internal-coordinate joints eliminate most constraints so n >> m, and only
coupled constraint blocks need be solved simultaneously.

### 4.1 Equations of motion

A Simbody system is a tree of **mobilized bodies**, each a body plus its unique inboard
internal-coordinate joint (**mobilizer**). "Mobilizer" avoids confusion with the biological term
"joint" (which might be a mobilizer, a constraint, force elements, or a combination). The i-th
mobilizer gives its body 0 <= nu_i <= 6 degrees of freedom relative to its parent, parameterized
by nu_i generalized speeds u_i (the basis for the equations of motion) and nq_i >= nu_i
generalized coordinates q_i (the pose). Coordinate rates relate to speeds by the kinematic
differential equation (Eq. 8, per-mobilizer; Eq. 9, whole system with block-diagonal N). When
nq_i > nu_i there are local constraints (almost always quaternion normalization) - introduced
only for numerical stability, producing no forces.

Constraints may be specified geometrically (point distance, non-penetration, non-slip) or
directly on speeds/coordinates (prescribed motion, couplers), plus linear acceleration-only
constraints. Simbody reduces these to holonomic (Eq. 10), nonholonomic (Eq. 11), and
acceleration-only (Eq. 12) algebraic relationships. Time derivatives add the holonomic velocity
and acceleration constraints (Eqs. 13, 14) and nonholonomic derivative (Eq. 15), with
`P = (∂p/∂q) N`, `V = ∂v/∂u`, and `b_p`, `b_v` collecting terms independent of udot. The
acceleration constraints are stacked (Eq. 16) with `G = [P V A]^T`, `b = [b_p b_v b_a]^T`.
G is the acceleration constraint Jacobian, generally poorly conditioned or singular due to
redundant constraints. The dynamic equations are the constrained equations of motion (Eq. 17)
and auxiliary ODEs (Eq. 18). Here `y = {q, u, z}` and c comprises `p`, `pdot`, `v`.

### 4.2 Solving for accelerations

For forward dynamics, Eqs. 16 and 17 are solved simultaneously for udot and lambda. When G has
full row rank (rank(G) = m) the solution is unique. When rank(G) < m (redundant constraints),
lambda is underdetermined although udot is still unique. Rather than dropping constraints (which
can give absurd results, as in SD/FAST), Simbody computes a **least-squares** solution for
underdetermined lambda, which often returns the infinite-stiffness limit of compliant elements
and always spreads the load so no redundant constraint carries zero.

Simbody does not normally form the n x n mass matrix M (that would be O(n^2), inverting O(n^3)).
The operator `M^{-1} v` for any vector v is available in O(n) via recursive spatial operators.
Eliminating udot from Eqs. 16, 17 gives Eq. 19 with `udot_0 = M^{-1}(f_applied - f_inertial)`
the unconstrained acceleration (O(n)); the RHS `g_0` (the unconstrained acceleration constraint
errors) is O(n+m). Defining `Y = G M^{-1} G^T` (Eq. 20), each column costs O(n+m) so forming Y
is O(mn + m^2). If Y has full rank, solve by LU with pivoting; to drop equations in the singular
case, use QR; for a least-squares solution use the pseudoinverse (Eq. 21) computed via complete
orthogonal factorization (QTZ, ~5x faster than SVD, via LAPACK). All factorizations are O(m^3),
so total cost is O(m^3 + mn + m^2). With multipliers known, `f_constraint = G^T lambda` (O(n))
and Eq. 17 gives the final constrained acceleration (a final O(n) M^{-1} application).

## 5. Contact modeling

Real contact forces arise from compliant-material deformation. Simbody provides two compliant
contact models: **Hertz** (analytic, accurate, simple geometry only) and the **Elastic
Foundation Model (EFM)** (meshes for arbitrary surfaces, simplified elastic model). Both use a
Hunt-Crossley dissipation model and Stribeck friction. Each contact element produces a force
with three effects (Eq. 22): stiffness, dissipation, and friction. f_stiffness differs between
the two models; dissipation and friction use the same method once stiffness is known.

### 5.1 Hertz stiffness, Hunt-Crossley dissipation, Stribeck friction

Hertz theory requires two linearly elastic materials in non-conforming contact, contact-patch
dimensions small vs curvatures and object size, contact initiating at a common point with
opposed surface normals, each surface well-approximated by a paraboloid there. The relative
separation is a paraboloid described by two principal curvatures. Simbody provides Hertz contact
for planes, spheres, and ellipsoids (and line/cylindrical contact).

Deformation is a scalar `x` (total deformation of the two surfaces along the contact normal,
x > 0 when contacting). The Hertz normal force magnitude is Eq. 23; R is a composite relative
radius of curvature, E* a composite elastic modulus, sigma an eccentricity factor (sigma = 1 for
circular contact, growing via elliptic integrals otherwise). The bracketed quantity is
independent of x, so the force-displacement law is nonlinear (x^{3/2}) purely from changing
geometry. Contact point P lies along the line between the initiation points; if materials are
equal, P is midway; if one is much stiffer, P is near the stiff (non-deforming) surface. Force
is applied at P along the normal in opposite directions, always pushing.

Hunt-Crossley dissipation (Eq. 24) uses an effective dissipation coefficient c*. Material
property c is measured as the negated slope of the coefficient-of-restitution-vs-impact-velocity
curve at low velocity: `e = 1 - c v` (e = measured restitution, v = impact speed; restitution is
not a material property, c is). f_HC is signed, creating empirical hysteresis. Since total force
`f_Hz + f_HC >= 0` under typical conditions (negative totals come from unmodeled losses like
"ringing"), Simbody clamps dissipation (Eq. 25) so the total normal force (Eq. 26) never goes
negative.

Friction: find the body stations coincident with contact point P, compute their relative
velocity v in the contact plane; with slip rate `v = |v|`, the friction magnitude is Eq. 27.
mu(v) is an effective friction coefficient depending only on slip velocity, parameterized by
static/dynamic/viscous friction coefficients and a transition speed at which static friction
peaks. Simbody uses a three-segment spline (first two segments quintic polynomials) for
C2-smooth transitions between static friction, the Stribeck transition, and sliding. The
transition velocity is set negligibly small; intermittent no-slip constraints can enforce exact
stiction, but the continuous method is robust without explicit stiction-event handling. A
drawback: very small transition velocity can make the system stiff, reducing explicit-integrator
efficiency; the implicit CPODES integrator handles stiffness for not-too-large systems.

### 5.2 Elastic Foundation Model

EFM treats contacting solids as rigid except for a thin elastic layer of thickness h at the
surfaces. A composite stiffness modulus E* is formed (with a different, linear combining rule
than Hertz). Each surface is approximated by a triangular mesh; at each triangle centroid a
spring of stiffness k (from triangle area, E*, and h) forms a "bed of springs." At run time,
when EFM body A contacts B, Simbody finds every triangle of A whose centroid is inside B; for
each, the closest point S on B's surface gives displacement x (centroid-to-S distance), a force
`k x` is applied at a contact point P along the centroid-S segment (P placed by relative
stiffness as in Hertz), a Hunt-Crossley-like dissipation `k x c* xdot` is added, and friction
(Eq. 27) is computed from relative velocity at P in the triangle plane. Repeated per overlapping
triangle (both directions if B is also EFM); contributions summed to net force/moment and center
of pressure. EFM is a discretization of a Winkler foundation; it does NOT account for
inter-element coupling and so does not converge to linear-elasticity results even at fine mesh
resolution, but can give good total-force agreement with FEM at much lower cost.

### 5.3 Rigid contact

An alternative treats contact objects as rigid using unilateral constraints to prevent
interpenetration and sliding, with impulsive collisions and a supplied coefficient of restitution
for dissipation. These are non-physical assumptions but useful in practice (e.g. muscle-induced
acceleration analysis, where constraint reaction forces are instantaneous). Simbody's
constraints/operators/event handling can implement this today; automated support is planned.

## 6. Future directions

Planned: rigid contacts + impulsive collisions with automated unilateral-constraint handling;
more analytical surfaces for Hertz contact; Jain and Rodriguez's speedup for prescribed motion;
Jain's method for eliminating locally-coupled constraints; a growing library of predefined
forces, constraints, and mobilizers.

## Appendix A. Contact model details

### A.1 Material-property combining rules; contact-point location

Given Young's modulus E_i and Poisson ratio nu_i (i = 1,2), form the plane-strain modulus
`E_i* = E_i / (1 - nu_i^2)`. Relative curvature is geometric: `R = R1 R2 / (R1 + R2)`.

### Derivation (not implemented): composite Hertz modulus E*

The literature suggests `E* = E1* E2* / (E1* + E2*)`, but that is inconsistent with the nonlinear
Hertz relation. Viewing body B1 meeting an infinitely rigid halfspace with its radius changed to
R gives `f1 = (4/3) sqrt(R) E1* x1^{3/2}`; symmetrically `f2 = (4/3) sqrt(R) E2* x2^{3/2}`; both
must equal `f = (4/3) sqrt(R) E* x^{3/2}` with `x = x1 + x2`. Hence
`E1* x1^{3/2} = E2* x2^{3/2} = E* (x1+x2)^{3/2}`, which raised to the 2/3 power gives
`E1*^{2/3} x1 = E2*^{2/3} x2 = E*^{2/3}(x1+x2)`, yielding the Hertz combining rule
`E* = ( E1*^{2/3} E2*^{2/3} / (E1*^{2/3} + E2*^{2/3}) )^{3/2}` and the deformation split
fractions s1, s2 (see equations.md, `eq:combining`, `eq:x-split`, `eq:cstar`). The time
derivatives xdot_1, xdot_2 split in the same ratios, giving `c* = c1 s1 + c2 s2`.

For the EFM (linear elements) the standard rule `E* = E1* E2* / (E1* + E2*)` is correct instead.

### A.2 Elliptical contact

The Hertz force (Eq. 23) uses eccentricity correction sigma (Eq. 28) with k = a/b >= 1,
m = 1 - (1/k)^2, and complete elliptic integrals K(m), E(m). The ellipse axis ratio k depends on
the principal semi-curvatures A, B (B >= A) of the separation paraboloid via Eq. 29, which is
solved numerically for k (Simbody uses approximations from ref. [70] giving sigma accurate to 5
decimal places; ref. [71] gives a machine-precision method used only for testing).
