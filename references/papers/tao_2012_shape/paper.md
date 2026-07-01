# Maintain rigid structures in Verlet based Cartesian molecular dynamics simulations

Peng Tao, Xiongwu Wu, Bernard R. Brooks. J. Chem. Phys. 137, 134110 (2012). DOI: 10.1063/1.4756796.

## Abstract

An algorithm is presented to maintain rigid structures in Verlet based Cartesian
molecular dynamics (MD) simulations. After each unconstrained MD step, the
coordinates of selected particles are corrected to maintain rigid structures
through an iterative procedure of rotation matrix computation. This algorithm,
named SHAPE and implemented in the CHARMM program suite, avoids the calculation
of Lagrange multipliers, so that the complexity of computation does not increase
with the number of particles in a rigid structure. The implementation does not
require significant modification of the propagation integrator, and can be
plugged into any Cartesian based MD integration scheme. A unique feature of the
SHAPE method is that it is interchangeable with SHAKE for any object that can be
constrained as a rigid structure using multiple SHAKE constraints. Unlike SHAKE,
SHAPE can be applied to large linear (three or more centers) and planar (four or
more centers) rigid bodies. Numerical tests with four model systems including two
proteins demonstrate that the accuracy and reliability of SHAPE are comparable to
SHAKE, but with much more applicability and efficiency.

## I. Introduction

Rigid body molecular dynamics is of increasing importance in computational
chemistry and physics. Generally, there are three categories of rigid body
dynamics methods:

1. Non-constrained dynamics integration is applied in each step before applying
   constrained corrections. SHAKE and RATTLE belong to this category.
2. Rigid body constraint forces are calculated to propagate each particle in the
   rigid bodies.
3. The rigid body motion is divided into translation of the center of mass and
   rotation about the center of mass. Translation and rotation are propagated
   based on the force on the center of mass and the torque about the center of
   mass.

Each of these has undesirable features, limitations, or restrictions. For a
nonlinear rigid structure containing M particles, `3M - 6` Lagrange equations
are typically solved to implement rigid structure constraints; this is
inconvenient for large M. This work presents SHAPE, an efficient algorithm for
simulation of rigid structures composed of an arbitrary number of particles.
SHAPE is fully consistent and interchangeable with the SHAKE method, and can be
combined with SHAKE so that rigid body constraints are solved consistently with
SHAKE for atoms involved in more than one type of constraint. For constant
pressure simulation, this method can use the same virial correction scheme as in
SHAKE.

## II. Theory

Assume that M atoms compose a rigid structure in an N-atom system. During each MD
step, the coordinates of the whole system with N atoms are integrated according
to the equation of motion. The coordinates of the M atoms are then corrected to
maintain the rigid structure.

For an MD integration at time `t` with `Δt` as step size and without any
constraint, the coordinates of the M atoms change from `r_i(t)` to
`r_i^non(t+Δt)`. Either with or without the rigid structure constraint, the
structure's momentum and angular momentum from `t` to `t+Δt` should be the same,
which lays the theoretical ground for maintaining rigid structures.

The centers of mass (COM) for two sets of coordinates are given by eq. (1)–(2).
The rigid structure motion of M atoms can be separated into translational and
rotational parts. The translational motion is represented by the COM of these M
atoms. The rotational motion is treated in a body-fixed local coordinate with its
origin at the COM of the M atoms and XYZ orientation the same as the global
coordinate. The body-fixed local coordinate is constructed by transforming the
origin of the global Cartesian coordinate to the COM without any rotation. The
coordinates of the M atoms at time `t` and `t+Δt` in the body-fixed local
coordinate are given by eq. (3) and (4).

The angular momentum `L^non` of the M atoms at time `t+Δt/2` without rigid body
constraints is given by eq. (5), where `p_i^b` is the momentum and `v_i^b` is the
velocity for atom `i` in the body-fixed local coordinate at time `t+Δt/2`. The
velocities in body-fixed local coordinate are finite differences of coordinates
calculated from eq. (3) and (4).

The rigid body motion of M atoms is achieved by correcting the coordinates of the
M atoms at time `t+Δt` from `r_i^non,b(t+Δt)`, to maintain the rigid structure,
and to have the angular momentum equal to `L^non(t+Δt/2)`. The correction is
computed in the body-fixed local coordinate through the iterative process below.

Assume the moment of inertia `I`, a 3×3 matrix, for the M atoms in body-fixed
local coordinates `r_i^b(t)` is given by eq. (6). The angular momentum of the
rigid structure can be calculated from `I` and the angular velocity vector
`ω^rig(t+Δt/2)` by eq. (7). Based on the equality of angular momentum (eq. 8),
plugging eq. (5) and (7) into eq. (8) gives eq. (9); rearranging with the inverse
of `I` gives the first estimate `ω^rig,(1)` (eq. 10).

Once vector `ω` is obtained, a skew-symmetric matrix `ω̂` is constructed (eq. 11).
The angle `θ` undertaken by the rigid body in time interval `Δt` is `‖ω‖Δt`, and
the corresponding skew-symmetric matrix `θ̂` is `ω̂Δt`. The rotation matrix
corresponding to angle `θ` can be expressed as the matrix exponential
`R = exp(θ̂)` (eq. 12). Following Rodrigues' formula, the matrix exponential is
given by eq. (13). Equation (13) faces numerical instability when `‖θ‖` is small;
in that situation a Taylor expansion is applied (eq. 14).

After obtaining rotation matrix `R`, the body-fixed local coordinates of the M
atoms after the current MD step are updated through rotational operation on the
coordinates in body-fixed local coordinate at time `t` (eq. 15), where superscript
`(n)` specifies the iteration number (n=1 for the first iteration). The rotational
operation `R^(n)` is applied to the M atoms around the center of mass.

From the updated coordinates, the angular momentum in the iteration, `L^rig,(n)`,
is calculated (eq. 16), where `ṙ_i^rig,b,(n)` is the velocity of atom `i` from the
updated coordinate in body-fixed local coordinate during the iteration.

If `L^rig,(n)(t+Δt/2)` differs from `L^non(t+Δt/2)` by less than a predefined
tolerance, the procedure is considered converged. Otherwise, the angular velocity
vector `ω^rig,(n)(t+Δt/2)` is updated for the next iteration.

Since the coordinates at time `t+Δt/2` are not readily available, they are
approximated using the square root of rotation matrix `R^(n)` (eq. 17).
Consequently, the angular momentum for both non-rigid and rigid structures can be
computed as in eq. (18) and (19). The updated angular velocity vector
`ω^rig,(n+1)(t+Δt/2)` is expected to satisfy eq. (20). Combining eq. (20) and (7)
gives eq. (21). Multiplying both sides by `I^-1` and plugging eq. (18) and (19)
into eq. (21) with rearrangement gives eq. (22).

Since `R^(n)` is close to the identity matrix for small time step `Δt`, we write
`R^(n) = 1 + δR` (eq. 23) with `δR` small. Hence the square root can be
approximated (eq. 24), and by neglecting the quadratic term, the square root of
`R^(n)` is approximated by eq. (25). Plugging eq. (25) into eq. (22) gives the
angular velocity updating scheme, eq. (26).

Using the updated angular velocity vector `ω^rig,(n+1)`, the procedure starting
from eq. (10) is repeated until the angular momentum `L^rig,(n)(t+Δt/2)`
converges to `L^non(t+Δt/2)`. Typically, three iterations (n=3) are sufficient to
converge to double precision accuracy. This iterative approach is far more
computationally efficient than any analytic approach.

After convergence, the coordinates in the global reference for the rigid body M
atoms are calculated by eq. (27). The final coordinates `r_i^rig,(n)(t+Δt)` are
the rigid structure coordinates for the desired MD trajectory.

The procedure iteratively solves eq. (8), or more explicitly eq. (28), where
`r_i^non,b(t+Δt)` and `r_i^rig,b(t+Δt)` are the body-fixed local coordinates
without and with rigid body constraints at time `t+Δt`, respectively.

### Sharing atoms between rigid bodies (outer loop)

One important aspect of this rigid body integration method is that different rigid
bodies may share one common atomic center. A rigid body may also share atoms
involved in SHAKE constraints. This is possible because, as with SHAKE, applying
the constraint to partial forces in multiple steps yields the same solution as
applying the full force once, to within numerical precision. When a given atom is
involved in two or more SHAKE constraints, the constraints are iterated in a
cyclic fashion until convergence. The same approach works for SHAPE: if an atom is
involved in multiple SHAPE constraints, these must also be applied cyclically
until convergence. As with SHAKE, each application of the constraint conserves
angular and linear momentum holonomically. An NVE simulation converged with this
approach also approximately conserves energy. An efficient implementation
converges both SHAKE and SHAPE constraints simultaneously within the same outer
iteration loop. A prospective use is making peptide bonds rigid and planar in a
protein: each Cα is involved in two constraints that must be solved to
consistency.

When using this method on a large-scale parallel machine where different atom
centers are integrated on different processors, large rigid bodies may span more
than one processor. It is necessary to communicate the components of `r_COM^non`
and `L^non` at the beginning of the cycle and distribute the `R` matrix at the
end. It is never necessary to distribute coordinates, since each processor can
apply the final `R` matrix to atoms in each rigid body it contains. Thus the
method can be efficiently applied even if whole proteins are made rigid.

The SHAPE rigid structure algorithm is implemented in the CHARMM program suite
within the SHAPE module.

### Derivation (not implemented): order of accuracy

The SHAPE method is based on the conservation of linear and angular momentum. The
body-fixed local coordinate system guarantees conservation of linear momentum
after rigid body correction of each MD step. Therefore the major error comes from
the angular momentum calculation. Using the leapfrog integration scheme, once
convergence is reached, we have eq. (29). The accuracy of integration is
determined by the angular momentum.

Expanding the angular momentum (eq. 30) shows the leapfrog half-step angular
momentum is reproduced to `O(Δt^2)`, where `τ` is the torque on the rigid
structure. The angular velocity at time `t+Δt/2` (eq. 31) and the rotational
angle (eq. 32) are computed to `O(Δt^2)` and `O(Δt^3)` respectively, and the
coordinates at time `t+Δt` (eq. 33) are accurate to `O(Δt^3)`. Therefore, for a
single time step this algorithm has third order of accuracy. Over a given time
period `t`, the number of steps is of order `Δt^-1`, so the global error is of
order `O(Δt^2)`.

## III. Numerical Tests

Four model systems were used as test cases:

- **System a**: a water box with 126 TIP3P water molecules. Only one water
  molecule was treated as a rigid structure using SHAPE. SHAKE for holonomic
  constraint was also applied to the same water molecule for comparison, and an
  unconstrained simulation was run as benchmark.
- **System b**: Trp-cage protein (PDB 1L2Y) in a box of 1169 water molecules. The
  side chain of residue Trp6 and the whole residue Pro20 were treated as separate
  rigid structures.
- **System c**: matrix metalloproteinase 2 (MMP2) with its substrate in a box of
  12 636 water molecules. The active site (side chains of His288, His292, His298,
  Glu289, the thiirane ring, methylene and sulfone groups in substrate, and the
  zinc ion, 43 atoms from 6 residues) is combined as a single rigid structure.
- **System d**: a nine-residue β-hairpin peptide in a box of 290 water molecules,
  to test the outer loop of SHAPE handling one atom shared by two rigid
  structures. All eight peptide planes are treated as individual rigid structures,
  so seven Cα carbons are shared by adjacent peptide bonds. SHAKE cannot impose
  planar constraints; in the comparison SHAKE constrains all main-chain chemical
  bonds without treating peptide planes as rigid structures.

A time step of 1 fs was used for systems a, b, c (also repeated with 1.5 and
2.0 fs for time-step analysis). System d was run for 100 ps with 1 fs time step.
The total energies reported include a high frequency correction based on the
expectation value of the symplectic shadow Hamiltonian.

For system a, the NVE simulation using SHAPE gave results essentially
interchangeable with SHAKE in terms of total energy and standard deviation along
the 1 ns trajectory. Both approximately conserve total energy at a level very
close to the unconstrained simulation. The average and standard deviation of RMSD
for the rigid water molecule over the 1 ns trajectory using SHAPE showed even
higher accuracy than SHAKE.

For systems b and c, the dynamics using SHAPE showed the same level of
conservation of total energy as the free MD simulation without constraints, and
maintained high accuracy of the rigid structures. For systems b and c it is not
practical to implement rigid structure constraints using SHAKE.

Time-step analysis: for system a, SHAPE showed very similar accuracy to SHAKE at
1.5 fs and only slightly worse at 2.0 fs. For systems b and c, SHAPE showed
comparable-to-slightly-better accuracy at 1.5 fs and marginally worse at 2.0 fs.
The time step had very little effect on maintaining the rigid structure. SHAPE is
robust across time steps.

For system d, the fluctuation of total energy (with high frequency correction)
from SHAPE is slightly higher but comparable to SHAKE. Both display observable
overall energy drift, with SHAKE slightly better. On average, the SHAPE outer loop
takes five more iterations to converge than SHAKE. In an additional SHAPE
simulation where each Cα belongs to the previous peptide plane only (no atoms
shared), the total energy fluctuation and energy drift are smaller than both
SHAPE (sharing) and SHAKE, and only one outer-loop iteration is needed.

## IV. Concluding Remarks

SHAPE implements rigid body constraints in MD through an iterative procedure of
constructing a rotation matrix for the rigid structure to preserve angular
momentum. It avoids computing Lagrange multipliers of individual holonomic
constraints. An arbitrary number of particles can form a single rigid structure,
and it is more efficient than SHAKE for systems with more than three atoms. An
arbitrary number of such rigid structures can be implemented. The iterative
procedure is independent from and conducted after each MD free integration step,
so no significant modification of an existing MD integrator is necessary. It can
be plugged into any simulation package that already includes SHAKE. Multiple
rigid structures can be defined, even those sharing one particle with neighboring
rigid structures.

Since this method is interchangeable with SHAKE, it suffers from the same failures
SHAKE exhibits at very large time steps or very high temperatures where the
equations cannot be solved. In that case, non-energy-conserving approximations are
used while temperatures are too large, as an alternative to program termination.

SHAPE conserves linear and angular momentum and is the only rigid body integration
algorithm fully consistent and interchangeable with SHAKE. It shares with SHAKE
both its holonomic and symplectic nature, and well conserves total energy for long
simulations with minimal energy drift. It extends SHAKE-like constraints to linear
systems with three or more atoms, planar systems with four or more atoms, and to
larger rigid structures where SHAKE is intractable.
