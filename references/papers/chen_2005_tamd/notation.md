# Notation - Chen, Im, Brooks 2005 (TAMD)

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $\theta$ | internal coordinates (torsion angles) | rad, $n\times1$ | generalized coords; $n\ll 3N-6$ |
| $\dot\theta,\ddot\theta$ | hinge angular velocity / acceleration | rad/s, rad/s$^2$, $n\times1$ | |
| $n$ | number of torsional DOF | int | typically ~1/10 of $3N-6$ |
| $N$ | total atom number | int | |
| $\Omega(\theta)$ | internal-coordinate mass matrix | $n\times n$ | non-diagonal, config-dependent; SPD |
| $C(\theta,\dot\theta)$ | Coriolis/centrifugal/gyroscopic + Cartesian-force term | $n\times1$ | |
| $T(\theta)$ | generalized (applied) force | $n\times1$ | =0 in eq:5 when all FF forces in $F^{(c)}$ |
| $V_k=\mathrm{Col}[\omega_k,v_k]$ | spatial velocity of cluster $k$ | $6\times1$ | angular $\omega$ stacked above linear $v$ |
| $\alpha_k$ | spatial acceleration of cluster $k$ | $6\times1$ | $=\dot V_k$ |
| $F_k$ | hinge spatial force between clusters $k{+}1,k$ | $6\times1$ | |
| $F_k^{(c)}$ | effective Cartesian spatial force on cluster $k$ | $6\times1$ | Col[torque about origin, net force] |
| $M_k$ | spatial inertia about cluster origin | $6\times6$ | SPD |
| $T_k$ | generalized hinge force along allowed DOF | scalar (single-DOF hinge) | |
| $\phi_{x,y}$ | spatial transform between frames with origins x,y | $6\times6$ | $\phi_{k+1,k}$ used in recursions |
| $H_k$ | hinge map | $1\times6$ (single DOF) | $H_k^{T}=\mathrm{Col}[\hat h_k,0,0,0]$ |
| $\hat h_k$ | unit vector along hinge $k$ axis | $3\times1$, dimensionless | |
| $a_k$ | spatial gyroscopic/Coriolis acceleration term | $6\times1$ | fn of velocities, inertia, hinge vel |
| $b_k$ | Coriolis force term | $6\times1$ | fn of velocities, inertia, hinge vel |
| $\Phi$ | spatial operator (stacked $\phi$) | $6n\times6n$ | lower-triangular |
| $H$ | stacked hinge operator | $n\times6n$ | block diagonal $\{H_1,\dots,H_n\}$ |
| $M$ | stacked spatial inertia | $6n\times6n$ | block diagonal |
| $K$ | innovations/Kalman-gain spatial operator | operator | from articulated-body recursion |
| $D$ | articulated-body factor | $n\times n$ | diagonal + SPD for single-DOF hinges |
| $q_{k,0}$ | position of origin atom of cluster $k$ | $3\times1$, Angstrom | hinge to parent originates here |
| $q_{k,b}$ | position of branching atom of cluster $k$ | $3\times1$, Angstrom | hinge to child originates here |
| $\delta E$ | relative total-energy fluctuation | kcal/mol | eq:11 accuracy metric |
| $E,E_k$ | total energy / total kinetic energy | kcal/mol | |
| $\phi,\psi$ | protein backbone dihedral angles | degrees | Ramachandran |
| $\chi_1$ | first side-chain dihedral angle | degrees | |

## Conventions
- **Base cluster** = designated cluster $(n{+}1)$; connected to inert frame by a virtual 6-DOF hinge (hinge matrix = $6\times6$ identity). Its orientation is integrated via quasi-coordinates (integrals of angular velocity) using Euler angles or a quaternion (Fincham's implicit quaternion algorithm).
- **Parent/child**: parent is on the path toward the base; tip clusters have no children.
- Single-DOF torsional hinges assumed throughout ($D$ diagonal, O(n) inversion).
- Forces from the Cartesian force field are used directly (per-atom spatial force = torque about cluster origin + net force); no explicit gradient of potential w.r.t. internal variables needed.
- Energy units: kcal/mol; length Angstrom; time fs. Reduced/rigid-body model: bond and angle DOF are frozen (covalent geometry fixed); their internal energy terms are NOT evaluated.
- Hydrogen mass in TAMD runs increased to 6.0 amu (mass-repartitioning-like trick to damp fast tip-group rotations).
- Softened vdW/electrostatics: interactions unchanged until energy exceeds a threshold, then switch to linear (soft-core) form.
