# Equations - Chen, Im, Brooks 2005 (TAMD / NEIMO + projected ICFF)

<!-- eq:1 -->
$$ \Omega(\theta)\ddot{\theta} + C(\theta,\dot{\theta}) = T(\theta) $$
- **what:** Equations of motion for an internal-variable (torsion-angle) system; solve for the internal accelerations $\ddot\theta$ given generalized forces.
- **symbols:** $\theta$ - $n$-vector of internal coordinates (torsion angles, rad); $\Omega(\theta)$ - $n\times n$ configuration-dependent mass matrix; $C(\theta,\dot\theta)$ - $n$-vector of Coriolis/centrifugal/gyroscopic + Cartesian-force terms; $T(\theta)$ - $n$-vector of generalized (applied) force; $n$ - number of torsional DOF ($\ll 3N-6$).

<!-- eq:2 -->
$$ V_{k} = \phi_{k+1,k}^{T} V_{k+1} + H_{k}^{T}\dot{\theta}_{k} $$
- **what:** Newton-Euler recursion for spatial velocity of cluster $k$ (base-to-tip), from parent velocity plus hinge contribution.
- **symbols:** $V_k=\mathrm{Col}[\omega_k,v_k]$ - 6-vector spatial velocity (angular $\omega_k$, linear $v_k$); $\phi_{k+1,k}$ - $6\times6$ spatial rigid-body transform between frames of cluster $k{+}1$ and $k$; $H_k^{T}=\mathrm{Col}[\hat h_k,0,0,0]$ - $6\times1$ hinge map, $\hat h_k$ unit vector along hinge $k$; $\dot\theta_k$ - hinge (torsion) rate; superscript $T$ = transpose.

<!-- eq:3 -->
$$ \alpha_k = \phi_{k+1,k}^{T}\alpha_{k+1} + H_k^{T}\ddot{\theta}_k + a_k $$
- **what:** Newton-Euler recursion for spatial acceleration of cluster $k$.
- **symbols:** $\alpha_k$ - 6-vector spatial acceleration (time derivative of $V_k$); $\ddot\theta_k$ - hinge angular acceleration; $a_k$ - spatial gyroscopic/Coriolis acceleration term (function of spatial velocities, inertia, hinge velocity).

<!-- eq:4 -->
$$ F_k = \phi_{k,k-1}F_{k-1} + M_k\alpha_k + b_k - F_k^{(c)} $$
- **what:** Newton-Euler recursion for hinge spatial force (tip-to-base), accumulating inertial and applied forces.
- **symbols:** $F_k$ - 6-vector hinge spatial force between clusters $k{+}1$ and $k$; $M_k$ - $6\times6$ spatial inertia about cluster origin; $b_k$ - Coriolis force term (function of velocities/inertia/hinge velocity); $F_k^{(c)}$ - effective Cartesian spatial force on cluster $k$ (torque about origin + net force from potential, combining all atoms).

<!-- eq:5 -->
$$ T_k = H_k F_k $$
- **what:** Projection of the hinge spatial force onto the allowed DOF of hinge $k$ gives the generalized torque; vanishes when all force-field forces are included in $F_k^{(c)}$.
- **symbols:** $T_k$ - generalized hinge force/torque along allowed DOF; $H_k$ - $1\times6$ hinge map (row).

<!-- eq:6 -->
$$ V = \Phi^{T} H^{T}\dot{\theta} $$
- **what:** Spatial-operator (stacked) form of the velocity recursion eq:2 over the whole chain.
- **symbols:** $\Phi$ - lower-triangular $6n\times6n$ spatial operator built from the $\phi_{k+1,k}$ blocks; $V$ - stacked $6n\times1$ spatial-velocity vector; $H^{T}=\mathrm{diag}\{H_1^{T},\dots,H_n^{T\}}$ - $6n\times n$ block-diagonal hinge operator; $\dot\theta$ - $n\times1$ hinge-rate vector.

<!-- eq:7 -->
$$ \Omega(\theta) \equiv H\,\Phi\,M\,\Phi^{T} H^{T} $$
- **what:** Newton-Euler operator factorization of the internal-coordinate mass matrix (from substituting the operator forms into eq:5).
- **symbols:** $M=\mathrm{diag}\{M_1,\dots,M_n\}$ - $6n\times6n$ block-diagonal spatial inertia; $H,\Phi$ as above; $\Omega$ - $n\times n$ mass matrix of eq:1.

<!-- eq:8 -->
$$ C(\theta,\dot{\theta}) \equiv H\Phi\big(M\Phi^{T}a + b - F^{(c)}\big) $$
- **what:** Operator expression for the velocity-dependent + applied-force term of eq:1.
- **symbols:** $a,b$ - stacked $6n$ gyroscopic/Coriolis vectors; $F^{(c)}$ - stacked $6n$ effective Cartesian spatial forces; result $C$ is $n\times1$.

<!-- eq:9 -->
$$ \Omega(\theta) = [\,1 + H\Phi K\,]\,D\,[\,1 + H\Phi K\,]^{T} $$
- **what:** Innovations operator factorization of the mass matrix into square, invertible factors; enables O(n) inversion.
- **symbols:** $K$ - Kalman-gain-like spatial operator from the recursion; $D$ - $n\times n$ articulated-body factor, positive-definite and diagonal when every hinge has a single DOF; $1$ - identity.

<!-- eq:10 -->
$$ \Omega^{-1} = [\,1 - H\Phi K\,]^{T} D^{-1} [\,1 - H\Phi K\,] $$
- **what:** Recursive O(n) inverse of the mass matrix (from the innovations factorization); solves eq:1 for $\ddot\theta$ without forming $\Omega$ explicitly.
- **symbols:** as eq:9; $D^{-1}$ trivial since $D$ is diagonal.
<!-- CHECK: raw eq (10) prints [1 - HΦK]^T on BOTH sides (double transpose), which is dimensionally inconsistent with eq:9. Standard spatial-operator-algebra (Jain/Rodriguez, ref 13) gives one transposed factor as written here: Ω^{-1} = [1 - HΦK]^T D^{-1} [1 - HΦK]. Verify factor ordering/transpose against ref 13 before implementing. -->

<!-- eq:11 -->
$$ \delta E = \frac{\langle E^{2}\rangle - \langle E\rangle^{2}}{\langle E_k\rangle} $$
- **what:** Relative total-energy-fluctuation metric used to assess integrator accuracy in NVE runs vs time step.
- **symbols:** $E$ - total energy; $E_k$ - total kinetic energy; $\langle\cdot\rangle$ - time (ensemble) average over the trajectory; $\delta E$ has units of energy (kcal/mol).

## Derivations note
The operator forms eq:6-eq:8 and the two factorizations eq:9-eq:10 are derived by stacking the per-cluster Newton-Euler recursions eq:2-eq:5 into spatial operators (Jain-Rodriguez spatial operator algebra, ref 12-14). Extension to branched trees: during tip-to-base recursion sum child contributions before proceeding; during base-to-tip recursion continue separately along each branch. These derivation steps are not reproduced here; see refs 12-14.
