# Equations — Tao et al. 2012, SHAPE rigid-structure constraint

Reduced/CHARMM MD units. `t` = time, `Δt` = time step, `M` = number of atoms in
the rigid structure. Superscript `b` = body-fixed local (COM-origin) coordinate,
`non` = unconstrained (free MD) value, `rig` = rigid-constrained value,
`(n)` = iteration index.

<!-- eq:1 -->
$$ \mathbf{r}_{COM}(t) = \frac{\sum_{i=1}^{M} m_i \mathbf{r}_i(t)}{\sum_{i=1}^{M} m_i} $$
- **what:** center of mass of the rigid structure at time t.
- **symbols:** r_COM - COM position (R^3); m_i - mass of atom i; r_i(t) - global position of atom i (R^3); M - number of atoms in rigid body.

<!-- eq:2 -->
$$ \mathbf{r}_{COM}^{non}(t+\Delta t) = \frac{\sum_{i=1}^{M} m_i \mathbf{r}_i^{non}(t+\Delta t)}{\sum_{i=1}^{M} m_i} $$
- **what:** COM of the unconstrained (free MD) positions at t+Δt.
- **symbols:** r_COM^non - COM of unconstrained positions; r_i^non(t+Δt) - free-MD updated global position of atom i.

<!-- eq:3 -->
$$ \mathbf{r}_{i}^{b}(t) = \mathbf{r}_{i}(t) - \mathbf{r}_{COM}(t) $$
- **what:** body-fixed local coordinate of atom i at time t (subtract COM, no rotation).
- **symbols:** r_i^b(t) - body-fixed position at t (R^3).

<!-- eq:4 -->
$$ \mathbf{r}_{i}^{non, b}(t + \Delta t) = \mathbf{r}_{i}^{non}(t + \Delta t) - \mathbf{r}_{COM}^{non}(t + \Delta t) $$
- **what:** body-fixed local coordinate of the unconstrained position at t+Δt.
- **symbols:** r_i^non,b(t+Δt) - body-fixed unconstrained position (R^3).

<!-- eq:5 -->
$$ L^{non}\left(t + \frac{\Delta t}{2}\right) = \sum_{i=1}^{M} \mathbf{r}_{i}^{non, b}\left(t + \frac{\Delta t}{2}\right) \times \left(m_{i}\frac{\mathbf{r}_{i}^{non, b}(t + \Delta t) - \mathbf{r}_{i}^{b}(t)}{\Delta t}\right) $$
- **what:** angular momentum of the unconstrained structure at half step t+Δt/2 (velocity by finite difference of body-fixed coords).
- **symbols:** L^non - unconstrained angular momentum (R^3); × - cross product; velocity v_i^b = (r_i^non,b(t+Δt) - r_i^b(t))/Δt.

<!-- eq:6 -->
$$ I = \begin{bmatrix} I_{xx} & I_{xy} & I_{xz} \\ I_{yx} & I_{yy} & I_{yz} \\ I_{zx} & I_{zy} & I_{zz} \end{bmatrix} $$
- **what:** moment-of-inertia tensor (3×3) of the M atoms in body-fixed local coords r_i^b(t).
- **symbols:** I - inertia tensor; components I_ab = Σ_i m_i (‖r_i^b‖^2 δ_ab − r_i,a^b r_i,b^b). <!-- CHECK: raw eq. (6) OCR showed I_xy and I_yz in the top row (I_xy,I_yz); corrected to symmetric I_xx,I_xy,I_xz assuming standard symmetric inertia tensor -->

<!-- eq:7 -->
$$ L^{rig}\left(t + \frac{\Delta t}{2}\right) = I \cdot \omega^{rig}\left(t + \frac{\Delta t}{2}\right) $$
- **what:** angular momentum of the rigid structure from inertia tensor and angular velocity.
- **symbols:** L^rig - rigid angular momentum (R^3); ω^rig - angular velocity vector = (ω_x, ω_y, ω_z) (R^3).

<!-- eq:8 -->
$$ L^{non}\left(t + \frac{\Delta t}{2}\right) = L^{rig}\left(t + \frac{\Delta t}{2}\right) $$
- **what:** the constraint solved: rigid angular momentum must equal the unconstrained angular momentum.
- **symbols:** equality of L^non (eq. 5) and L^rig (eq. 7).

<!-- eq:10 -->
$$ \omega^{rig,(1)}\left(t + \frac{\Delta t}{2}\right) = \mathbf{I}^{-1} \cdot \sum_{i=1}^{M} \mathbf{r}_{i}^{non, b}\left(t + \frac{\Delta t}{2}\right) \times \left(m_{i} \frac{\mathbf{r}_{i}^{non, b}(t + \Delta t) - \mathbf{r}_{i}^{b}(t)}{\Delta t}\right) $$
- **what:** first estimate of the rigid angular velocity = I^-1 · L^non.
- **symbols:** ω^rig,(1) - first-iteration angular velocity; I^-1 - inverse inertia tensor.

<!-- eq:11 -->
$$ \hat{\omega} = \begin{bmatrix} 0 & -\omega_z & \omega_y \\ \omega_z & 0 & -\omega_x \\ -\omega_y & \omega_x & 0 \end{bmatrix} $$
- **what:** skew-symmetric matrix built from angular velocity vector ω (hat map).
- **symbols:** ω̂ - skew-symmetric 3×3; ω_x,ω_y,ω_z - components of ω. Satisfies ω̂ v = ω × v.

<!-- eq:12 -->
$$ \mathbf{R} = e^{\hat{\theta}} $$
- **what:** rotation matrix = matrix exponential of the skew-symmetric angle matrix; θ̂ = ω̂ Δt, θ = ‖ω‖Δt.
- **symbols:** R - rotation matrix (SO(3)); θ̂ - skew matrix of rotation-angle vector θ = ω Δt.

<!-- eq:13 -->
$$ e^{\hat{\theta}} = \mathbf{1} + \frac{\hat{\theta}}{\|\theta\|} \sin(\|\theta\|) + \frac{\hat{\theta}^2}{\|\theta\|^2} (1 - \cos(\|\theta\|)) $$
- **what:** Rodrigues' formula for the rotation matrix (numerically unstable when ‖θ‖ small).
- **symbols:** 1 - identity matrix; ‖θ‖ - Euclidean norm of θ = rotation angle; θ̂^2 - matrix square of θ̂.

<!-- eq:14 -->
$$ e^{\hat{\theta}} = \mathbf{1} + \left(1 - \frac{\|\theta\|^2}{3!} + \frac{\|\theta\|^4}{5!} - \cdots\right)\hat{\theta} + \left(\frac{1}{2!} - \frac{\|\theta\|^2}{4!} + \frac{\|\theta\|^4}{6!} - \cdots\right)\hat{\theta}^2 $$
- **what:** Taylor expansion of the rotation matrix for small ‖θ‖ (numerically stable form of eq. 13).
- **symbols:** coefficients are the series for sin(‖θ‖)/‖θ‖ and (1−cos‖θ‖)/‖θ‖^2.

<!-- eq:15 -->
$$ \mathbf{r}_{i}^{rig,\,b(n)}(t+\Delta t) = \mathbf{R}^{(n)} \cdot \mathbf{r}_{i}^{b}(t) $$
- **what:** update body-fixed rigid coordinates at t+Δt by rotating time-t body-fixed coords with R^(n).
- **symbols:** r_i^rig,b(n)(t+Δt) - rigid body-fixed position, iteration n; R^(n) - rotation matrix at iteration n.

<!-- eq:16 -->
$$ L^{rig, (n)}\left(t + \frac{\Delta t}{2}\right) = \sum_{i=1}^{M} \mathbf{r}_{i}^{rig, b, (n)}\left(t + \frac{\Delta t}{2}\right) \times \left(m_{i}\frac{\mathbf{r}_{i}^{rig, b, (n)}(t + \Delta t) - \mathbf{r}_{i}^{b}(t)}{\Delta t}\right) $$
- **what:** angular momentum of the rigid structure computed from updated coordinates at iteration n.
- **symbols:** L^rig,(n) - iteration-n rigid angular momentum; ṙ_i^rig,b,(n) = (r_i^rig,b,(n)(t+Δt) − r_i^b(t))/Δt.

<!-- eq:17 -->
$$ \mathbf{r}_{i}^{rig, b, (n)} \left( t + \frac{\Delta t}{2} \right) = (\mathbf{R}^{(n)})^{\frac{1}{2}} \mathbf{r}_{i}^{b}(t) $$
- **what:** approximate the half-step (t+Δt/2) coordinate using the matrix square root of R^(n).
- **symbols:** (R^(n))^{1/2} - matrix square root of the rotation matrix (half rotation).

<!-- eq:18 -->
$$ L^{non}\left(t + \frac{\Delta t}{2}\right) \approx (\mathbf{R}^{(n)})^{\frac{1}{2}} \cdot \mathbf{L}'^{non}, \qquad \mathbf{L}'^{non} = \sum_{i=1}^{M} \mathbf{r}_{i}^{non, b}(t) \times \left(m_{i} \frac{\mathbf{r}_{i}^{non, b}(t + \Delta t) - \mathbf{r}_{i}^{b}(t)}{\Delta t}\right) $$
- **what:** factor the half-step rotation out of the unconstrained angular momentum; L'^non uses time-t body-fixed coords.
- **symbols:** L'^non - reduced (unrotated) unconstrained angular momentum (R^3).

<!-- eq:19 -->
$$ L^{rig,(n)}\left(t + \frac{\Delta t}{2}\right) \approx (\mathbf{R}^{(n)})^{\frac{1}{2}} \cdot \mathbf{L}'^{rig,(n)}, \qquad \mathbf{L}'^{rig,(n)} = \sum_{i=1}^{M} \mathbf{r}_{i}^{rig,(n), b}(t) \times \left(m_{i} \frac{\mathbf{r}_{i}^{rig,(n), b}(t + \Delta t) - \mathbf{r}_{i}^{b}(t)}{\Delta t}\right) $$
- **what:** factor the half-step rotation out of the rigid angular momentum; L'^rig,(n) is the reduced (unrotated) form.
- **symbols:** L'^rig,(n) - reduced rigid angular momentum at iteration n (R^3).

<!-- eq:20 -->
$$ \boldsymbol{I} \cdot \omega^{rig,(n+1)} \left( t + \frac{\Delta t}{2} \right) = \boldsymbol{L}^{non} \left( t + \frac{\Delta t}{2} \right) $$
- **what:** target condition for the next iteration's angular velocity (reproduce L^non).
- **symbols:** ω^rig,(n+1) - next-iteration angular velocity.

<!-- eq:21 -->
$$ I \cdot \omega^{rig,(n+1)} \left( t + \frac{\Delta t}{2} \right) = I \cdot \omega^{rig,(n)} \left( t + \frac{\Delta t}{2} \right) + L^{non} \left( t + \frac{\Delta t}{2} \right) - L^{rig,(n)} \left( t + \frac{\Delta t}{2} \right) $$
- **what:** combine eq. (20) and (7): fixed-point update of the angular momentum residual.
- **symbols:** residual = L^non − L^rig,(n) drives the correction.

<!-- eq:22 -->
$$ \omega^{rig,(n+1)}\left(t + \frac{\Delta t}{2}\right) = \omega^{rig,(n)}\left(t + \frac{\Delta t}{2}\right) + \mathbf{I}^{-1} \cdot (\mathbf{R}^{(n)})^{\frac{1}{2}} \cdot \left(\mathbf{L}'^{non} - \mathbf{L}'^{rig,(n)}\right) $$
- **what:** angular velocity iteration update using the matrix square root of R^(n).
- **symbols:** exact update before the square-root approximation (eq. 25).

<!-- eq:23 -->
$$ \mathbf{R}^{(n)} = \mathbf{1} + \delta \mathbf{R}, \qquad \delta \mathbf{R} = \mathbf{R}^{(n)} - \mathbf{1} $$
- **what:** small-rotation decomposition of R^(n) for small time step.
- **symbols:** δR - small deviation of R^(n) from identity.

<!-- eq:25 -->
$$ (\mathbf{R}^{(n)})^{\frac{1}{2}} \approx \frac{1}{2} (\mathbf{1} + \mathbf{R}^{(n)}) $$
- **what:** first-order approximation of the matrix square root (neglect O((δR)^2)).
- **symbols:** used to make the update eq. (26) cheap.

<!-- eq:26 -->
$$ \omega^{rig,(n+1)}\left(t + \frac{\Delta t}{2}\right) = \omega^{rig,(n)}\left(t + \frac{\Delta t}{2}\right) + I^{-1} \cdot \frac{(\mathbf{1} + \mathbf{R}^{(n)})}{2} \cdot (\mathbf{L}'^{non} - \mathbf{L}'^{rig,(n)}) $$
- **what:** practical angular velocity updating scheme (eq. 22 with square-root approximation eq. 25). Iterate to convergence; ~3 iterations gives double precision.
- **symbols:** ω^rig,(n+1) - updated angular velocity for iteration n+1.

<!-- eq:27 -->
$$ \mathbf{r}_{i}^{rig, (n)}(t + \Delta t) = \mathbf{r}_{i}^{rig, b(n)}(t + \Delta t) + \mathbf{r}_{COM}^{non}(t + \Delta t) $$
- **what:** convert converged body-fixed rigid coordinates back to the global frame (add COM back).
- **symbols:** r_i^rig,(n)(t+Δt) - final global rigid coordinate = desired MD trajectory position.

<!-- eq:28 -->
$$ \sum_{i=1}^{M} \mathbf{r}_{i}^{b}\left(t + \frac{\Delta t}{2}\right) \times \left(m_{i}\frac{\mathbf{r}_{i}^{non,b}(t+\Delta t) - \mathbf{r}_{i}^{b}(t)}{\Delta t}\right) = \sum_{i=1}^{M} \mathbf{r}_{i}^{b}\left(t + \frac{\Delta t}{2}\right) \times \left(m_{i}\frac{\mathbf{r}_{i}^{rig,b}(t+\Delta t) - \mathbf{r}_{i}^{b}(t)}{\Delta t}\right) $$
- **what:** explicit form of the constraint eq. (8) that the iteration solves — equal angular momentum of unconstrained and rigid displacements.
- **symbols:** r_i^non,b, r_i^rig,b - body-fixed positions at t+Δt without/with rigid constraint.
