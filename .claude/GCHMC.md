# The Gibbs scan

## Definitions

### **HAMILTONIAN GEOMETRY**

The space of all possible states `(p,q)` (position and momentum) is called the phase space. The Hamiltonian `Hamiltonian(q,p)` is the total energy written as a function of position and momentum and it dictates the complete motion through 2 equations:

- `dq/dt =  d(H)/dp`: Position is steered by how energy varies with momentum.
- `dp/dt = -d(H)/dq`: Momentum by how energy varies with position, but with a flipped sign.

The system does not move towards lower energy, but stays along iso-Hamiltonian levels. At every phase point, the symplectic form maps the energy gradient to a flow that conserves `H`.

### **MARKOV CHAIN**

The sampler produces `s0 -> s1 -> s2 -> s3 -> s...`. The next state depends on the current one: `prob(s(n+1)|s(n),s(n-1),...) = prob(s(n+1)|s(n))`. This is the Markov property. The Markov chain is the central object of the entire theory.

### **STATIONARY DISTRIBUTION**

A distribution `pi` is stationary if one step leaves it unchanged: `pi*TransitionKernel = pi`. If states are already distributed according to `pi`, another sampling step still gives `pi`.

### **DETAILED BALANCE**

Let $a$ and $b$ be two states and $T$ the transition kernel between them.
The detailed balance condition requires:

$$\pi(a) \cdot T(a \to b) = \pi(b) \cdot T(b \to a)$$

**Interpretation:** The probability current flowing from `a` to `b` equals the current flowing backward.
Detailed balance implies stationarity. Stationarity does not require detailed balance.

---

#### Expanded form

Each side decomposes as:

$$\pi(a) = p(q_a) \cdot \sqrt{|G(q_a)|}$$

$$\pi(b) = p(q_b) \cdot \sqrt{|G(q_b)|}$$

$$T(a \to b) = p(u_a) \cdot \hat{\phi}^t(a \to b)$$

$$T(b \to a) = p(u_b) \cdot \hat{\phi}^t(a \to b) \cdot A(a \to b)$$

where $A(a \to b)$ is the **accepting probability** of the move from `a` to `b`.

So detailed balance becomes:

$$p(q_a) \cdot \sqrt{|G(q_a)|} \cdot p(u_a) \cdot \hat{\phi}^t
\;=\;
p(q_b) \cdot \sqrt{|G(q_b)|} \cdot p(u_b) \cdot \hat{\phi}^t \cdot A(a \to b)$$

Note: the Hamiltonian flow $\hat{\phi}^t$ cancels on both sides by Liouville's theorem (Hamiltonian flow is volume-preserving).

---

#### Mass metric tensor $\sqrt{|G(q)|}$

$\sqrt{|G(q)|}$ is the square root of the determinant of the mass-weighted metric tensor:

$$\sqrt{|G(q)|} = \sqrt{\left|J_{q}^T(\mathbf{r})\, M\, J_{q}(\mathbf{r})\right|}$$

Its form depends on the coordinate system.

**Case 1 — Cartesian coordinates:**

$$\sqrt{|G(q)|} = 1$$

**Case 2 — Internal (BAT) coordinates:**

$$Q_\text{BAT}(\beta) = \Xi \frac{(2\pi kT)^{3N/2}}{h^{3N}}
\int_{\mathcal{B}} e^{-\beta U(q)} \cdot
\left| J_{q}^T(\mathbf{r})\, M\, J_{q}(\mathbf{r}) \right|^{1/2}
\,d^{3N}q$$

**Case 3 — Constrained internal coordinates:**

$$Q(\beta) = \Xi \frac{(2\pi kT)^{\mathcal{N}/2}}{h^{\mathcal{N}}}
\int_{\mathcal{B}_f} e^{-\beta U(q)} \cdot
\frac{\left| J_{q_f}^T(\mathbf{r})\, M\, J_{q_f}(\mathbf{r}) \right|^{1/2}}
     {\left| J_{q}^T(\mathbf{r})\, M\, J_{q}(\mathbf{r}) \right|^{1/2}}
\,d^{\mathcal{N}}q$$

Where $$U(q)$$ is the potential energy function calculated on the full-atom molecular system.

### **ERGODICITY (OPTIONAL)**

In Robosample, ergodicity implies sampling all Q's i.e degrees of freedom available for the entire molecular system.

### **HAMILTONIAN**

Let the Hamiltonian energy function be `Hamiltonian(q,p) = U(q) + KE(q,p) + U_Fixman(q) + U_Fixman_External(q)`.

For generalized coordinates:

- `U(q)` is the potential energy computed on full-atom model.
- `KE(q,p) = 1/2* p^T *M(q)^-1* p` is the kinetic energy.
- `U_Fixman(q)` is the artificial movement of the metric tensor inside the exponential: -kt * log(G(q)). The Fixman term for the external degrees of freedom is called `U_Fixman_External`. This term is null in Cartesian coordinates.
- `U_Fixman_External(q)` is `ln( sin^2(gamma2) )` summed over all `Free`-rooted molecules. Note that `gamma2` is the pitch of each root's absolute orientation. This term is null in Cartesian coordinates.

For cartesian coordinates, `KE = 1/2 * p^T * M^-1 * p` where `p^T` is the transpose momenta vector.
The Hamiltonian determines the equilibrium distribution: `pi(q,p) proportional to e^(-beta*H(q,p))`

### **ACCEPTANCE**

A proposal is a candidate move: `s -> s'`. Example: draw momenta, integrate dynamics for 100 steps. The proposal is not yet accepted.

The Acceptance probability (Metropolis Hastings): `A=min(1,e^(-beta*deltaH))`. If `deltaH == 0`, then accept with probability `1`. If energy increased, accept probabilistically. Acceptance restores detailed balance.

Simple HMC example: draw momenta, integrate Hamilton's equations, compute `deltaH`, accept or reject. This produces a Markov kernel `T(s->s')`

### **REVERSIBILITY**

A proposal is reversible if the reverse path exists:

- Forward: `(q, p)` -> `(q', p')`.
- Reverse: `(q',-p)` -> `(q, -p)`.

Without reversibility, the backward transition probability cannot be computed.

### **VOLUME PRESERVATION**

The proposal must preserve phase-space volume `(dq,dp)`. No compression. No expansion. Otherwise the proposal density changes. The Jacobian would appear in the acceptance probability.

### **SYMPLECTICITY**

Symplectic integrators preserve Hamiltonian geometry. Consequences: phase-space volume preservation, near-energy conservation, shadow Hamiltonian. Energy errors remain bounded: `H(t) - H(0) = O(dt^2)`. This gives high acceptance rates. Symplecticity improves efficiency. It is not strictly required for correctness.

### **GIBBS SAMPLING**

A Gibbs step updates only part of the variables. Examples:

- `q=(q1, q2, q3)`. Update only `q2`. The other coordinates remain fixed.
- torsional dynamics: only torsions move; bond lengths and bond angles are fixed. This constitutes a torsional Gibbs block.

A Gibbs block defines: mobile coordinates, frozen coordinates, mass matrix, Hamiltonian, timestep.

### **GENERALIZED COORDINATES**

Instead of Cartesian coordinates `(x,y,z)`, we use `q=(bond_length, bond_angle, torsion/dihedral)`. The mass matrix then becomes `M(q)`. Based on the accessible degrees of freedom, the space can be:

- Flat (linear Taylor):
  - `Pin` torsions live on the circle `S^1`, which has **zero intrinsic curvature** (a circle is a line made periodic, locally indistinguishable from `R`) and is integrated as an unwrapped real.
  - `Translation` and `Free`  coordinates lie in the flat, zero-curvature space `R^3`.
- Curved path (exponential map): The **orientation** of a `Free` or `Ball` body, stored as a unit quaternion living on **S^3, the unit 3-sphere in R^4** which is a compact manifold of **constant positive curvature** that **double-covers the rotation group SO(3)** (`q` and `-q` are the same physical rotation).

From a physics point of view, bonds and angles vary rapidly and contribute little to the overall RMSD. On the contrary, torsions are softer modes and heavily influence conformational transitions, but these transitions are dependent on all degrees of freedom. For example, torsional dynamics cannot sample alone cis to trans isomerization of alanine dipeptide which is gated by a 1-4 clash that requires angle relaxation. The transition occurs only when this torsions are sampled together with bond lengths.

### **EQUIPARTITION**

The correct momentum distribution is `p ~ N(0, R*T*M(q))`, not `p ~ N(0, I)`, called the Maxwell-Boltzmann distribution.

### **CONSTRAINTS**

Constraints define a manifold `sigma(q)=0`. Examples: ring closure. SHAKE projects positions. RATTLE projects velocities.

Blocks are executed in a chosen order (a deterministic systematic scan) or randomized (random scan). The order is a free choice, not a constraint of the method: each block individually leaves `pi` invariant, and a composition of `pi`-invariant kernels leaves `pi` invariant regardless of order. A systematic scan is generally not reversible (its time-reversal is the reversed block order). A random scan that picks a block from a fixed distribution is reversible.  Either way `pi` is preserved. **Stationarity (`pi*TransitionKernel = pi`), not reversibility, is the requirement**.

Although torsional dynamics is the most efficient way to exploit the internal-coordinate representation -- and the reason the representation exists -- it is **not** mandatory; a scan may mix Cartesian, torsional, and other generalized-coordinate blocks freely.

- **Overlap.** The restricted spaces of different blocks may overlap. Redundant updates do not break stationarity, provided each block individually leaves pi invariant.
- **Ergodicity.** Across one full round (all blocks), every mobile coordinate in the block set has been **offered** a move. Ergodicity is over the **sampled (restricted) subspace** and requires both that the blocks collectively cover that subspace and that each block achieves nonzero acceptance. Being offered a move is **not** the same as accepting one: a block stuck near zero acceptance contributes nothing even though aggregate statistics look healthy. Coordinates that are frozen in **every** block (e.g. bonds and angles in torsion-only operation) are not sampled; the sampled distribution is then the (Fixman-corrected) marginal over the union of the blocks' restricted spaces.
- **Per-block diagnostics.** Acceptance rates are reported **per block** specifically to detect an under-sampled subspace, which would otherwise be masked by healthy global energy and acceptance statistics.

## The HMC move inside a block

### Momentum resampling and the internal-coordinate equipartition principle

At the start of each block the momenta are **fully** resampled (not partial): `p ~ N(0, R*T*M(q))`. Equivalently, `u ~ N(0, R*T*M(q)^-1)`.

This is a Gibbs step on the momenta (always accepted). It realizes the **internal-coordinate equipartition principle** : in Cartesian coordinates every DOF carries `(1/2) * R * T` of kinetic energy, but in generalized coordinates with a configuration-dependent metric `M(q)` this does not hold per coordinate. The *generalized* equipartition theorem, `< p_i (dH/dp_j) > = R * T * delta_ij`, still holds, which forces the momentum covariance to be`< p p^T > = R * T * M(q)`. Drawing `p ~ N(0, R*T*M(q))` is exactly this distribution, so temperature is assigned consistently despite the non-constant metric. The draw is matrix-free via the articulated `sqrt(M(q)^-1)` operator (drawing velocities; momenta follow as `p = M(q)*u`), so `det M(q)` never appears explicitly in the proposal; its configuration dependence is bookkept by the Fixman term in the acceptance.

In a ring-closure block the freshly drawn momenta are then RATTLE-projected onto the velocity constraint surface `G * M^-1 * p = 0` before integration begins, so the initial velocity is consistent with the closed ring. That projection is exactly what injects the `det(G M^-1 G^T)^(-1/2)` factor into the configurational marginal, the factor the loop-closure term of `F` cancels. For acyclic molecules there is no projection and no such factor.

### Propagation under the bare potential `U`

Dynamics are integrated with the **fixed-step Verlet integrator** (internal and Cartesian alike), using forces derived from the bare potential `U` only. Concretely, the integrator's forces are the per-atom Cartesian forces `-grad(U)` from OpenMM, reduced to per-body spatial forces for the articulated solver; the gradients of the Fixman potential `U_Fixman` and the Jacobian `U_Jacobian` (the "torques") are not added to these forces. So `U` is what generates the trajectory; `U_Fixman` and `U_Jacobian` enter only the acceptance. After each step, constraints (if any) are projected onto the manifold.

The integrator is a second-order, semi-explicit predictor-corrector, not an explicit kick-drift-kick leapfrog. This distinction is load-bearing and is the source of the integrator's actual conservation properties: the position is advanced by an explicit second-order Taylor step using the start-of-step acceleration, and the velocity is advanced by an **implicit trapezoidal corrector** solved by functional iteration. The two coincide with textbook velocity Verlet only for a separable, constant-mass Hamiltonian (the Cartesian case); for the configuration-dependent metric `M(q)` and velocity-dependent (Coriolis/gyroscopic) forces of internal coordinates they do not, which is exactly where the symplecticity caveat applies.

The integrator is run with a **fixed step**, making it symplectic and reversible for separable systems and the step is **taken unconditionally** since there is no step-size adaptation. Thus, a non-converged corrector does not shrink `dt` or reject the step. It is the trajectory-level Metropolis test, not per-step control, that supplies correctness.

The per-step sequence (one fixed-step Verlet step of size `dt`; the SHAKE/RATTLE projections are no-ops on acyclic molecules). `qdot0 = N(q0) * u0` for quaternion DOF and `u0` otherwise; `qDotDot(q)` is the corresponding q-acceleration; `udot = M^-1 * (f - f_bias)` is one O(n) forward-dynamics sweep, with `f_bias` carrying the gyroscopic and Coriolis terms:

#### 1. Position (2nd-order Taylor expansion)

- `q <- q0 + dt*qdot0 + (dt^2/2)*qdotdot0` - not applied to quaternions.

---

#### 2. Quaternion Rotation (exponential map update)

- `q <- advanceQuatExp(wHalf, dt) (x) q0`
- where `wHalf = u0 + (dt/2)*udot0`
- then normalize:
  - `q <- q / |q|`

---

#### 3. SHAKE (position constraint projection)

- Project positions to satisfy constraints `sigma(q) = 0`
- Followed by kinematic refresh after projection

---

#### 4. Velocity Predictor (forward Euler step)

- `u <- u0 + dt*udot0` where forces are evaluated from `-grad(U) at q1`

---

#### 5. Velocity Corrector (trapezoidal iteration)

- Iterate up to 10 times:

  - `u <- u0 + (dt/2)*(udot0 + udot1)`
  - recompute `udot1`

- Convergence criterion:
  - `||du|| / ||u|| <= tol`

- where:
  - `tol = min(1e-4, 0.1*accuracy)`
  - where `accuracy` is a user-defined parameter.

---

#### 6. RATTLE (velocity constraint projection)

The RATTLE velocity correction uses fixed-point iteration to enforce `G * M^-1 p = 0`. If the iteration fails to converge within the allowed iterations, the best available velocity estimate is accepted. Free-body quaternions are not advanced by the linear Taylor formula of step 1; they are advanced by an **exact exponential-map rotation** (step 2): the increment quaternion `advanceQuatExp(wHalf, dt) = (cos(theta), sin(theta) * wHalf/|wHalf|)` with `theta = (1/2)|wHalf| dt` and the midpoint angular velocity `wHalf = u0 + (dt/2) udot0`, applied as a left Hamilton product onto q0, followed by renormalization that mops up only rounding. This preserves unit norm and avoid dependence on `qddot`. The resulting constrained integrator remains second-order accurate globally, with third-order local error estimates that are computed but unused in the fixed-step sampler.

## Acceptance

The proposal is accepted or rejected by an MH test on the full Hamiltonian function H, using the **exact** dH between the pre- and post-trajectory states: `probability min(1, exp(-beta * dH))` where `dH = H_new - H_old`.

HMC (Duane, Kennedy, Pendleton & Roweth 1987) admits a **distinct guidance and acceptance Hamiltonian**. The *guidance* Hamiltonian generates the trajectory (here: `U` forces); the *acceptance* Hamiltonian is the one whose Boltzmann distribution is actually sampled, used in the MH test. Duane et al.'s key observation is that the trajectory generator need **not** equal the target: as long as the proposal map is reversible and volume-preserving and the MH test uses the *exact* acceptance Hamiltonian, the move samples the acceptance distribution exactly, whatever the guidance was. Two consequences used here:

- The Fixman torque is unnecessary for correctness because `U_Fixman` enters only the acceptance and the integrator need not compute `grad(U_Fixman)`
- An efficiency-oriented guidance potential is admissible: a cheaper or smoother potential (e.g. an internal-coordinate force field) may be used to *guide* proposals while the exact atomistic `H` is retained for acceptance, with no bias.

The holonomic constrained proposal map must be time-reversible and volume-preserving with respect to the constrained measure (or carry its Jacobian in the ratio); with exact `dH` of the acceptance Hamiltonian, the move is then exactly `pi`-stationary. Symplecticity is a strengthening of volume preservation that additionally yields a conserved shadow Hamiltonian and hence bounded energy error (an efficiency property, not a correctness one). For the separable Cartesian case, exact-arithmetic Verlet is reversible, volume-preserving, and symplectic. For the non-separable internal-coordinate case the truncated corrector and non-converged RATTLE acceptance can break reversibility and volume preservation themselves (the correctness properties) unless the truncation is shown to be a time-symmetric, volume-preserving involution. Until then, π-stationarity of the internal-coordinate move is an assumption, not a theorem.

A larger `dt` is a coarser Verlet map. Provided the map is reversible and volume-preserving and the **exact dH** is used in the MH test, the move is exactly `pi`-stationary at any `dt`. Assigning small `dt` to stiff blocks and large `dt` to soft blocks therefore introduces **no bias**, only a change in efficiency. This is the sense in which Robosample is a multiscale HMC sampler.

Reversibility is configuration dependent. It is **not** the curvature of S^3 that varies since that is **constant**, so the orientation manifold's intrinsic geometry is identical everywhere. What varies is the **mass metric M(q)** and the **local force stiffness d^2(U) / d(q)^2** (a near-clash or compressed region is far stiffer than an open one). A `dt` that is perfectly symplectic in an open conformation can become non-symplectic and begin pumping energy when the chain visits a stiffer geometry.

Any method that samples in a reduced (e.g. torsional) space while holding bonds/angles fixed must confront one fact: **the equilibrium distribution of the soft coordinates in a model with rigid bonds/angles is not the same as their marginal in the fully flexible model**. It has two distinct pieces, frequently conflated:

1. **A measure (kinetic / metric) piece.** Marginalizing the Gaussian momenta of a constrained system leaves a configuration-dependent factor `det(M(q))^(1/2)`. Uncorrected, the sampler targets `exp(-beta*U)*det(M)^(1/2)` instead of `exp(-beta*U)`. This is what the **Fixman potential** removes.
2. **A potential-of-mean-force (PMF) piece.** In a flexible molecule, bond angles *relax* in response to the torsional configuration (an angle widens to relieve a 1-4 clash at a given torsion). Freezing the angle removes that relaxation. This is an energetic effect: the location of the bond/angle minimum moves with torsion -- and **no determinant correction captures it**. It is recovered by **mobilizing the coupled hard DOF** through the
**fully flexible Cartesian world**, which alone covers every DOF and so suffices for correctness. Robosample relaxes them in **separate Gibbs worlds** corrected by an exact MH test, so no independence assumption is needed.

Running torsion-angle dynamics with rigid covalent geometry **grossly distorts** the torsional energy surface of Cartesian-parameterized force fields, but that can be repaired by building a specialized **internal-coordinate force field (ICFF)** (modified torsion (CMAP-style) terms plus softened van der Waals and
electrostatics) that approximates the source Cartesian field *without* a compensating potential. Robosample needs no such correction, for two independent reasons:

- **The proposal is corrected by an exact acceptance test.** Moreover, in Robosample, the rigid-geometry dynamics is only a **proposal**; the MH test uses the full atomistic `U` (plus Fixman), so the sampled distribution is exact regardless of how distorted the proposal forces are (the distinct guidance/acceptance design). At worst a poor proposal lowers acceptance; it cannot bias the result.
- **The frozen-geometry PMF distortion is integrated out.** The residual piece-2 distortion from holding bonds/angles fixed is removed not by patching `U` but by **mobilizing those coordinates** in the `Cartesian` worlds of the scan.
