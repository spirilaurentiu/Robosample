# Equations — Jain, Vaidehi, Rodriguez 1993, O(N) recursive MD dynamics

Serial-chain convention: clusters numbered $1..n$ from **tip to base**; cluster
$k{+}1$ is the parent of cluster $k$; hinge $k$ couples clusters $k{+}1$ and $k$.
Inertial frame index is $n{+}1$; $f(0)=0$. Shorthand $x(k)\equiv
x(\mathcal{O}_k^-)$. `*` = matrix transpose (adjoint of a spatial operator).

<!-- eq:1.1 -->
$$ \mathcal{M}(\theta)\ddot{\theta} + \mathcal{C}(\theta, \dot{\theta}) = T(\theta) $$
- **what:** internal-variable equations of motion for a tree-topology molecular system.
- **symbols:** $\theta$ - $\mathcal{N}$-vector of generalized (internal) coordinates; $\ddot\theta$ - generalized accelerations; $T$ - $\mathcal{N}$-vector of generalized forces; $\mathcal{M}(\theta)$ - $\mathcal{N}\times\mathcal{N}$ symmetric positive-definite mass matrix (non-diagonal, nonlinear in $\theta$); $\mathcal{C}(\theta,\dot\theta)$ - $\mathcal{N}$-vector of Coriolis/velocity-dependent nonlinear forces; $\mathcal{N}$ - number of dof.

<!-- eq:PE -->
$$ PE = \mathcal{P}[\theta] + \mathcal{P}[x] $$
- **what:** potential energy split into internal-coordinate and Cartesian-coordinate parts.
- **symbols:** $\mathcal{P}[\theta]$ - internal-coordinate potential; $\mathcal{P}[x]$ - Cartesian-coordinate potential; $x$ - Cartesian atom positions.

<!-- eq:2.1 -->
$$ \mathcal{M}_c \ddot{x} = f $$
- **what:** free-atom Cartesian Newton equations of motion.
- **symbols:** $\mathcal{M}_c$ - $3n\times 3n$ diagonal mass matrix (atomic masses); $x\in\Re^{3n}$ - Cartesian atom coordinates; $f$ - $3n$-vector of inter-atomic forces; $n$ - number of atoms.

<!-- eq:2.2 -->
$$ A\dot{x} = B $$
- **what:** $m$ hard constraints as instantaneous linear constraints on atomic velocities.
- **symbols:** $A\in\Re^{m\times 3n}$ - constraint matrix (config-dependent); $B\in\Re^m$ - constraint vector (config-dependent); $m$ - number of constraints.

<!-- eq:2.3 -->
$$ \mathcal{M}_r \ddot{x}_r = f_r $$
- **what:** dimensionally reduced ODE after eliminating constrained velocity components.
- **symbols:** $\mathcal{M}_r$ - $\mathcal{N}\times\mathcal{N}$ reduced mass matrix (non-diagonal); $x_r$ - reduced generalized coordinates; $\mathcal{N}\triangleq 3n-m$ - remaining dof; $f_r$ - reduced force vector.

<!-- eq:4.2 -->
$$ V(k) = \phi^*(k+1, k)\, V(k+1) + H^*(k)\, \dot{\theta}(k) $$
- **what:** recursive spatial velocity of frame $\mathcal{O}_k^-$ (base-to-tip sweep).
- **symbols:** $V(k)$ - 6-dim spatial velocity at $\mathcal{O}_k^-$; $\phi^*(k{+}1,k)$ - transpose of $6\times 6$ spatial transform from parent frame; $H^*(k)$ - $6\times m$ hinge matrix; $\dot\theta(k)$ - $m$-vector hinge velocity.

<!-- eq:4.3 -->
$$ \alpha(k) = \dot{V}(k) = \phi^*(k+1, k)\, \alpha(k+1) + H^*(k)\, \ddot{\theta}(k) + a(k) $$
- **what:** recursive spatial acceleration of frame $\mathcal{O}_k^-$ (time derivative of eq:4.2).
- **symbols:** $\alpha(k)$ - 6-dim spatial acceleration at $\mathcal{O}_k^-$; $\ddot\theta(k)$ - hinge acceleration; $a(k)$ - Coriolis spatial acceleration (eq:4.4).

<!-- eq:4.4 -->
$$ a(k) = \begin{pmatrix} 0 \\ \tilde{\omega}(k+1)\,[v(k) - v(k+1)] \end{pmatrix} + \begin{pmatrix} \tilde{\omega}(k) & 0 \\ 0 & \tilde{\omega}(k) \end{pmatrix} H^*(k)\, \dot{\theta}(k) $$
- **what:** Coriolis spatial acceleration term $a(k)$. Since $H(k)$ is constant in $\mathcal{O}_k^\pm$, $\omega(k{+}1)$ may replace $\omega(k)$ in the first block.
- **symbols:** $\tilde\omega$ - $3\times 3$ cross-product (skew) tensor of angular velocity $\omega$; $v(k)$ - linear velocity part of $V(k)$; the second matrix is $6\times 6$ block-diagonal of $\tilde\omega(k)$.

<!-- eq:4.5 -->
$$ T(k) \triangleq \nabla_{\theta(k)} \mathcal{P}[\theta] $$
- **what:** generalized (hinge) force from the internal-coordinate potential.
- **symbols:** $T(k)$ - $m$-vector generalized force at hinge $k$; $\nabla_{\theta(k)}$ - gradient w.r.t. hinge coordinates.

<!-- eq:4.6 -->
$$ \hat{f}_c(k) \triangleq \begin{pmatrix} \sum_{i=1}^{r(k)} \tilde{l}\left[ \mathcal{O}_k^-, x_i(k) \right] \hat{f}_i(k) \\ \sum_{i=1}^{r(k)} \hat{f}_i(k) \end{pmatrix} $$
- **what:** reduce per-atom Cartesian forces of a cluster to one effective 6-dim spatial force at the hinge frame (torque part on top, net force on bottom). This is the per-body spatial force Robosample builds from OpenMM per-atom forces.
- **symbols:** $\hat f_c(k)$ - 6-dim effective spatial Cartesian force on cluster $k$ about $\mathcal{O}_k^-$; $r(k)$ - number of atoms in cluster $k$; $\hat f_i(k) = \nabla_{x_i(k)}\mathcal{P}[x]$ - 3-dim Cartesian force on atom $i$; $\tilde{l}[\mathcal{O}_k^-, x_i(k)]$ - skew tensor of the vector from hinge frame $\mathcal{O}_k^-$ to atom $i$.

<!-- eq:4.7 -->
$$ f(k) = \phi(k, k-1)\, f(k-1) + M(k)\, \alpha(k) + b(k) + \hat{f}_c(k) $$
$$ T(k) = H(k)\, f(k) $$
- **what:** tip-to-base force recursion and hinge force projection. Incorporating $T(k)$ (eq:4.5) and $\hat f_c(k)$ (eq:4.6) separately avoids computing $\nabla_{\theta(k)}\mathcal{P}[x]$.
- **symbols:** $f(k)$ - 6-dim spatial interaction force at $\mathcal{O}_k^-$ between clusters $k{+}1$ and $k$; $\phi(k,k{-}1)$ - spatial transform from child; $M(k)$ - $6\times 6$ spatial inertia about $\mathcal{O}_k^-$ (eq:A.4); $b(k)$ - gyroscopic spatial force (eq:A.6); $H(k)$ - $m\times 6$ (transpose of $H^*(k)$).

<!-- eq:4.8 -->
$$ \begin{aligned}
&V(n+1) = 0,\quad \alpha(n+1) = 0 \\
&\textbf{for } k = n \cdots 1: \\
&\quad V(k) = \phi^*(k+1,k)\, V(k+1) + H^*(k)\, \dot{\theta}(k) \\
&\quad \alpha(k) = \phi^*(k+1,k)\, \alpha(k+1) + H^*(k)\, \ddot{\theta}(k) + a(k) \\
&f(0) = 0 \\
&\textbf{for } k = 1 \cdots n: \\
&\quad f(k) = \phi(k,k-1)\, f(k-1) + M(k)\, \alpha(k) + b(k) + \hat{f}_c(k) \\
&\quad T(k) = H(k)\, f(k)
\end{aligned} $$
- **what:** full Newton-Euler recursive EOM (inverse dynamics): base-to-tip sweep for $V,\alpha$, then tip-to-base sweep for $f,T$. <!-- CHECK: index ranges read from prose; base-to-tip means k=n..1, tip-to-base means k=1..n given tip=1, base=n numbering. -->
- **symbols:** as above; loops run over all clusters.

<!-- eq:4.9 -->
$$ \mathscr{E}_{\phi} \triangleq \begin{pmatrix} 0 & 0 & \cdots & 0 & 0 \\ \phi(2,1) & 0 & \cdots & 0 & 0 \\ 0 & \phi(3,2) & \cdots & 0 & 0 \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \cdots & \phi(n,n-1) & 0 \end{pmatrix} \in \Re^{6n \times 6n} $$
- **what:** shift operator holding first-subdiagonal spatial transforms; strictly lower block-bidiagonal.
- **symbols:** $\phi(k{+}1,k)$ - $6\times 6$ spatial transform (eq:A.2); zero blocks are $6\times 6$.

<!-- eq:4.10 -->
$$ V = \mathcal{E}_{\phi}^{*} V + H^{*}\dot{\theta} $$
- **what:** stacked operator form of the velocity recursion.
- **symbols:** $V=[V^*(1)\cdots V^*(n)]^*$ - $6n$-stacked spatial velocity; $H=\mathrm{diag}\{H(k)\}$; $\dot\theta$ - stacked hinge velocities.

<!-- eq:4.11 -->
$$ \phi \triangleq [I - \mathcal{E}_{\phi}]^{-1} = \begin{pmatrix} I & 0 & \cdots & 0 \\ \phi(2,1) & I & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ \phi(n,1) & \phi(n,2) & \cdots & I \end{pmatrix} \in \Re^{6n \times 6n}, \qquad \phi(i,j) \triangleq \phi(i,i-1)\cdots\phi(j+1,j)\ \text{for } i>j $$
- **what:** lower-triangular spatial propagation operator (inverse of $I-\mathcal{E}_\phi$); composite transforms have the semigroup property.
- **symbols:** $I$ - $6\times 6$ identity blocks; $\phi(i,j)$ - composite transform from frame $j$ to frame $i$.

<!-- eq:4.12 -->
$$ V = \phi^* H^* \dot{\theta} $$
- **what:** closed operator form of stacked spatial velocities.
- **symbols:** as above.

<!-- eq:4.13 -->
$$ \begin{aligned}
V &= \phi^* H^* \dot{\theta} \\
\alpha &= \phi^* (H^* \ddot{\theta} + a) \\
f &= \phi (M\alpha + b + \hat{f}_c) = \phi M \phi^* H^* \ddot{\theta} + \phi (M \phi^* a + b + \hat{f}_c) \\
T &= H f = H \phi M \phi^* H^* \ddot{\theta} + H \phi (M \phi^* a + b + \hat{f}_c)
\end{aligned} $$
- **what:** operator-level EOM assembling velocity, acceleration, force, and generalized force.
- **symbols:** $M=\mathrm{diag}\{M(k)\}$; $a,b,\hat f_c$ - stacked Coriolis, gyroscopic, Cartesian spatial forces.

<!-- eq:4.14 -->
$$ T = \mathcal{M}(\theta)\ddot{\theta} + \mathcal{C}(\theta, \dot{\theta}) $$
- **what:** operator EOM (same as eq:1.1, sign of $\mathcal{C}$ moved).
- **symbols:** as eq:1.1.

<!-- eq:4.15a -->
$$ \mathcal{M}(\theta) \triangleq H\phi M\phi^* H^* \in \Re^{\mathcal{N} \times \mathcal{N}} $$
- **what:** Newton-Euler operator factorization of the mass matrix (equivalent to eq:4.8 recursion).
- **symbols:** $H,\phi,M$ - operators above.

<!-- eq:4.15b -->
$$ \mathcal{C}(\theta, \dot{\theta}) \triangleq H\phi(M\phi^*a + b + \hat{f}_c) \in \Re^{\mathcal{N}} $$
- **what:** Coriolis/centrifugal/gyroscopic/Cartesian force vector in operator form.
- **symbols:** as above.

<!-- eq:5.1 -->
$$ \begin{aligned}
&P^{+}(0) = 0 \in \Re^{6 \times 6} \\
&\textbf{for } k = 1 \cdots n: \\
&\quad P(k) = \phi(k, k-1)\, P^{+}(k-1)\, \phi^{*}(k, k-1) + M(k) \\
&\quad D(k) = H(k)\, P(k)\, H^{*}(k) \\
&\quad G(k) = P(k)\, H^{*}(k)\, D^{-1}(k) \\
&\quad K(k+1, k) = \phi(k+1, k)\, G(k) \\
&\quad \bar{\tau}(k) = I - G(k)\, H(k) \\
&\quad P^{+}(k) = \bar{\tau}(k)\, P(k) \\
&\quad \psi(k+1, k) = \phi(k+1, k)\, \bar{\tau}(k)
\end{aligned} $$
- **what:** tip-to-base discrete Riccati recursion defining the articulated-body quantities. This is the core of the O(N) solver.
- **symbols:** $P(k)$ - $6\times 6$ articulated-body spatial inertia (Riccati); $D(k)$ - $m\times m$ hinge inertia $H P H^*$; $G(k)$ - $6\times m$ Kalman gain $PH^*D^{-1}$; $K(k{+}1,k)$ - $6\times m$ shifted gain; $\bar\tau(k)$ - $6\times 6$ projector $I-GH$; $P^+(k)$ - $6\times 6$ propagated inertia; $\psi(k{+}1,k)$ - $6\times 6$ articulated-body transform.

<!-- eq:5.2 -->
$$ \begin{aligned}
D &\triangleq HPH^* = \operatorname{diag}\{D(k)\} \in \Re^{\mathcal{N} \times \mathcal{N}} \\
G &\triangleq PH^*D^{-1} = \operatorname{diag}\{G(k)\} \in \Re^{6n \times \mathcal{N}} \\
K &\triangleq \mathscr{E}_{\phi}G \in \Re^{6n \times \mathcal{N}} \\
\bar{\tau} &\triangleq I - GH = \operatorname{diag}\{\bar{\tau}(k)\} \in \Re^{6n \times 6n}
\end{aligned} $$
- **what:** operator-level definitions of the Riccati quantities. Only nonzero blocks of $K$ are $K(k{+}1,k)$ on the first subdiagonal.
- **symbols:** $P=\mathrm{diag}\{P(k)\}$; others as eq:5.1.

<!-- eq:5.3 -->
$$ \mathscr{E}_{\psi} \triangleq \mathscr{E}_{\phi} \bar{\tau} = \begin{pmatrix} 0 & 0 & \cdots & 0 & 0 \\ \psi(2,1) & 0 & \cdots & 0 & 0 \\ 0 & \psi(3,2) & \cdots & 0 & 0 \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \cdots & \psi(n,n-1) & 0 \end{pmatrix} \in \Re^{6n \times 6n} $$
- **what:** shift operator built from articulated-body transforms $\psi(k{+}1,k)$.
- **symbols:** as eq:5.1.

<!-- eq:5.4 -->
$$ \psi \triangleq [I - \mathscr{E}_{\psi}]^{-1} = \begin{pmatrix} I & 0 & \cdots & 0 \\ \psi(2,1) & I & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ \psi(n,1) & \psi(n,2) & \cdots & I \end{pmatrix} \in \Re^{6n \times 6n} $$
- **what:** lower-triangular articulated-body propagation operator (analog of $\phi$).
- **symbols:** as eq:5.1.

<!-- eq:5.5 -->
$$ \psi(i,j) \triangleq \psi(i,i-1)\cdots\psi(j+1,j) \quad \text{for } i>j $$
- **what:** semigroup composition of $\psi$ transforms.
- **symbols:** as above.

<!-- eq:5.6 -->
$$ \mathcal{M} = [I + H\phi K]\, D\, [I + H\phi K]^* $$
- **what:** Lemma 5.1 — innovations operator factorization; closed-form block $LDL^*$ decomposition of $\mathcal{M}$. Factor $[I+H\phi K]$ is square, block lower-triangular, nonsingular; $D$ block-diagonal.
- **symbols:** as above.

<!-- eq:5.7 -->
$$ [I + H\phi K]^{-1} = [I - H\psi K] $$
- **what:** Lemma 5.2 — closed-form inverse of the triangular factor.
- **symbols:** as above.

<!-- eq:5.8 -->
$$ \mathcal{M}^{-1} = [I - H\psi K]^* D^{-1} [I - H\psi K] $$
- **what:** Lemma 5.3 — closed-form mass-matrix inverse ($LDL^*$ of $\mathcal{M}^{-1}$).
- **symbols:** as above.

<!-- eq:5.9 -->
$$ \ddot{\theta} = [I - H\psi K]^* D^{-1} \left[ T - H\psi \{KT + Pa + b + \hat{f}_c\} \right] - K^*\psi^* a $$
- **what:** Lemma 5.4 — closed operator expression for generalized accelerations from hinge forces $T$ and Cartesian spatial forces $\hat f_c$.
- **symbols:** as above; $a,b$ - stacked Coriolis/gyroscopic terms.

<!-- eq:5.10 -->
$$ \begin{aligned}
z &= \psi [KT + Pa + b + \hat{f}_c] \\
\varepsilon &= T - Hz \\
v &= D^{-1}\varepsilon \\
\alpha &= \psi [H^*v + a] \\
\ddot{\theta} &= v - K^*\alpha
\end{aligned} $$
- **what:** decomposition of eq:5.9 into a computable operator sequence.
- **symbols:** $z$ - $6n$ residual force; $\varepsilon$ - $\mathcal{N}$ innovation/residual; $v$ - $\mathcal{N}$ scaled residual; $\alpha$ - $6n$ spatial acceleration.

<!-- eq:5.11a -->
$$ \begin{aligned}
&z^{+}(0) = 0 \\
&\textbf{for } k = 1 \cdots n: \\
&\quad z(k) = \phi(k, k-1)\, z^{+}(k-1) + P(k)\, a(k) + b(k) + \hat{f}_c(k) \\
&\quad \varepsilon(k) = T(k) - H(k)\, z(k) \\
&\quad v(k) = D^{-1}(k)\, \varepsilon(k) \\
&\quad z^{+}(k) = z(k) + G(k)\, \varepsilon(k)
\end{aligned} $$
- **what:** O(N) tip-to-base recursion for residual forces and innovations (mergeable with eq:5.1).
- **symbols:** as eq:5.1/eq:5.10, per-cluster.

<!-- eq:5.11b -->
$$ \begin{aligned}
&\alpha(n+1) = 0 \\
&\textbf{for } k = n \cdots 1: \\
&\quad \alpha^{+}(k) = \phi^{*}(k+1, k)\, \alpha(k+1) \\
&\quad \ddot{\theta}(k) = v(k) - G^{*}(k)\, \alpha^{+}(k) \\
&\quad \alpha(k) = \alpha^{+}(k) + H^{*}(k)\, \ddot{\theta}(k) + a(k)
\end{aligned} $$
- **what:** O(N) base-to-tip recursion computing the generalized accelerations $\ddot\theta(k)$.
- **symbols:** $\alpha^+(k)$ - propagated spatial acceleration before adding hinge contribution; others as above. Note $K^*(k{+}1,k)=G^*(k)\phi^*(k{+}1,k)$ so $G^*(k)\alpha^+(k)$ realizes the $K^*\alpha$ term.

<!-- eq:A.1 -->
$$ V(\mathcal{O}) \triangleq \begin{pmatrix} \omega(\mathcal{O}) \\ v(\mathcal{O}) \end{pmatrix}, \qquad f(\mathcal{O}) \triangleq \begin{pmatrix} N(\mathcal{O}) \\ F(\mathcal{O}) \end{pmatrix} $$
- **what:** spatial velocity (angular over linear) and spatial force (moment over force).
- **symbols:** $\omega$ - angular velocity; $v$ - linear velocity; $N$ - moment/torque; $F$ - force; all 3-dim.

<!-- eq:A.2 -->
$$ \phi(\mathcal{O}_x, \mathcal{O}_y) \triangleq \begin{pmatrix} I_3 & \tilde{l}(\mathcal{O}_x, \mathcal{O}_y) \\ 0_3 & I_3 \end{pmatrix} $$
- **what:** $6\times 6$ spatial (rigid) transformation matrix between two frames.
- **symbols:** $I_3,0_3$ - $3\times 3$ identity/zero; $\tilde{l}(\mathcal{O}_x,\mathcal{O}_y)$ - skew tensor of the vector from origin of $\mathcal{O}_x$ to origin of $\mathcal{O}_y$.

<!-- eq:A.3 -->
$$ f(\mathcal{O}_x) = \phi(\mathcal{O}_x, \mathcal{O}_y)\, f(\mathcal{O}_y), \qquad V(\mathcal{O}_y) = \phi^*(\mathcal{O}_x, \mathcal{O}_y)\, V(\mathcal{O}_x) $$
- **what:** how the spatial transform maps forces (forward) and velocities (transpose) between frames.
- **symbols:** `*` - transpose.

<!-- eq:A.4 -->
$$ M(\mathcal{O}) \triangleq \begin{pmatrix} \mathcal{I}(\mathcal{O}) & m\tilde{p} \\ -m\tilde{p} & mI_3 \end{pmatrix} $$
- **what:** $6\times 6$ spatial inertia of a rigid body about frame $\mathcal{O}$.
- **symbols:** $m$ - body mass; $p$ - vector from $\mathcal{O}$ to center of mass; $\tilde p$ - its skew tensor; $\mathcal{I}(\mathcal{O})$ - $3\times 3$ inertia matrix about $\mathcal{O}$.

<!-- eq:A.5 -->
$$ f(\mathcal{O}) = M(\mathcal{O})\, \alpha(\mathcal{O}) + b(\mathcal{O}) $$
- **what:** rigid-body spatial equations of motion about frame $\mathcal{O}$.
- **symbols:** $f$ - effective spatial force; $\alpha$ - spatial acceleration; $b$ - gyroscopic spatial force.

<!-- eq:A.6 -->
$$ b(\mathcal{O}) = \begin{pmatrix} \tilde{\omega}(\mathcal{O})\, \mathcal{I}(\mathcal{O})\, \omega(\mathcal{O}) \\ m\tilde{\omega}(\mathcal{O})\, \tilde{\omega}(\mathcal{O})\, l(\mathcal{O}, \mathcal{O}_{CM}) \end{pmatrix} $$
- **what:** gyroscopic (velocity-dependent) spatial force of a rigid body.
- **symbols:** $\tilde\omega$ - skew of angular velocity; $\mathcal{I}$ - inertia; $l(\mathcal{O},\mathcal{O}_{CM})$ - vector from $\mathcal{O}$ to center of mass; $m$ - mass.

## Derivations (not implemented)

Appendix B proofs of Lemmas 5.1, 5.2, 5.4 are pure operator algebra; key
identities are eq:B.1 ($\bar\tau P\bar\tau^*=\bar\tau P$ / Riccati rewrite),
eq:B.2 ($[I+H\phi K]^{-1}=I-H\phi[I+KH\phi]^{-1}K$), eq:B.3
($\psi^{-1}=\phi^{-1}+KH$), eq:B.5 ($[I-H\psi K]H\phi=H\psi$), eq:B.7 ($\psi
M\phi^*=\psi P+P\tilde\phi^*$). Not required for a port; see paper.md.
