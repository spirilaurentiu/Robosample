# Equations — Jain 1997, Compensating Mass Matrix Potential

Implementable equations for computing the Fixman compensating potential
$\mathcal{V}_c$ and its gradient (compensating torque $T_c$) via $O(\mathcal{N})$
spatial-operator / articulated-body-inertia recursions.

Pure proof-algebra steps (the trace manipulations proving Lemmas 4.2–4.5) are
kept in `paper.md` under `Derivation (not implemented)` markers and are not
reproduced here.

---

<!-- eq:2.3 -->
$$ \text{K.E.} = \tfrac{1}{2}\beta^* \mathcal{M}(\boldsymbol{\theta})\beta = \tfrac{1}{2}p^* \mathcal{M}^{-1}(\boldsymbol{\theta})p $$
- **what:** Kinetic energy of the constrained (internal-coordinate) system in velocity and momentum form.
- **symbols:** $\beta$ - internal velocity coordinates vector ($\mathbb{R}^{\mathcal{N}}$); $\mathcal{M}(\boldsymbol{\theta})$ - configuration-dependent mass matrix / metric tensor ($\mathbb{R}^{\mathcal{N}\times\mathcal{N}}$); $p$ - conjugate momenta ($\mathbb{R}^{\mathcal{N}}$); $\boldsymbol{\theta}$ - internal (hinge) coordinates.

<!-- eq:2.4 -->
$$ p = \mathcal{M}(\boldsymbol{\theta})\beta $$
- **what:** Conjugate momentum from generalized velocity.
- **symbols:** as above.

<!-- eq:2.9 -->
$$ \mathcal{V}'(\boldsymbol{\theta}) \triangleq \mathcal{V}(\boldsymbol{\theta}) + \mathcal{V}_c(\boldsymbol{\theta}), \qquad \mathcal{V}_c(\boldsymbol{\theta}) \triangleq \tfrac{1}{2} \ln \det \{ \mathcal{M}(\boldsymbol{\theta}) \} $$
- **what:** Fixman modified potential = standard potential + compensating mass-matrix (metric-tensor) potential. Adding $\mathcal{V}_c$ removes the $\det\{\mathcal{M}^{1/2}\}$ bias from constrained-MD ensemble averages.
- **symbols:** $\mathcal{V}$ - standard potential energy; $\mathcal{V}'$ - modified potential used in constrained MD; $\mathcal{V}_c$ - compensating potential; $\ln\det$ - natural log of determinant.

<!-- eq:2.10 -->
$$ T' = \nabla_{\theta} \mathcal{V}'(\boldsymbol{\theta}) = T + T_c, \qquad T \triangleq \nabla_{\theta} \mathcal{V}(\boldsymbol{\theta}), \quad T_c \triangleq \nabla_{\theta} \mathcal{V}_c(\boldsymbol{\theta}) $$
- **what:** Total generalized (hinge) torque used for forces = standard torque + compensating torque.
- **symbols:** $T'$ - total hinge torque vector; $T$ - standard-potential torque; $T_c$ - compensating torque; $\nabla_\theta$ - gradient w.r.t. internal coordinates.

<!-- eq:2.11 -->
$$ T_c(k) = \frac{\partial \mathcal{V}_c(\boldsymbol{\theta})}{\partial \boldsymbol{\theta}(k)} = \frac{1}{2} \frac{\partial \ln \det\{\mathcal{M}(\boldsymbol{\theta})\}}{\partial \boldsymbol{\theta}(k)} $$
- **what:** $k$th component of the compensating torque as a derivative of $\mathcal{V}_c$.
- **symbols:** $T_c(k)$ - compensating torque on hinge $k$; $\boldsymbol{\theta}(k)$ - $k$th internal coordinate.

<!-- eq:2.14 -->
$$ T_c(k) = \frac{1}{2}\operatorname{Trace}\!\left\{\mathcal{M}^{-1}(\boldsymbol{\theta})\frac{\partial\mathcal{M}(\boldsymbol{\theta})}{\partial\boldsymbol{\theta}(k)}\right\} = \frac{1}{2}\operatorname{Trace}\!\left\{\mathcal{M}^{-1}(\boldsymbol{\theta})\,\mathcal{M}_{\boldsymbol{\theta}(k)}(\boldsymbol{\theta})\right\} $$
- **what:** Compensating torque as a trace of $\mathcal{M}^{-1}$ times the mass-matrix sensitivity (uses $\partial\ln\det X/\partial X = (X^*)^{-1}$). The naive $O(\mathcal{N}^3)$ baseline; Lemmas 4.2–4.4 reduce it to $O(\mathcal{N})$.
- **symbols:** $\mathcal{M}_{\theta(k)} \equiv \partial\mathcal{M}/\partial\theta(k)$ - mass-matrix sensitivity; $\operatorname{Trace}$ - matrix trace.

<!-- eq:3.4a -->
$$ \mathcal{M}(\boldsymbol{\theta}) \triangleq H \phi M \phi^* H^* \in \mathbb{R}^{\mathcal{N} \times \mathcal{N}} $$
- **what:** Newton–Euler operator factorization of the system mass matrix.
- **symbols:** $H$ - joint-map (hinge axis) operator ($\mathcal{N}\times 6n$); $\phi$ - rigid-body transformation propagation operator ($6n\times 6n$, lower triangular); $M$ - block-diagonal spatial (link) inertia operator ($6n\times 6n$, blocks $M(k)\in\mathbb{R}^{6\times6}$); $*$ - transpose (adjoint).

<!-- eq:3.5a -->
$$ \mathcal{M} = [I + H\phi K]\,D\,[I + H\phi K]^* $$
- **what:** Innovations (block $LDL^*$) factorization of the mass matrix. Factor $[I+H\phi K]$ is square, block lower triangular, nonsingular, with identity blocks on the diagonal; $D$ block diagonal.
- **symbols:** $K$ - Kalman gain operator ($6n\times \mathcal{N}$); $D$ - block-diagonal hinge inertia ($\mathcal{N}\times\mathcal{N}$); $I$ - identity.

<!-- eq:3.5b -->
$$ [I + H\phi K]^{-1} = [I - H\psi K] $$
- **what:** Closed-form inverse of the innovations factor.
- **symbols:** $\psi$ - articulated-body transformation propagation operator (block lower triangular).

<!-- eq:3.5c -->
$$ \mathcal{M}^{-1} = [I - H\psi K]^*\, D^{-1}\,[I - H\psi K] $$
- **what:** Closed-form block $LDL^*$ decomposition of $\mathcal{M}^{-1}$.
- **symbols:** as above.

<!-- eq:3.6 -->
$$
\begin{aligned}
&P^{+}(0) = 0 \\
&\textbf{for } k = 1 \cdots \mathcal{N}: \\
&\quad P(k) = \phi(k, k-1)\,P^{+}(k-1)\,\phi^{*}(k, k-1) + M(k) \\
&\quad D(k) = H(k)\,P(k)\,H^{*}(k) \\
&\quad G(k) = P(k)\,H^{*}(k)\,D^{-1}(k) \\
&\quad K(k+1, k) = \phi(k+1, k)\,G(k) \\
&\quad \tau(k) = G(k)\,H(k) \\
&\quad \bar{\tau}(k) = I - \tau(k) \\
&\quad P^{+}(k) = \bar{\tau}(k)\,P(k) \\
&\quad \psi(k+1, k) = \phi(k+1, k)\,\bar{\tau}(k)
\end{aligned}
$$
- **what:** Algorithm 3.1 — tip-to-base articulated-body-inertia (Riccati / Kalman) recursion. $P(k)$ is the articulated body inertia outboard of hinge $k$. Supplies $P(k), D(k), G(k), K, \bar\tau, \psi$ used by the acceleration solve and by the compensating-torque computation.
- **symbols:** $P(k)\in\mathbb{R}^{6\times6}$ - articulated body inertia; $P^+(k)$ - inertia after projecting out hinge $k$; $D(k)\in\mathbb{R}^{d_k\times d_k}$ - hinge inertia ($d_k$ = dof of hinge $k$, =1 for single-dof rotational hinge); $G(k)$ - Kalman gain block; $\tau,\bar\tau$ - projection operators; $\phi(k,k-1)$ - rigid transform between adjacent links.

<!-- eq:3.8 -->
$$ \ddot{\boldsymbol{\theta}} = [I - H\psi K]^*\,D^{-1}\,[T - H\psi\{KT + Pa + b + \hat{f}_c\}] - K^*\psi^* a $$
- **what:** Lemma 3.2 — closed-form generalized accelerations from hinge forces $T$ and Cartesian spatial forces $\hat f_c$.
- **symbols:** $\ddot{\boldsymbol{\theta}}$ - generalized accelerations; $a$ - Coriolis/velocity-dependent acceleration terms; $b$ - gyroscopic spatial-force terms; $\hat f_c$ - applied Cartesian spatial forces.

<!-- eq:3.9 -->
$$
\begin{aligned}
&z^{+}(0) = 0 \\
&\textbf{for } k = 1 \cdots n: \\
&\quad z(k) = \phi(k, k-1)\,z^{+}(k-1) + P(k)\,a(k) + b(k) + \hat{f}_{c}(k) \\
&\quad \varepsilon(k) = T(k) - H(k)\,z(k) \\
&\quad \nu(k) = D^{-1}(k)\,\varepsilon(k) \\
&\quad z^{+}(k) = z(k) + G(k)\,\varepsilon(k)
\end{aligned}
$$
- **what:** Tip-to-base residual-force sweep of the $O(\mathcal{N})$ acceleration solve. Can be fused with the Eq. (3.6) sweep.
- **symbols:** $z(k), z^+(k)$ - residual spatial forces; $\varepsilon(k)$ - innovation (residual hinge force); $\nu(k)$ - $D^{-1}$-weighted innovation.
<!-- CHECK: raw OCR had "\nu(k = D^{-1}(k)\varepsilon(k)" — repaired to \nu(k)=D^{-1}(k)\varepsilon(k) -->

<!-- eq:3.10 -->
$$
\begin{aligned}
&\alpha(n+1) = 0 \\
&\textbf{for } k = n \cdots 1: \\
&\quad \alpha^{+}(k) = \phi^{*}(k+1, k)\,\alpha(k+1) \\
&\quad \ddot{\boldsymbol{\theta}}(k) = \nu(k) - G^{*}(k)\,\alpha^{+}(k) \\
&\quad \alpha(k) = \alpha^{+}(k) + H^{*}(k)\,\ddot{\boldsymbol{\theta}}(k) + a(k)
\end{aligned}
$$
- **what:** Base-to-tip acceleration sweep completing the $O(\mathcal{N})$ solve for $\ddot\theta(k)$.
- **symbols:** $\alpha(k)$ - spatial acceleration of link $k$; $\alpha^+(k)$ - propagated parent acceleration.

<!-- eq:4.11 -->
$$ \mathcal{M}_{\boldsymbol{\theta}_{i}} = H \phi \left[\mathbb{H}_{\delta}^{i}\,\phi M - M \phi^{*}\,\mathbb{H}_{\delta}^{i}\right] \phi^{*} H^{*} $$
- **what:** Lemma 4.1 — closed-form spatial-operator expression for the mass-matrix sensitivity w.r.t. hinge angle $\theta_i$.
- **symbols:** $\mathbb{H}_\delta^i\in\mathbb{R}^{6n\times6n}$ - all-zero except a single $6\times6$ block $\mathbb{H}(i)$ at the $i$th diagonal location; other operators as in Eq. (3.4a).

<!-- eq:4.12 -->
$$ \mathbb{H}(i) = \begin{pmatrix} \tilde{h}(i) & 0 \\ 0 & \tilde{h}(i) \end{pmatrix}, \qquad \tilde{v} \triangleq \begin{pmatrix} 0 & -z & y \\ z & 0 & -x \\ -y & x & 0 \end{pmatrix} \text{ for } v=\begin{pmatrix}x\\y\\z\end{pmatrix} $$
- **what:** The $6\times6$ diagonal block from the hinge rotational-axis unit vector $h(i)$; $\tilde v$ is the $3\times3$ skew (cross-product) matrix of a 3-vector.
- **symbols:** $h(i)\in\mathbb{R}^3$ - hinge $i$ rotational axis unit vector; $\tilde h(i)\in\mathbb{R}^{3\times3}$ - its skew matrix.

<!-- eq:4.13 -->
$$ T_{c}(i) = \operatorname{Trace}\{P\,\Omega\,\mathbb{H}_{\delta}^{i}\}, \qquad \Omega \triangleq \psi^{*} H^{*} D^{-1} H \psi \in \mathbb{R}^{6n \times 6n} $$
- **what:** Lemma 4.2 — compensating torque as a trace over articulated-body operators. Note factor of $\tfrac12$ from Eq. (2.14) is absorbed here (the $\mathcal{M}_{\theta_i}$ of Eq. 4.11 is anti-symmetrized, yielding $2\times$ the symmetric contribution).
- **symbols:** $\Omega$ - articulated-body operator; $P$ - block-diagonal articulated body inertia operator.

<!-- eq:4.15 -->
$$ \Omega = Y + \tilde{\psi}^{*} Y + Y \tilde{\psi}, \qquad \tilde{\psi} \triangleq \psi - I $$
- **what:** Decomposition of $\Omega$ into a block-diagonal $Y$ and strictly-lower/upper triangular corrections.
- **symbols:** $Y\in\mathbb{R}^{6n\times6n}$ - block diagonal, blocks $Y(k,k)\in\mathbb{R}^{6\times6}$; $\tilde\psi=\psi-I$ - strictly lower triangular.

<!-- eq:4.16 -->
$$
\begin{aligned}
&Y(n+1) = 0 \\
&\textbf{for } k = n \cdots 1: \\
&\quad Y(k) = \psi^{*}(k+1, k)\,Y(k+1)\,\psi(k+1, k) + H^{*}(k)\,D^{-1}(k)\,H(k)
\end{aligned}
$$
- **what:** Base-to-tip recursion for the block-diagonal $Y(k)$ needed by the compensating-torque formula (Eq. 4.17/4.18).
- **symbols:** $Y(k)\in\mathbb{R}^{6\times6}$ - block-diagonal element; $\psi(k+1,k)$ - articulated-body transform (from Eq. 3.6).
<!-- CHECK: index label — text says the loop runs k=n..1 (base-to-tip); Y(n+1)=0 is the tip boundary condition -->

<!-- eq:4.17 -->
$$ T_c(i) = \operatorname{Trace}\{P(i)\,Y(i)\,\mathbb{H}(i)\} $$
- **what:** Lemma 4.3 — compensating torque reduced to a single $6\times6$ block trace at hinge $i$.
- **symbols:** $P(i), Y(i), \mathbb{H}(i)\in\mathbb{R}^{6\times6}$.

<!-- eq:4.18 -->
$$ T_c(i) = -h^{*}(i)\,\mathcal{F}[Q_{11} + Q_{22}], \qquad P(i)Y(i) = \begin{pmatrix} Q_{11} & Q_{12} \\ Q_{21} & Q_{22} \end{pmatrix} $$
- **what:** Lemma 4.4 — final closed-form scalar compensating torque per single-dof rotational hinge. Partition $P(i)Y(i)$ into $3\times3$ blocks; take the top-left + bottom-right, apply the axial map $\mathcal{F}$, dot with the (negated) hinge axis.
- **symbols:** $Q_{jk}\in\mathbb{R}^{3\times3}$ - blocks of $P(i)Y(i)$; $h(i)\in\mathbb{R}^3$ - hinge axis; $\mathcal{F}[\cdot]:\mathbb{R}^{3\times3}\to\mathbb{R}^3$ - axial-vector map defined below.

<!-- eq:4.18b -->
$$ v = \mathcal{F}[A] \iff \tilde{v} = A - A^{*}, \qquad \operatorname{Trace}\{A\tilde{v}\} = -v^{*}\mathcal{F}[A] $$
- **what:** Definition of the axial-vector map $\mathcal{F}$ (inverse-skew of the antisymmetric part) and the trace identity used to derive Eq. (4.18).
- **symbols:** $A\in\mathbb{R}^{3\times3}$; $v\in\mathbb{R}^3$; $\tilde v$ - skew matrix of $v$.

<!-- eq:4.20 -->
$$ \mathcal{V}_c(\boldsymbol{\theta}) = \frac{1}{2} \sum_{i=1}^{n} \ln \det \{ D(i) \} $$
- **what:** Lemma 4.5 — closed-form compensating potential as a sum of log-dets of the (small, block-diagonal) hinge inertias $D(i)$, since $\det\{\mathcal{M}\}=\det\{D\}=\prod_i\det\{D(i)\}$. This is the cheap way to evaluate $\mathcal{V}_c$ from Algorithm 3.1.
- **symbols:** $D(i)\in\mathbb{R}^{d_i\times d_i}$ - hinge inertia block from Eq. (3.6); $n$ - number of hinges/clusters.
