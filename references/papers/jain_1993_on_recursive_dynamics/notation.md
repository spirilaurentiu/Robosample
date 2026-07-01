# Notation — Jain, Vaidehi, Rodriguez 1993

## Conventions

- **Cluster numbering:** clusters $1..n$ ordered **tip (1) to base (n)**; the
  inertial frame has index $n{+}1$. Parent of cluster $k$ is $k{+}1$ (toward
  base); child is $k{-}1$. `f(0)=0` (no force inboard of the tip).
- **Sweeps:** "base-to-tip" = decreasing $k$ ($n\to 1$); "tip-to-base" =
  increasing $k$ ($1\to n$). (Because tip=1, base=n.)
- **Spatial quantities are 6-dimensional**, stacked **[angular; linear]** for
  velocity/acceleration and **[moment/torque; force]** for forces. Consistent
  with eq:4.6: torque block on top, net force block on bottom.
- **`*` denotes matrix transpose** (adjoint of a spatial operator), NOT complex
  conjugate.
- **$\tilde{x}$** is the $3\times 3$ skew-symmetric cross-product tensor of a
  3-vector $x$ ($\tilde x\, y = x\times y$).
- **$H^*$ is $6\times m$**, $H$ is $m\times 6$; hinge matrix components are
  **constant** in the $\mathcal{O}_k^\pm$ frames (enforced by pseudocluster
  decomposition for multi-dof hinges).
- Kinematic map $\dot\theta=\Pi(\theta)\beta$ assumed identity ($\beta=\dot\theta$)
  throughout, WLOG.
- Units: not fixed to reduced units — standard MD units; masses/inertias physical.

## Symbols

| symbol | meaning | shape / dtype | convention |
|---|---|---|---|
| $\mathcal{N}$ | number of degrees of freedom | scalar int | $=3n-m$ |
| $n$ | number of atoms (Cartesian) / clusters (internal) | scalar int | context-dependent |
| $m$ | number of constraints; or hinge dof | scalar int | |
| $\theta$ | generalized (internal) coordinates | $\mathcal{N}$-vec | per hinge $\theta(k)\in\Re^m$ |
| $\dot\theta,\ddot\theta$ | generalized velocity, acceleration | $\mathcal{N}$-vec | |
| $\mathcal{M}(\theta)$ | system mass matrix | $\mathcal{N}\times\mathcal{N}$ | symmetric PD, non-diagonal |
| $\mathcal{C}(\theta,\dot\theta)$ | Coriolis/velocity-dependent forces | $\mathcal{N}$-vec | |
| $T$ | generalized (hinge) forces | $\mathcal{N}$-vec; $T(k)\in\Re^m$ | $T(k)=\nabla_{\theta(k)}\mathcal{P}[\theta]$ |
| $\mathcal{P}[\theta],\mathcal{P}[x]$ | internal / Cartesian potential | scalar | $PE=\mathcal{P}[\theta]+\mathcal{P}[x]$ |
| $x$ | Cartesian atom positions | $\Re^{3n}$ | |
| $\mathcal{M}_c$ | Cartesian mass matrix | $3n\times 3n$ diagonal | atomic masses |
| $f$ | inter-atomic Cartesian forces | $\Re^{3n}$ | |
| $A,B$ | constraint matrix / vector | $\Re^{m\times 3n}$, $\Re^m$ | $A\dot x=B$ |
| $V(k)$ | spatial velocity at $\mathcal{O}_k^-$ | 6-vec [ω; v] | |
| $\omega,v$ | angular, linear velocity | 3-vec each | |
| $\alpha(k)$ | spatial acceleration $=\dot V$ | 6-vec | |
| $a(k)$ | Coriolis spatial acceleration | 6-vec | eq:4.4 |
| $f(k)$ | spatial interaction force at $\mathcal{O}_k^-$ | 6-vec [N; F] | between clusters $k{+}1,k$ |
| $N,F$ | moment/torque, force | 3-vec each | |
| $b(k)$ | gyroscopic spatial force | 6-vec | eq:A.6 |
| $\hat f_i(k)$ | Cartesian force on atom $i$ of cluster $k$ | 3-vec | $=\nabla_{x_i(k)}\mathcal{P}[x]$ |
| $\hat f_c(k)$ | effective spatial Cartesian force on cluster | 6-vec | eq:4.6 |
| $r(k)$ | atoms in cluster $k$ | scalar int | |
| $H^*(k)$ | hinge matrix | $6\times m$ | constant in $\mathcal{O}_k^\pm$ |
| $M(k)$ | spatial inertia of cluster $k$ about $\mathcal{O}_k^-$ | $6\times 6$ | eq:A.4 |
| $\mathcal{I}(\mathcal{O})$ | $3\times 3$ inertia matrix about frame | $3\times 3$ | |
| $p$ | frame-to-CM vector | 3-vec | in eq:A.4 |
| $\phi(\mathcal{O}_x,\mathcal{O}_y)$ | spatial transformation matrix | $6\times 6$ | eq:A.2, $[[I,\tilde l],[0,I]]$ |
| $l(\mathcal{O}_x,\mathcal{O}_y)$ | vector between frame origins | 3-vec | |
| $\phi(i,j)$ | composite transform, $i>j$ | $6\times 6$ | semigroup eq:4.11 |
| $\mathcal{E}_\phi$ | shift operator | $6n\times 6n$ | strictly lower block-bidiag eq:4.9 |
| $\phi$ (operator) | $[I-\mathcal{E}_\phi]^{-1}$ | $6n\times 6n$ | lower block-triangular |
| $M$ (operator) | $\mathrm{diag}\{M(k)\}$ | $6n\times 6n$ | |
| $H$ (operator) | $\mathrm{diag}\{H(k)\}$ | $\mathcal{N}\times 6n$ | |
| $P(k)$ | articulated-body spatial inertia (Riccati) | $6\times 6$ | eq:5.1 |
| $P^+(k)$ | propagated inertia $\bar\tau P$ | $6\times 6$ | |
| $D(k)$ | hinge inertia $HPH^*$ | $m\times m$ | invertible |
| $G(k)$ | Kalman gain $PH^*D^{-1}$ | $6\times m$ | |
| $K(k{+}1,k)$ | shifted gain $\phi G$ | $6\times m$ | first-subdiagonal of $K$ |
| $\bar\tau(k)$ | projector $I-GH$ | $6\times 6$ | |
| $\psi(k{+}1,k)$ | articulated-body transform $\phi\bar\tau$ | $6\times 6$ | |
| $\psi$ (operator) | $[I-\mathcal{E}_\psi]^{-1}$ | $6n\times 6n$ | analog of $\phi$ |
| $z(k),z^+(k)$ | residual spatial force | 6-vec | eq:5.11a |
| $\varepsilon(k)$ | innovation/residual | $m$-vec | eq:5.11a |
| $v(k)$ | scaled residual $D^{-1}\varepsilon$ | $m$-vec | eq:5.10 |
| $\mathcal{O}_k^+,\mathcal{O}_k^-$ | hinge frames on the two coupled clusters | frame | |
| $\Pi(\theta)$ | kinematic map $\dot\theta=\Pi\beta$ | $m\times m$ | assumed $=I$ |
