# Simbody - notation

Units are SI (m, kg, s, N, Pa, rad). Generalized coordinates need not equal generalized
speeds in number (quaternion mobilizers have nq > nu).

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| t | time (independent variable) | s | scalar |
| y | continuous state vector | R^dim | y = {q, u, z} |
| d | discrete state variables | - | fixed during a continuous interval |
| q | full generalized coordinates | R^nq | pose of all bodies vs parents |
| u | full generalized speeds | R^n | basis for equations of motion |
| z | auxiliary continuous variables | R^nz | muscle/controller state, first-order ODE |
| q_i, u_i | coordinates/speeds of mobilizer i | R^nq_i, R^nu_i | 0 <= nu_i <= 6 dofs per mobilizer |
| N, N_i | kinematic map qdot = N u | block-diagonal, N_i is nq_i x nu_i | invertible when nq_i = nu_i |
| M(q) | internal-coordinate mass matrix | R^{n x n}, symmetric PD | never formed explicitly (O(n) recursive M^{-1} v) |
| f_applied | applied generalized forces | N (generalized) | gravity + body forces/torques + direct generalized forces |
| f_inertial(q,u) | Coriolis/gyroscopic generalized forces | N (generalized) | velocity-dependent |
| f_constraint | generalized constraint force = G^T lambda | N (generalized) | |
| lambda | Lagrange multipliers (constraint forces) | R^m | least-squares solution if G rank-deficient |
| p | holonomic (position) constraints | R^{mp} | p(t;q) = 0 |
| v | nonholonomic (velocity) constraints | R^{mv} | v(t,q;u) = 0 |
| a | acceleration-only constraints | R^{ma} | a = A udot - b_a = 0 |
| P | holonomic velocity Jacobian | (∂p/∂q) N | |
| V | nonholonomic Jacobian | ∂v/∂u | |
| A | acceleration constraint Jacobian block | | |
| G | stacked acceleration constraint Jacobian | R^{m x n} = [P V A]^T | generally poorly conditioned/singular (redundant constraints) |
| b | acceleration constraint residual | R^m = [b_p b_v b_a]^T | terms independent of udot |
| Y | condensed constraint matrix G M^{-1} G^T | R^{m x m} | cost O(mn+m^2) to form |
| Y^+ | pseudoinverse of Y | | via complete orthogonal factorization (QTZ), ~5x faster than SVD |
| n | number of generalized speeds/coordinates | int | O(n) recursive solver |
| m | number of constraint equations | int | typically n >> m (internal-coord joints eliminate constraints); O(m^3) term |
| alpha | accuracy request | dimensionless | ~ relative tolerance; alpha = 10^{-n_digits}; real-time uses 1-10% |
| W, T | diagonal weighting matrices | | map heterogeneous errors (length/angle/velocity) to unit error |
| epsilon_y | integrator absolute error estimate for y | | RMS norm used in step acceptance |
| x | normal contact deformation | m | x > 0 when surfaces contacting |
| xdot | normal approach velocity | m/s | |
| R, R_i | composite / per-body relative radius of curvature | m | R = R1 R2 / (R1+R2) |
| E_i | Young's modulus of material i | Pa | |
| nu_i | Poisson ratio of material i | dimensionless | |
| E_i* | plane-strain modulus = E_i/(1-nu_i^2) | Pa | |
| E* | composite elastic modulus | Pa | Hertz: 2/3-power rule; EFM: E1*E2*/(E1*+E2*) |
| sigma | eccentricity correction factor | dimensionless | =1 for circular contact, else Eq. 28 |
| s_1, s_2 | deformation split fractions (s1+s2=1) | dimensionless | locate contact point P |
| c_i, c* | per-material / effective dissipation coefficient | s/m | c* = c1 s1 + c2 s2 (Hertz weighting) |
| f_Hz | Hertz stiffness force | N | (4/3 sigma R^{1/2} E*) x^{3/2} |
| f_HC | Hunt-Crossley dissipation force | N | signed; = (3/2) f_Hz c* xdot |
| mu(v) | effective friction coefficient | dimensionless | Stribeck curve, C2 spline (2 quintics + slide) |
| v (friction) | slip speed in contact plane | m/s | v = |v_rel| |
| k, h (EFM) | per-triangle spring stiffness, elastic layer thickness | N/m, m | k from triangle area, E*, h |
| a, b (ellipse) | semi-major, semi-minor axes of contact ellipse | m | k = a/b >= 1 |
| A, B (paraboloid) | principal semi-curvatures of separation paraboloid | 1/m | B >= A |
| K(m), E(m) | complete elliptic integrals of 1st, 2nd kind | dimensionless | m = 1 - (1/k)^2 |
