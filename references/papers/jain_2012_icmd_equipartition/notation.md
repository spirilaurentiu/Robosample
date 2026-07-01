# Notation — Jain 2012, Equipartition Principle for ICMD

Convention: `x*` denotes the **transpose** of vector/matrix x (not conjugate).
`≜` means "is defined as". `δ_{i=j}` in the paper is the **Kronecker delta**
(1 if i=j, else 0) despite being verbally called a "Dirac delta". Statistical
mechanics uses `k` = Boltzmann constant, `T` = temperature, so `kT` is an energy;
each independent quadratic momentum/modal DOF carries mean `kT/2` of kinetic
energy.

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| n | number of atoms | scalar int | Cartesian model has 3n DOF |
| N | number of ICMD configuration DOF | scalar int | BAT: N=3n; TAMD/constrained: N<3n |
| k | Boltzmann constant | energy/temperature | statistical mechanics |
| T | thermodynamic temperature | temperature | simulation temperature |
| kT | thermal energy scale | energy | mean kinetic energy per DOF = kT/2 |
| h | Planck's constant | action | in partition-function prefactor 1/h^N |
| q, p | generalized coords / conjugate momenta | R^n each | CANONICAL pair (general theory) |
| H(q,p) | Hamiltonian | energy | eq:2, eq:10 |
| L(q,q̇) | Lagrangian | energy | eq:1 |
| Z(T) | canonical partition function | dimensionless | eq:2, eq:11, eq:20 |
| ⟨f⟩ | canonical ensemble average of f | same as f | eq:3, eq:21, eq:23 |
| y_i, y_j | any phase-space coordinate pair | scalar | eq:4 (canonical only) |
| δ_{ij} | Kronecker delta | 0/1 | =1 iff i=j |
| I_m | m×m identity matrix | R^{m×m} | I_{2n}, I_N, I_{3n} |
| x | Cartesian atom positions | R^{3n} | absolute coordinates |
| ẋ | Cartesian atom velocities | R^{3n} | eq:28 |
| m_i | mass of atom i | mass | Cartesian |
| M | Cartesian mass matrix | R^{3n×3n}, diagonal, constant | atom masses on diagonal |
| M^{1/2} | sqrt of Cartesian mass matrix | R^{3n×3n}, diagonal | sqrt of atom masses |
| U(x), U(θ) | potential energy | energy | eq:6, eq:10 |
| θ | ICMD generalized coordinates | R^N | BAT/TAMD internal coords |
| θ̇ | ICMD generalized velocities | R^N | physical velocities |
| p (ICMD) | ICMD conjugate momenta | R^N | p = M(θ)θ̇, eq:9 |
| M(θ) (script M) | ICMD mass matrix | R^{N×N}, SPD, dense, config-dependent | M = J* M J, eq:29 |
| M^{-1}(θ) | inverse ICMD mass matrix | R^{N×N}, SPD | eq:10, eq:16 |
| K_e | ICMD kinetic energy | energy | eq:8 |
| m(θ) | mass-matrix factor, M = m m* | R^{N×N}, invertible | NOT necessarily lower-triangular/Cholesky |
| l(θ) | inverse factor, l = m^{-1} | R^{N×N} | M^{-1} = l* l, eq:16 |
| ν (nu) | modal velocity coordinates | R^N | ν = m*θ̇ = l p, eq:17; NONCANONICAL |
| U'(θ) | corrected/effective potential | energy | U + U_c, eq:22 |
| U_c(θ) | Fixman-type correction | energy | ½ ln det{M(θ)}, eq:22 |
| J(θ) | Cartesian↔ICMD Jacobian | R^{3n×N} | ẋ = J θ̇, eq:28; square/invertible only for BAT |
| P | constrained velocity mapping | R^{N×3n} | θ̇ = P ẋ, eq:37; valid choice eq:39 |
| H (script H) | hinge-articulation spatial operator | block operator | eq:31; per-body hinge map H(k) |
| φ (phi) | rigid-body propagation operator | block lower-triangular | φ*(k+1,k) transposed transform, Table 1 |
| M (script, in eq:31) | link spatial-inertia operator | block-diagonal | Newton–Euler factorization |
| K (script K) | spatial gain operator | block operator | in [I+HφK] |
| ψ (psi) | articulated-body propagation operator | block operator | in [I−HψK], eq:31 |
| D (script D) | articulated-body inertia | block-DIAGONAL, R^{N-ish blocks} | D^{1/2}, D^{-1/2} block-wise |
| G(k) | articulated-body gain (per body) | R^6-ish | Table 1 recursion; from K/ψ (ref 19/25) |
| V(k) | combined angular+linear velocity of frame k | R^6 | intermediate in Table 1 recursion |
| V^+(k) | propagated velocity from child | R^6 | Table 1 |
| α_i, γ_i | integration limits for coordinate i | coordinate units | geometry-determined |
