# Notation - Jain et al. 2013, Fixman potential for branched molecules

Conventions: reduced/consistent MD units. Temperature enters as kT (k = Boltzmann constant). `*` superscript = matrix transpose throughout (NOT complex conjugate). "triangle-equals" = "defined as". Torsion angles measured in radians for pdf normalization (uniform pdf = 1/2pi); histograms reported in degrees (bin dalpha = 7.2 deg). Spatial operators act on 6-dim spatial (rotation+translation) vectors per rigid cluster; base cluster carries 6 rigid-body DOF.

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| n | number of atoms | integer | 3n total BAT coordinates |
| N (script N) | number of unconstrained coordinates (torsions) | integer | alpha has N entries; (3n-N) constrained |
| alpha | unconstrained BAT coordinates (assumed torsions) | rad | the free DOF in constrained model |
| beta_i | ith torsion angle | rad | flexible-model uniform pdf 1/2pi |
| q | coordinates to be constrained (bond lengths, bond angles) | Angstrom / rad | frozen at q_0 in constrained model |
| q_0 | fixed constrained-coordinate values | - | hard holonomic constraint values |
| p | canonical momenta (full BAT set) | - | conjugate to (alpha,q) |
| H (Hamiltonian) | Hamiltonian | energy | eq:1 |
| U | force-field potential energy | kcal/mol | bond+angle(+nonbonded) |
| U_f | Fixman compensating potential | kcal/mol (energy) | eq:11/13/17 |
| U' | Fixman-corrected potential | energy | U(alpha,q0)+U_f |
| T(i) (script T) | Fixman torque on coordinate i | energy/rad | = -dU_f/dalpha_i |
| M_B | BAT-coordinate mass matrix (unconstrained) | mass units | R^{3n x 3n}, dense, config-dependent |
| M (script M) | constrained-model mass matrix | mass units | R^{N x N}, alpha-alpha sub-block of M_B |
| M_c | Cartesian mass matrix | mass units | R^{3n x 3n} diagonal, atom masses on diagonal |
| J_B | Jacobian BAT->Cartesian | - | R^{3n x 3n} |
| m_i (fraktur m) | mass of ith atom | amu | beads: 14 amu in serial-chain tests |
| d_i | bond lengths | Angstrom | (n-1) of them; 1.54 A in serial tests |
| theta_i | bond angles | rad | (n-2) of them; 90 deg serial, 109 deg C5 |
| theta_ex | overall molecule orientation coordinate | rad | in Go-Scheraga det J_B |
| Z(T) | partition function | - | eqs:2,3,9 |
| rho | probability density | 1/rad | torsion pdf, target uniform 1/2pi |
| c_1,c_2,c_3,c_5,c_f | constants (masses/bond-length/normalization dependent) | - | absorb config-independent factors |
| k | Boltzmann constant | energy/temperature | kT is thermal energy |
| T | temperature | K | simulations at 300 K |
| H (operator) | hinge articulation operator | block-diagonal | torsional axes per DOF |
| phi | rigid-body force-propagation operator | lower-triangular | phi(i,j): force cluster j -> cluster i |
| M (operator) | link spatial inertia operator | block-diagonal | 6x6 spatial inertia per cluster |
| D, D(i) | articulated-body (AB) inertia | block-diagonal | D(i) scalar for i>=1, D(0) is 6x6 |
| K | AB gain operator | - | in [I+HphiK] factorization |
| psi | AB force propagation (phi with feedback) | - | in M^{-1} factorization |
| P(i) | AB matrix per cluster (GNEIMO by-product) | 6x6 | in torque eq:22 |
| Upsilon(i) | matrix from recursive scatter algorithm | 6x6 | extra recursion from base cluster |
| H-tilde_omega(i) | rotational part of hinge map | 6x6 | in torque eq:22 |
| h(i) | torsion hinge axis | 3-vector | unit rotation axis |
| F(A) | "uncross" operator | 3x3 -> R^3 | returns v s.t. (A-A^*) is skew(v) |
| Q_{jk} | 3x3 sub-blocks of P(i)Upsilon(i) | R^{3x3} | eq:23 |
| gamma | Langevin damping coefficient | 1/fs | 0.01/fs in simulations |
| dW | Wiener process increments | - | independent per DOF |
| R_flex, R_tor, R_Fix | RMS deviation of torsion pdf from uniform | dimensionless | FLEXIBLE / TORSIONAL / FIXMAN metrics |
