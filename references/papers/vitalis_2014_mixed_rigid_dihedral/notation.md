# Notation - Vitalis & Pappu 2014

Reduced/convention notes:
- `beta = 1/(k_B T)`. Temperature in K; energies in kcal/mol; lengths in Angstrom; time in fs/ps; mass in Da (atomic mass units).
- Kinetic energy has TWO definitions here: the true Cartesian KE `(1/2) p^T M^{-1} p` and the scheme's **modified** KE `(1/2) omega^T I_D omega` (eq:5). They differ because I_D is the diagonal of the MMT, not the full MMT. Temperatures follow the same split (`<T_c>` from Cartesian KE, `<T>` from modified KE).
- Generalized coordinate space is MIXED: 6 rigid-body DOF per molecule (3 translation + 3 rotation) plus one dihedral hinge per rotatable torsion. Total K <= N_at.
- Recursions are O(N_at): inward (tip->base) for F_phi and I_D (eq:B1,B2); outward (base->tip) for Cartesian velocities (eq:B3).
- Half-integer time labels: current conformation at t_{1.5}; velocities at t_1; forces evaluated once per step at t_{1.5}; positions advance to t_{2.5}. Force-explicit (one U/grad evaluation per step).

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| N_at | number of atoms | int | up to system size |
| N_mol | number of molecules | int | |
| M (molecule count) | Z-matrix atom count per molecule | int | M-1 bonds, M-2 angles, M-3 dihedrals |
| K | number of free (generalized) DOF | int | K <= N_at; 6 rigid-body + hinges per molecule |
| r, r_i | Cartesian positions | R^{3N_at}, R^3 | space-fixed global frame |
| p | Cartesian momenta | R^{3N_at} | conjugate to r |
| M | mass matrix (Cartesian) | diagonal 3N_at x 3N_at | m_i on diagonal |
| m_i | atomic mass | Da | |
| U(r) | potential energy | kcal/mol | separable Hamiltonian |
| F_{r,i} | Cartesian force on atom i | kcal/mol/A | = -grad_i U |
| phi, phi_k | generalized coordinates (dihedrals + rigid-body) | R^K | tree topology, backward dependency |
| omega, omega_k | generalized velocities | rad/time (angular) or A/time (translation) | conjugate to phi via I_D |
| p_phi | generalized momenta | = G omega | |
| J | Jacobian dr/dphi | 3N_at x K | J_{kl} = dr_k/dphi_l |
| Y | reduced Jacobian (free DOF only) | 3N_at x K | eq:A2; sign of a_k depends on base |
| G = J^T M J | mass-metric tensor (MMT) | K x K | full (non-diagonal) |
| G_S | MMT of the free subsystem | K x K | constrained subsystem |
| I_D | diagonal effective-mass matrix | diagonal K x K | KEY approximation = diag(G) |
| I_kk, I_ii | diagonal effective-mass element of DOF k | Da (transl.) or Da*A^2 (rotational inertia) | time-dependent for hinges |
| I_{x/y/z} | rigid-body rotational inertia about lab axes through COM | Da*A^2 | |
| I_{phi i} | effective inertia of dihedral i | Da*A^2 | |
| F_phi, F_{phi,k} | projected/generalized force | eq:B1 | axis-projected net torque |
| a_k | unit rotation axis of DOF k | R^3, unit | bond vector for a hinge |
| b_k | reference point on rotation axis of DOF k | R^3, A | for rigid rotation: molecular COM |
| c_1 | unit normal in Z-matrix placement | R^3, unit | eq:A1 |
| alpha | bond angle (Z-matrix) | rad | across r_i, r_j, r_l |
| Omega_k | mass-weighted generalized velocity | = omega_k I_kk^{1/2} | dynamical variable in eq:17 |
| delta_t, dt | integration time step | fs | |
| Lambda | number of velocity sub-steps in eq:11 | int (multiple of 2) | typical 4 |
| tau_Lambda | sub-step size | = delta_t / Lambda | |
| t_1, t_{1.5}, t_2, t_{2.5} | half-integer time labels | | forces at t_{1.5}, velocities at t_1 |
| q_rot | rigid-body rotation quaternion | R^4 unit | [w, x, y, z], eq:12 |
| beta | inverse temperature | 1/(k_B T) | |
| k_B | Boltzmann constant | | |
| T, <T> | temperature from modified KE omega^T I_D omega | K | <T> = [k_B(6N_mol-3)]^{-1} <omega^T I_D omega> for water |
| T_c, <T_c> | temperature from Cartesian KE p^T M^{-1} p | K | <T_c> = [k_B(6N_mol-3)]^{-1} <p^T M^{-1} p> |
| alpha_T | global thermostat rescaling factor (Bussi) | scalar | replaces omega_k(t1) -> alpha_T omega_k(t1) |
| tau_T | thermostat coupling time | ps | |
| C_v | heat capacity at constant volume | cal/mol/K | |
| D | translational diffusion coefficient | 10^{-5} cm^2/s | |
| tau_rot | rotational autocorrelation time | ps | |
| eps_r | relative dielectric constant | dimensionless | from dipole fluctuations |
| N_h, N_s, N_1 | # alpha-helical H-bonds / segments / isolated helical residues | counts | helix-coil order params |
| Delta T | rotation-translation temperature imbalance | K | 2[3 k_B N_mol]^{-1}(<K_rot> - <K_trans>) |
