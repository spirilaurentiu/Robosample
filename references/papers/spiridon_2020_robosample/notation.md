# Notation - Robosample (Spiridon et al. 2020)

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| N | number of atoms | scalar int | - |
| m | number of constraints | scalar int | - |
| N_f | number of generalized (free) coordinates | scalar int | N_f = 3N - m |
| φ_f | free / flexible generalized coordinates | vector R^{N_f} | evolve in reduced dynamics |
| φ_v | constrained ("velocity") generalized coordinates | vector | frozen in reduced dynamics |
| {φ} | full generalized coordinate set | R^{3N×3N} transform block | split into φ_f (free) and φ_v (constrained) |
| p_f, p_v | conjugate momenta to φ_f, φ_v | vectors | - |
| H | Hamiltonian (kinetic + potential [+ Fixman]) | energy | - |
| U | potential energy | energy (kJ/mol in sims) | force field (ff14SB, GLYCAM06) |
| U' | Fixman correcting potential | energy | eq:7 |
| M | diagonal Cartesian mass matrix | 3N×3N diagonal | atom masses |
| M_tot | full (Cartesian-space) mass matrix tensor | N_f×N_f (or 3N×3N) | M_tot = J^T M J |
| M(φ_f) | reduced-coordinate mass matrix tensor | N_f×N_f | body-frame articulated inertia |
| \|M\| | determinant of a mass matrix tensor | scalar | appears as \|M\|^{-1/2} in marginals |
| J | Jacobian of Cartesian↔generalized transform | 3N×N_f | J = {dx_i/dφ_f} |
| x_i | Cartesian position of atom i | R^3 | - |
| ρ | probability density (micro/config state) | - | canonical ensemble |
| β | inverse temperature | 1/energy | β = 1/(kT) |
| kT, k | thermal energy, Boltzmann constant | energy | - |
| T | temperature (also used for HMC integration time) | K (or fs) | context-dependent |
| ε | MD / HMC integration time step | fs | - |
| L | HMC trajectory length (number of integration steps) | int | - |
| T = L·ε | HMC integration time | fs | - |
| φ, ψ | alanine dipeptide backbone torsions | degrees | Ramachandran angles |

## Conventions / gotchas
- **Reduced units none stated for energy**: simulation energies reported in kJ/mol; free-energy surfaces in kJ/mol.
- **Kinetic-energy form**: eqs 2/5/8 print kinetic energy as `p^T M p` (no 1/2, no explicit inverse). Physically KE = (1/2) p^T M^{-1} p; when implementing, treat the `M` symbol in the momentum-quadratic term as the inverse mass-matrix operator. See CHECK notes in equations.md.
- **Ramachandran dynamics** = only φ and ψ torsions mobile; all other DOFs constrained.
- **Gibbs "world"** = a defined set of rigid bodies + joint types mapped onto the same molecular graph; simulation alternates worlds so every DOF is sampled by at least one move (ergodicity). At least one fully-flexible world per cycle ensures ergodicity.
- **Joint types**: Cartesian (fully flexible), Torsional (TD), Cylindrical (Cyl, 2 DOF: translation + rotation about that axis), Spherical/Ball (3 rotational DOF, Euler angles or quaternions), Weld, Slider, Universal, Free.
- **Mobilizer** = joint connecting a child rigid body to a parent; described by 4 reference frames (one per body + two joint frames, parent-fixed and child-mobile). Robosample default: joint frames placed along bonds at atom centers (random placement degrades transition rates).
- **Fixman potential/torque** requires |M_tot| and |M| determinants; computed via O(N) Spatial Operator Algebra rather than O(N^3) direct inversion.
- **CDHMC** = Constrained Dynamics Hamiltonian Monte Carlo used as a Gibbs sampling move.
- Optimal HMC acceptance rate target: 0.651.
