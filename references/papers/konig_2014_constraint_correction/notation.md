# Notation - Constraint free-energy correction (Konig & Brooks 2014)

Units: energies kcal/mol; lengths Angstrom (A); angles radians (unless stated
degrees in a table); force constants kcal/mol/A^2 (bonds) or kcal/mol/rad^2 (angles).
Thermodynamic convention: `β = 1/(k_B T)`, T = 300 K in all simulations.

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| q, Δq | internal coordinate; deviation from its energy minimum | A (bond) or rad (angle), scalar | Δq = 0 at minimum |
| K | harmonic force constant of a bond or angle | kcal/mol/A^2 or kcal/mol/rad^2 | force field parameter |
| U, U_0 | potential energy; zero-point/reference energy | kcal/mol | U_0 = value at minimum |
| U'(Δq_cons), U''(Δq_cons) | 1st/2nd derivative (gradient, Hessian) of U at constrained conf | kcal/mol/A, kcal/mol/A^2 | from numerical differentiation |
| β | inverse temperature 1/(k_B T) | mol/kcal | T = 300 K |
| Z, Z^{h.o.}, Z_cons, Z_add | partition functions: total, harmonic-osc, constrained, released-DOF | dimensionless | |
| G, G^{h.o.}, G_cons | absolute free energies | kcal/mol | |
| ΔG_cons | free energy of IMPOSING a constraint | kcal/mol | = ΔH + ΔG_harm (+ΔG_Jacobian in Cartesian) |
| ΔG_rls | free energy of RELEASING a constraint = -ΔG_cons | kcal/mol | |
| ΔH | enthalpic contribution (potential-energy rise to minimum) | kcal/mol | >= 0 always |
| ΔG_harm | vibrational-entropy contribution | kcal/mol | <= 0 always (no entropy gained adding a bond) |
| ΔG_Jacobian | internal->Cartesian Jacobian contribution | kcal/mol | can be + or - |
| Δq_cons | constrained value of the internal coordinate (deviation) | A or rad | where the constraint is set |
| n | number of atoms | int | |
| m | number of constrained DOFs to correct | int | |
| C | reduced orthonormal basis set | 3n x m matrix | mass-weighted, normalized; columns c_i |
| c_i | basis vector for the i-th constrained DOF | length-3n vector | ∂x/∂r_i (bond) or ∂x/∂θ_i (angle) |
| F | full mass-weighted energy 2nd-derivative (Hessian) | 3n x 3n | |
| H | reduced Hessian = C† F C | m x m | off-diagonals = coupling between constraints |
| g | reduced-basis force/gradient = C† f | length-m vector | |
| f | full-system force vector | length-3n vector | |
| Λ | diagonal eigenvalue matrix of H (force constants of modes) | m x m diagonal | from Λ = U† H U |
| U (eq:18) | unitary eigenvector matrix of H | m x m | NOTE: overloaded symbol - NOT potential energy |
| ‖...‖_1 | sum over all elements (L1-style) | scalar | |
| J_r | bond Jacobian factor | A^2 | = r_ik^2 |
| J_θ | angle Jacobian factor (linear chain) | dimensionless | = sin θ |
| J_θ' | branching-angle Jacobian factor | dimensionless | = 1/sin θ' (reciprocal sine); NOT implemented in code |
| r_ik | bond length between atoms i and k | A | |
| θ_jkl | bond angle at vertex k (atoms j,k,l) | rad | |
| <...>_cons | ensemble average over the constrained trajectory | - | Zwanzig exponential average |
| λ | alchemical coupling parameter (ethane<->methanol) | dimensionless [0,1] | dual-topology mixing |

## Sign / regime conventions (implementation-critical)

- ΔH >= 0: constrained conformation potential energy cannot be below the minimum.
- ΔG_harm <= 0 always.
- ΔG_Jacobian may be positive (bond stretched, e.g. polar H attracted to Cl-) or negative.
- Method is valid ONLY for **hard** DOFs (bonds, angles) near equilibrium; NOT for soft DOFs (dihedrals, translations, rotations).
- Correction ordering: correct **bonds (CBND) before angles (CANG)**, else the angle basis vectors use stale bond lengths.
- Hessian may be approximated by the force-field force constant (fast path), which neglects vibrational-entropy change (~0.1% error reduction loss).
- CHARMM implementations: VIBRAN REDUce = full C (keeps coupling); RAYLeigh = per-vector (neglects inter-constraint coupling).
