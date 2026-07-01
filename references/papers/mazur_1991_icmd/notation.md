# Notation — Mazur 1991, explicit internal-coordinate MD

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| θ_k, θ_i, θ_m | generalized internal coordinate (bond length, valence/planar angle, torsion/phase angle) | Å (bond) or rad (angle) | free (non-frozen) coordinates only |
| θ̇_k | generalized velocity | Å/fs or rad/fs | dot = d/dt |
| θ̈_i | generalized acceleration | Å/fs² or rad/fs² | solved for from eq:18 |
| n, n_var | total number of free generalized coordinates | int | global numbering 1..n_var |
| n_α | number of variables in chain V_α that determine atom α's position | int | local per-atom count |
| V_α | ordered chain of generalized coordinates determining position of atom α | index list | node order per Fig. 2 |
| d_k | set of atoms whose position depends on θ_k | atom-index set | used as summation range |
| α | atom index | int | |
| m_α | mass of atom α | amu | |
| **r**_α | Cartesian position of atom α | Å, R^3 | global frame |
| ṙ_α | Cartesian velocity of atom α | Å/fs, R^3 | |
| **r**_i^θ | radius vector of the node of variable i | Å, R^3 | |
| **r**_{α/k} | relative position **r**_α − **r**_k^θ (atom α relative to node k) | Å, R^3 | |
| **e**_i, **e**_i^α | unit vector of variable i | dimensionless, R^3, ‖e‖=1 | rotation axis (angle) or translation direction (bond) |
| S_i | variable-type indicator | {0,1} | 1 = angle variable (rotation), 0 = bond length (translation) |
| φ | torsion angle | rad | first in node order |
| Φ | phase angle (dihedral added to torsion to branch from a node) | rad | shares node with φ |
| ω | planar / valence angle | rad | after torsion/phase in node order |
| b | bond length | Å | last in node order |
| L | Lagrangian, L = T − U | kcal/mol | |
| T (energy context) | kinetic energy T(θ, θ̇) | kcal/mol | |
| U | potential energy U(θ) | kcal/mol | ECEPP + CHARMM; intra-rigid-body terms discarded |
| ∂U/∂θ_k | generalized force on coordinate k | kcal/mol/(Å or rad) | RHS of eqs. 2,16,18 |
| a_ki | mass-matrix element (coeff. of θ̈_i in eq. k) | mass-weighted | symmetric, positive-definite |
| b_ki | coefficient of θ̇_i² (centrifugal term) | | |
| c_kim | coefficient of θ̇_m θ̇_i (Coriolis cross term) | | |
| T (thermo) | temperature | K | from eq:temperature |
| N | number of internal degrees of freedom | int | used in T = 2⟨K⟩/(N k_B), NOT 3×atoms |
| k_B | Boltzmann constant | kcal/mol/K | |
| K | instantaneous kinetic energy | kcal/mol | |
| δ_E | relative RMS energy-conservation error | dimensionless | eq:deltaE |
| h | integration time step | fs | |

## Conventions

- **BKS-tree**: molecular system = single unified tree of rigid bodies rooted at the
  global frame origin; cycles disconnected; virtual atoms/bonds connect molecules
  and enforce tree topology. Robosample's "robot" (kinematic tree of rigid bodies)
  is the same object.
- **Rigid body** = group of atoms depending on the same set of internal variables;
  interactions internal to a rigid body are always discarded from U.
- **Variable ordering at a node (Fig. 2), small index first:** torsion φ and phase
  Φ, then planar angle ω, then bond length b. Ordering is load-bearing: eqs. 6-16
  assume that if variable i influences the unit vector of variable j then i < j.
- **Indicator convention:** S = 1 ⇒ rotational DOF (angle), contributes **e** × **r**
  terms; S = 0 ⇒ translational DOF (bond length), contributes plain **e** terms.
- Vector notation in eq:16: adjacent parenthesized vectors are dot products;
  **a** × **b** × **c** ≡ **a** × (**b** × **c**).
- Reduced/units: energies in kcal/mol, lengths in Å, time in fs (consistent with the
  0.5 fs step and kcal/mol energies reported).
