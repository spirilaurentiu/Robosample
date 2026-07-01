# Notation - Vaidehi & Jain 2015, ICMD / GNEIMO

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $\theta$ | generalized internal coordinates (torsion / BAT angles) | rad, vector R^N | reduced (constrained) dof set for torsional MD |
| $\dot\theta$ | generalized velocities | rad/s, R^N | |
| $\ddot\theta$ | generalized accelerations | rad/s^2, R^N | solved from eq 1 / eq 2 |
| $\mathcal{M}(\theta)$ | mass matrix (moment-of-inertia tensor) in reduced coords | N x N, symmetric positive-definite | configuration dependent; dense (cubic naive inverse) |
| $\mathcal{M}_B(\theta,q_0)$ | mass matrix in full BAT coords with frozen dofs at $q_0$ | full-dim square | used in Fixman potential (eq 5) |
| $C(\theta,\dot\theta)$ | Coriolis / centrifugal (velocity-dependent) force | generalized force, R^N | |
| $\mathcal{T}(\theta)$ | generalized forces / torques | R^N | from force field, reduced to hinge torques |
| $\mathcal{F}(\theta,\dot\theta)$ | Nosé-Hoover thermostat frictional force | generalized force, R^N | canonical ensemble term (eq 3) |
| $\eta$ | Nosé-Hoover thermostat dynamic variable | dimensionless | not a canonical coordinate |
| $\tau$ | thermostat mass parameter | time | set $\tau = 10\,\Delta t$ |
| $T$ | instantaneous temperature | K | |
| $T_B$ | thermostat / bath (target) temperature | K | |
| $k$ | Boltzmann constant | energy/K | |
| $\mathcal{U}_f(\theta)$ | Fixman compensating potential | energy | gradient = Fixman torque |
| $q_0$ | frozen (constrained) degree-of-freedom coordinates | | held fixed |
| $\mathcal{H}$ | joint/hinge map operator | spatial operator | GNEIMO factorization (Jain 2010) |
| $\psi$ | rigid-body transformation propagation operator | spatial operator | |
| $\mathcal{K}$ | Kalman-gain-like operator | spatial operator | |
| $\mathcal{D}$ | articulated hinge inertia (block diagonal) | invertible per hinge | |
| $\mathcal{P}$ | articulated-body inertia operator | spatial operator | |
| $a$ | Coriolis/gyroscopic spatial acceleration bias | spatial vector | |
| $b$ | gyroscopic spatial-force bias | spatial vector | |

## Conventions

- **Model:** tree-topology macromolecule = rigid "clusters" (single atom up to a
  whole domain) connected by flexible **hinges**; hinge has 1-6 dof (1 dof = pure
  torsion; 6 dof = stretch + bend + torsion about the connecting bond).
- **Torsional MD** = special case of ICMD/GNEIMO with all bond lengths and bond
  angles rigid (dof = torsions only).
- **Cost:** naive dense-matrix inverse of $\mathcal{M}$ is O(N^3); GNEIMO
  recursive spatial-operator solution (eq 2) is O(N) in number of dof.
- **Fixman:** corrects bias only for stiff, uncoupled bond angles treated as
  rigid; NOT for soft bond angles coupled to torsions/nonbond. GNEIMO-Fixman adds
  ~24% compute over the base solver.
- Angles reported in degrees in figures (e.g. bin size $d\theta = 18^\circ$).
