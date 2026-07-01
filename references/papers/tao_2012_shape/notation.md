# Notation — Tao et al. 2012, SHAPE

Units follow CHARMM/leapfrog MD conventions (length Å, mass amu, energy kcal/mol,
time fs). The algorithm assumes a leapfrog (Verlet) integrator producing
half-step velocities.

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| N | total atoms in system | int | — |
| M | atoms in a single rigid structure | int | SHAPE cost independent of M |
| t | current time | fs | leapfrog step at t |
| Δt | integration time step | fs | 1.0 fs default; tested 1.0/1.5/2.0 |
| m_i | mass of atom i | amu | — |
| r_i(t) | global Cartesian position of atom i at t | R^3 | pre-step positions |
| r_i^non(t+Δt) | free-MD (unconstrained) global position at t+Δt | R^3 | output of the unconstrained integrator |
| r_i^rig,(n)(t+Δt) | rigid-constrained global position, iteration n | R^3 | final desired trajectory (eq. 27) |
| r_COM(t) | center of mass at t | R^3 | mass-weighted mean (eq. 1) |
| r_COM^non(t+Δt) | COM of unconstrained positions at t+Δt | R^3 | eq. 2 |
| superscript b | body-fixed local frame (origin at COM, no rotation) | — | r^b = r − r_COM (eq. 3,4) |
| r_i^b(t) | body-fixed position at t | R^3 | reference for rotation |
| r_i^non,b(t+Δt) | body-fixed unconstrained position at t+Δt | R^3 | eq. 4 |
| v_i^b, p_i^b | body-fixed velocity / momentum at t+Δt/2 | R^3 | finite difference (r^b(t+Δt)−r^b(t))/Δt |
| L^non(t+Δt/2) | unconstrained angular momentum at half step | R^3 | eq. 5; conserved target |
| L^rig(t+Δt/2) | rigid-structure angular momentum | R^3 | I·ω (eq. 7) |
| L'^non | reduced (unrotated) unconstrained ang. mom. | R^3 | uses time-t coords (eq. 18) |
| L'^rig,(n) | reduced rigid ang. mom., iteration n | R^3 | eq. 19 |
| I | moment of inertia tensor | R^{3×3}, symmetric | body-fixed, computed from r_i^b(t) (eq. 6) |
| I^-1 | inverse inertia tensor | R^{3×3} | eq. 10, 22, 26 |
| ω, ω^rig | angular velocity vector | R^3, rad/fs | ω = (ω_x,ω_y,ω_z) |
| ω^rig,(n) | angular velocity at iteration n | R^3 | fixed-point iterate |
| ω̂ | skew-symmetric (hat) matrix of ω | R^{3×3} | ω̂ v = ω × v (eq. 11) |
| θ | rotation-angle vector = ω Δt | R^3, rad | ‖θ‖ = rotation angle |
| θ̂ | skew matrix of θ = ω̂ Δt | R^{3×3} | eq. 12 |
| R, R^(n) | rotation matrix = exp(θ̂) | SO(3), R^{3×3} | eq. 12–14 |
| (R^(n))^{1/2} | matrix square root of R^(n) | R^{3×3} | half rotation; approx ½(1+R^(n)) (eq. 25) |
| δR | small deviation R^(n) − 1 | R^{3×3} | small-Δt regime (eq. 23) |
| n | outer/inner iteration index | int | ~3 iterations to double precision |
| τ | torque on the rigid structure | R^3 | derivation of accuracy (eq. 30) |
| 1 | identity matrix | R^{3×3} | — |
| tolerance | convergence tolerance on ‖L^rig,(n) − L^non‖ | scalar | 1e-7 used in tests |

Conventions / notes:
- Linear momentum is conserved exactly by working in the body-fixed COM frame;
  only angular momentum requires iteration.
- `×` is the vector cross product; `·` is matrix–vector / matrix–matrix product.
- The half-step velocity is the leapfrog finite difference of body-fixed
  coordinates, `(r^b(t+Δt) − r^b(t))/Δt`.
- Convergence normally reached in 3 iterations (single rigid body). With atoms
  shared between rigid bodies an outer cyclic loop over constraints is required.
