# Notation - Sohl-Dickstein 2014, LAHMC

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| x | position / sample | real vector, R^N | target variable |
| v | momentum (auxiliary) | real vector, R^N | mass set to 1 (Eq 6) |
| N | dimension of state space | positive integer | |
| ζ (zeta) | extended state {x, v} | pair of R^N vectors | phase-space point |
| E(x) | potential energy function | scalar | target = exp(-E)/Z (unit temperature, β_thermo=1) |
| Z | partition function / normalizer | scalar | |
| p(x), p(v), p(ζ) | densities | scalar | Boltzmann form |
| H(ζ) | Hamiltonian = E(x) + (1/2)v^T v | scalar | total energy |
| F | momentum flip operator | involution, v -> -v | own inverse; unit Jacobian; preserves H |
| L, L(ε, M) | leapfrog integrator operator | deterministic map | M steps of size ε; unit Jacobian; L^-1 = FLF |
| L^a | a-fold application of L | deterministic map | a in 1..K |
| ε (epsilon) | leapfrog step length | scalar | experiments: ε = 1 |
| M | leapfrog steps per L application | positive integer | experiments: M = 10 |
| K | max number of L applications in LAHMC | positive integer | experiments: K = 4 |
| R(β) | momentum randomization operator | stochastic map | partial refresh (Eq 15) |
| β (beta) | momentum noise fraction | scalar in [0, 1] | β=1 full refresh, β=0 none; not thermodynamic β |
| n | Gaussian noise vector | R^N | n ~ N(0, I) |
| α (alpha) | fraction of momentum randomized per unit sim time | scalar | β = α^(1/(εM)) |
| π_accept | MH acceptance probability (standard HMC) | scalar in [0,1] | Eq 18 |
| π_{L^a}(ζ) | prob of transition ζ -> L^a ζ (LAHMC) | scalar in [0,1] | greedy, Eq 25 |
| π_F(ζ) | prob of momentum-flip transition (LAHMC) | scalar in [0,1] | residual, Eq 27 |
| ζ^(t,s) | state at sampling step t, substep s | phase-space point | |
| d | (future-work) direction indicator | {-1, +1} | p(d=1) = 1/2 |

Conventions:
- Reduced/natural units: mass = 1, temperature implicit (energy in exp(-E)).
- p(ζ) is invariant under F because H(Fζ) = H(ζ), so p(Fζ) = p(ζ).
- Density ratios reduce to energy differences: p(ζ')/p(ζ) = exp(H(ζ) - H(ζ')).
- LAHMC has NO Metropolis accept/reject step; the greedy π's are the sampler.
