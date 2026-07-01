# Checks / fixtures - Duane et al. 1987, Hybrid Monte Carlo

The paper is primarily algorithmic; concrete numbers come from the compact-QED test runs. These are qualitative/coarse benchmarks (read from figures), not high-precision tables, so treat tolerances loosely. The most valuable "checks" for a molecular HMC port are the algorithmic invariants, listed first.

## Algorithmic invariants (implementation regression tests)

- **Detailed balance / correct stationary distribution:** for ANY guidance Hamiltonian $H'$ (even $H' \neq H$), the sampler converges to $P_S(\phi) \propto \exp[-S(\phi)]$. Given a leapfrog integrator that is reversible (eq:14) and area-preserving, the accept rule eq:12 makes the chain exact regardless of step size $\delta\tau$. Test: sample a 1-D/2-D Gaussian or double well; histogram must match $\exp[-S]$ within MC error independent of $\delta\tau$.
- **Energy conservation limit:** when $H' = H$ and integration is exact, $\delta H = 0$ and acceptance $P_A = 1$ (reduces to the standard hybrid/MD algorithm).
- **Leapfrog reversibility:** integrating $n$ steps forward then flipping momentum and integrating $n$ steps must return to the start (to round-off). Map is area-preserving (Jacobian = 1) for any $\delta\tau$.
- **Leapfrog error orders:** half-steps err at $\mathcal{O}(\delta\tau^2)$, intermediate steps at $\mathcal{O}(\delta\tau^3)$ per step; global energy error $\langle\delta H\rangle = \mathcal{O}(\delta\tau^2)$.
- **Acceptance vs energy:** $\langle\exp(-\delta H)\rangle = 1$ (fluctuation identity) is a strong regression check for a correct HMC implementation.

## Compact-QED numerical fixtures

Setup: quenched (pure gauge) and dynamical-electron compact QED, Wilson gauge action + staggered fermions, 4-D lattices.

| quantity | setup | value |
|---|---|---|
| optimal leapfrog step $\delta\tau$ | $8^4$ lattice, quenched, coupling $\beta = 0.97$ | $\approx 0.1$ |
| effective step size $a\,\delta\tau$ at optimum | same | $\approx 0.06$ |
| max "safe" standard-hybrid step size | for comparison | $0.01$-$0.02$ |
| step-size scaling for constant acceptance | lattices $4^4$, $8^4$, $12^4$ | $\delta\tau \propto 1/L$ |
| acceptance-rate model | vs lattice size and step | $a \propto \exp(-L^2\delta\tau^2)$ |
| dynamical-fermion runs (Fig 3) | $8^4$ lattice, electron mass $= 0.25$ (lattice units), $\beta = 0.8$ | plaquette & chiral condensate $\langle\bar\psi\psi\rangle$ vs $\delta\tau$; step-size dependence similar to pure gauge |
| $\beta - \beta'$ acceptance peak (Fig 4) | $8^4$ quenched, $\beta = 0.97$, $\delta\tau = 0.15$ | acceptance peaks at $\beta \neq \beta'$ (guidance $\neq$ acceptance improves acceptance) |

Notes:
- $\beta = 0.97$ is the QED coupling constant in the test (not inverse temperature).
- Pure-gauge $8^4$, $\beta = 0.97$ results agree with prior high-statistics data (ref [10], Lang 1986).
- The new HMC needs no extrapolation to $\delta\tau \to 0$ (bias-free at finite step), unlike standard Langevin/MD hybrids.
