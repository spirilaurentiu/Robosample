# Notation - Metropolis et al. 1953

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| N | number of particles | int (up to several hundred; 224 or 56 used) | 2D system in a unit square |
| E | total configurational potential energy | energy | Eq. (1); half sum over ordered pairs, i≠j |
| V | pair potential between molecules | energy | spherically symmetric, two-body only |
| $d_{ij}$ | minimum-image distance between i and j | length | periodic (nearest image over surrounding squares) |
| $d_{AB}$ | minimum distance between particles A and B | length | shortest of A to any periodic image of B |
| k | Boltzmann constant | energy/temperature | classical statistics assumed |
| T | absolute temperature | temperature | canonical ensemble |
| $kT$ | thermal energy | energy | Boltzmann factor is $\exp(-E/kT)$; no $\beta$ used |
| F | observable of interest | varies | e.g. pressure, density |
| $\bar{F}$ | ensemble / Monte Carlo average of F | varies | Eq. (2) ensemble; Eq. (4) estimator |
| M | number of Monte Carlo moves (samples) | int | every attempt counts, accepted or not |
| $\alpha$ | maximum allowed particle displacement (step size) | length | tuned so ~half of moves are rejected; set to $(d-d_0)$ in runs |
| $\xi_1, \xi_2$ | random numbers for x,y displacement | float | uniform in $[-1, 1]$ |
| $\xi_3$ | acceptance random number | float | uniform in $[0, 1]$; accept if $\xi_3 < \exp(-\Delta E/kT)$ |
| $\Delta E$ | energy change of a trial move | energy | $E_{\text{new}} - E_{\text{old}}$; for rigid spheres 0 or +∞ |
| $\nu_r$ | number of ensemble systems in state r | int | target: $\propto \exp(-E_r/kT)$ |
| $E_r$ | energy of state r | energy | a state = a point in configuration space |
| $P_{rs}$ | a priori transition probability r->s | float | symmetric: $P_{rs}=P_{sr}$ (before accept test) |
| P | pressure | force/length (2D) | virial theorem, Eq. (7) |
| A | area (2D analog of volume) | area | fixed = 1 in the actual calculation |
| $A_0$ | close-packed area | area | reference for reduced density $A_0/A$ |
| $E_{\text{kin}}$ | total kinetic energy | energy | $= Nm\bar{v}^2/2 = NkT$ |
| m | particle mass | mass | Maxwellian velocity distribution |
| $\bar{v}^2$ | mean-square speed | (length/time)^2 | isotropic |
| $d_0$ | collision diameter / forbidden distance | length | twice actual sphere radius; Eq. (11a) |
| d | initial trigonal lattice nearest-neighbor spacing | length | $d = 1/14$ (unit square, 14 per row) |
| $\bar{n}$ | average contact number density at particle surface | 1/area | drives EOS, Eq. (10) |
| $\nu$ | auxiliary parameter | float, 0..7 | sets $d_0$ and $A/A_0$; NOT $\nu_r$ (state count) |
| K | outer-radius multiplier for RDF histogram | float >1 | chosen per $\nu$ |
| $N_m$ | number of pairs in histogram bin m | int | m = 1..64 equal-area zones |
| $N(r^2)$ | radial distribution function | count | binned in $\pi r^2$ |
| $N_{\frac{1}{2}}$ | RDF extrapolated to contact $r^2=d_0^2$ | count | $\propto \bar{n}$; differs by const factor of N, K |
| $\Delta A^2$ | histogram zone area | area | Eq. deltaA2 |
| $C_1..C_4$ | virial coefficients | dimensionless | Eq. (13); final numeric values in Eq. (14) |
| $A_{i,k}$ | cluster integral, i particles, k overlap bonds | dimensionless | normalized configuration-space volume; primes = distinct bond topologies |
| $f(r_{ij})$ | overlap bond indicator | {0,1} | 1 if $r_{ij}<d$ else 0 |

## Conventions

- **Reduced/derived units:** the actual rigid-sphere calculation fixes area $A=1$
  (unit square) and lattice spacing $d = 1/14$; the density is varied by changing
  $d_0$ via $\nu$, not by changing the box.
- **Boltzmann factor** written directly as $\exp(-E/kT)$; the paper does not use
  $\beta = 1/kT$.
- **Symmetric proposal:** the uniform-square proposal makes $P_{rs}=P_{sr}$, which
  is what makes the simple $\min$-style acceptance sample the Boltzmann
  distribution. A modern (Metropolis-Hastings) implementation with an asymmetric
  proposal would need the proposal-ratio correction, absent here.
- **Rejected moves are counted:** the current configuration is re-sampled on
  rejection (critical for correct averages).
- **Rigid spheres:** $\Delta E \in \{0, +\infty\}$, so the acceptance test reduces
  to "reject any move that creates an overlap."
- **Naming collision:** $\nu$ (auxiliary density parameter, 0..7) is distinct from
  $\nu_r$ (ensemble population of state r).
