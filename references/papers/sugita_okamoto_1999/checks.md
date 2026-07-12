# Checks - Sugita & Okamoto 1999 (Replica-Exchange MD)

Numeric fixtures for a port. The test system is the penta-peptide Met-enkephalin
(Tyr-Gly-Gly-Phe-Met, N-terminus acetyl-blocked, C-terminus N-methyl-blocked) in
gas phase with the all-atom AMBER force field, dielectric constant 1.

## Simulation parameters (test setup)

- Given: unit MD timestep = 0.5 fs.
- Given: run length = $2\times10^6$ timesteps = 1.0 ns per replica.
- Given: thermalization = 100 ps regular canonical MD per temperature, then 100 ps
  replica-exchange, before data collection.
- Given: $M = 8$ replicas/temperatures.
- Given: temperatures (exponentially spaced): 700, 585, 489, 409, 342, 286, 239,
  200 K.
- Given: replica exchange attempted every 10 fs (= every 20 timesteps).
- Given: $N_{\text{sim}} = 10^5$ measurements per replica (data stored just before
  each exchange).
- Given: per exchange event, 4 pairs of neighboring temperatures swapped; pairing
  alternated between the two possible choices each event.
- Reweighting: $R = M = 8$; $g_m = \text{const}$ assumed; $n_m = N_{\text{sim}} =
  10^5$. Self-consistent WHAM (eq:21, eq:22) converged in 10-100 iterations.

## Acceptance ratios (Table 1) - primary regression fixture

Given the above setup, expect the following replica-exchange acceptance ratios for
neighboring temperature pairs:

| pair (K) | expected acceptance ratio |
|---|---|
| 200 <-> 239 | 0.160 |
| 239 <-> 286 | 0.149 |
| 286 <-> 342 | 0.143 |
| 342 <-> 409 | 0.139 |
| 409 <-> 489 | 0.142 |
| 489 <-> 585 | 0.146 |
| 585 <-> 700 | 0.146 |

- Sanity criteria stated by the paper: ratios should be roughly uniform (all about
  15%) and each greater than 0.1 for the temperature ladder to be adequate.

## Energy fixture

- Given $T = 200$ K, comparing the highest-probability conformation:
  - regular canonical MD: average potential energy about -141 kcal/mol.
  - replica-exchange MD: average potential energy about -143 kcal/mol
    (about 2 kcal/mol lower than the trapped canonical result).
- At $T = 700$ K, canonical and replica-exchange potential-energy distributions
  agree (high-temperature canonical sampling is accurate).

## Acceptance-criterion unit test (from equations)

- Given two neighboring replicas $i,j$ with potential energies $E_i, E_j$ and
  inverse temperatures $\beta_m,\beta_n$, expect
  $\Delta = (\beta_n-\beta_m)(E_i - E_j)$ and accept with probability
  $\min(1, e^{-\Delta})$ (eq:15, eq:17).
- Degenerate check: if $E_i = E_j$ then $\Delta = 0$ and the swap is always
  accepted. If $\beta_m = \beta_n$ then $\Delta = 0$.
- Momentum check: after an accepted swap moving replica $i$ from $T_m$ to $T_n$,
  its velocities are multiplied by $\sqrt{T_n/T_m}$ (eq:12); this preserves
  $\langle K\rangle = \tfrac{3}{2}Nk_BT$ (eq:4) at the new temperature.
