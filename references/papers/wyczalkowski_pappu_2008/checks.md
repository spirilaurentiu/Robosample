# Checks / fixtures - Wyczalkowski & Pappu (2008)

## Acetamide free energy of hydration (Table I(a), current work)

BAR estimator, bootstrap statistical errors, 21 replicas.

| condition | $\Delta F$ (kcal/mol) | $\epsilon_{rms}$ (kcal/mol) |
|---|---|---|
| No replica exchange | $-8.35 \pm 0.051$ | 0.120 |
| Replica exchange | $-8.14 \pm 0.053$ | 0.023 |

- Fixture: given the same 21-replica multicanonical protocol, expect replica
  exchange to lower $\epsilon_{rms}$ by roughly an order of magnitude
  ($0.120 \to 0.023$, factor $\approx 5$) while leaving the bootstrap
  statistical error on $\Delta F$ essentially unchanged ($0.051 \to 0.053$).
- Fixture: RE reaches the same $\epsilon_{rms}$ magnitude about 5x faster
  (stated range 4-8x shorter simulations).

## Literature / experimental comparison (Table I(b))

All computational values use OPLS-AA for acetamide.

| source | $\Delta G$ (kcal/mol) | water model / estimator |
|---|---|---|
| MacCallum & Tieleman | $-8.25 \pm 0.26$ | TIP4P, TI |
| Shirts et al. | $-8.20 \pm 0.03$ (no long-range vdW correction) | TIP3P, TI |
| Chang et al. | $-8.54 \pm 0.1{-}0.3$ | TIP4P, BAR |
| Udier-Blagovic et al. | $-9.65 \pm 0.3{-}0.5$ | TIP4P, FEP |
| Experimental | $-9.54$ | - |

- Fixture: this work's $\Delta F$ ($-8.14$ to $-8.35$) should fall within the
  cluster of TIP4P/TIP3P computational values ($-8.2$ to $-8.5$), not the
  experimental $-9.54$.

## Analytic limits of the linearized swap probability (eq:18)

- Given $\delta_\lambda = 0$: expect $\langle p_{swap} \rangle = 1/2$.
- Given small $\delta_\lambda$: expect $\langle p_{swap} \rangle = \tfrac12 -
  \tfrac14 \beta^2 \delta_\lambda^2 C_\lambda$, monotonically decreasing in
  $C_\lambda$ and in $\delta_\lambda^2$.
- Consistency: a positive spike in $C_\lambda$ coincides with a downward spike
  in $\langle p_{swap} \rangle$ and a large $\epsilon_H$ (near $\lambda \sim
  0.1{-}0.3$, cavitation region peaks near $\lambda \sim 0.15$).

## Fluctuation-theorem / FEP identities (exact, for unit tests)

- eq:19: given a converged $U_0$ ensemble, $\langle \exp(-\beta W_D^F)
  \rangle_0 = 1$ (Jarzynski identity form). Expect deviation from 1 to scale
  with sampling of the $W_D^F < 0$ tail.
- eq:4/eq:11: given any $\Gamma_0, \Gamma_1$, check
  $\Delta U_{swap} = U_0(\Gamma_1)+U_1(\Gamma_0)-U_0(\Gamma_0)-U_1(\Gamma_1)
  = W^F(\Gamma_0)+W^R(\Gamma_1)$.
- eq:7: $W_D^F(\Gamma) + W_D^R(\Gamma) = 0$ for the same configuration.
- eq:12b: $\rho_N'(\gamma) = \rho_N(\gamma')$ (swap = pair reflection).

## Simulation setup parameters (reproduction fixtures)

| parameter | value |
|---|---|
| replicas | 21, $\lambda \in [0,1]$ step 0.05 |
| solute / solvent | 1 rigid acetamide (ACE), 343 TIP4P waters |
| force field | OPLS-AA (solute), TIP4P (water) |
| ensemble | N-V-T, 298 K, 21.8 Angstrom cubic box |
| sampler | Metropolis MC (translation+rotation of one water per move) |
| cycle | 343 MC moves |
| production | $2\times10^6$ cycles, first $10^5$ discarded |
| avg acceptance | 31% |
| electrostatic cutoff | 10.5 Angstrom (group based) |
| vdW cutoff | 10 Angstrom (atom based) |
| long-range correction | none |
| RE sim-round length | Normal(mean 500, sd 50) cycles |
| energy autocorrelation time | ~500 cycles |
| swap round | $21^2$ attempts, random replica pairs (non-neighbor allowed) |
| logging | native/foreign energies and $dU/d\lambda_C, dU/d\lambda_{LJ}$ every 10 cycles |
| bootstrap | 10000 resamples; $n^* = 1900$ (one independent obs per 2 autocorr times) |

## Soft-core potential parameters (Appendix B)

| parameter | value | used in |
|---|---|---|
| $\alpha_C$ | 1.5 Angstrom | eq:B1 (Coulomb soft core) |
| $\alpha_{LJ}$ | 0.5 Angstrom | eq:B3a (LJ soft core) |
| $b$ (= $a$) | 4 | eq:B3a exponent |
| $k$ | 1 | eq:B3b well-depth scaling |

- Boundary check (eq:B3b): $B(\lambda_{LJ}=0) = 0$ and $B(\lambda_{LJ}=1) =
  4\epsilon$.
- Boundary check (eq:B1): $U_C(\lambda_C=0) = 0$; $U_C(\lambda_C=1) =
  q_i q_j / r$ (unscaled Coulomb).
