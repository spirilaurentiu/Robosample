# Checks / fixtures - Pear & Weiner 1979

Concrete numbers and closed forms an implementation can regress against. Reduced units:
$m=l=1$, $T_{\mathrm{ref}}=600$ K $=1$, $\beta=1/(k_B T)$.

## Unit / time-step conversions

- Reference frequency: `omega = 1.33e13 sec^-1`.
- Time step `dt = 1` (reduced) `= 7.51e-14 sec` (cgs). Given the reduced-unit system
  ($m$ carbon mass, $l$ C-C bond, $T_{\mathrm{ref}}=600$ K all unity).
- Time steps used in simulations: `dt = 0.1/omega`, `0.05/omega`, `0.025/omega`.
  Consistency check: `0.1/omega = 0.1/1.33e13 s = 7.52e-15 s`.
- Barrier height relation: given stiffness `k`, expect `E_b = k/4` (from continuity of
  the piecewise-quadratic potential Eq. 4.1).

## Metric determinant, three-bond 90-degree chain (Eq. 3.6)

Closed form: `g(phi) = C_g * (a + b*cos(phi) + c*cos^2(phi) + cos^4(phi))` with, for
end/interior mass ratio `alpha`:

```
a = (1 + 7*alpha + 15*alpha^2 + 10*alpha^3 + 2*alpha^4) / alpha^4
b = 2*(1 + alpha) / alpha^2
c = -(3*alpha^2 + 7*alpha + 6) / alpha^2
```

Fixtures (evaluate coefficients):
- `alpha = 1`: `a = (1+7+15+10+2)/1 = 35`; `b = 2*2/1 = 4`; `c = -(3+7+6)/1 = -16`.
  So `g ∝ 35 + 4*cos(phi) - 16*cos^2(phi) + cos^4(phi)`.
- `alpha = 10`: `a = (1+70+1500+10000+20000)/10000 = 31571/10000 = 3.1571`;
  `b = 2*11/100 = 0.22`; `c = -(300+70+6)/100 = -376/100 = -3.76`.
  So `g ∝ 3.1571 + 0.22*cos(phi) - 3.76*cos^2(phi) + cos^4(phi)`.

Behavioural check: for `alpha=10` the rigid-model dihedral distribution
`rho_R(phi) = C_R * sqrt(g(phi))` is markedly nonuniform; the flexible model gives the
uniform `rho_F = 1/(2*pi)`; adding the Fixman potential `U = kB*T*ln(sqrt(g))` restores
`rho_RF = 1/(2*pi)` (Figs. 3a-c).

## End-to-end distance (three-bond 90-degree chain)

- `R(phi) = l*(3 + 2*cos(phi))^(1/2)`.
  - `phi = 0`  -> `R = l*sqrt(5) ≈ 2.2360679*l`.
  - `phi = pi` -> `R = l*sqrt(1) = l`.
  - `phi = pi/2` -> `R = l*sqrt(3) ≈ 1.7320508*l`.

## Transition-state rate formulas (regression targets)

- Uncompensated rigid (Eq. 4.10):
  `f_TS = (I_c / (2*sqrt(pi))) * sqrt(g(phi_b)*G33*kB*T) * exp(-E_b/(kB*T))`,
  `I_c = 1 / integral_0^{2pi} sqrt(g(phi)) * exp(-beta*V(phi)) dphi`.
  Predicts SPLITTING of rate between `phi_b=0` and `phi_b=pi` via `g(phi_b)`.
- Fixman-compensated (Eq. 4.11):
  `f_TS = (I_c / (2*sqrt(pi))) * sqrt(kB*T*G33) * exp(-E_b/(kB*T))`,
  `I_c = 1 / integral_0^{2pi} exp(-V(phi)/(kB*T)) dphi`.
  `g(phi_b)` drops out -> rate INDEPENDENT of barrier position (curves for `phi_b=0` and
  `phi_b=pi` coalesce); note `G33` has the same value at `phi_3=0` and `phi_3=pi`.
- Applied stress (Eq. 5.9):
  `f_TS(sigma) = (I_c(sigma)/(2*sqrt(pi))) * sqrt(kB*T*G33) * sinh(sigma*R(phi_b)/(kB*T)) * exp(-E_b/(kB*T)) / R(phi_b)`,
  `I_c(sigma) = 1 / integral_0^{2pi} (1/R(phi)) * sinh(sigma*R(phi)/(kB*T)) * exp(-V(phi)/(kB*T)) dphi`.
  Limit check: as `sigma -> 0`, `sinh(x)/R -> sigma/(kB*T)` and Eq. (5.9) must reduce to
  Eq. (4.11).

## Simulation parameter sets (from figure captions)

| Fig | model | alpha | dt | eta/(m*omega) | sigma*l/E_b | phi_b | note |
|---|---|---|---|---|---|---|---|
| 5 | rigid, no Fixman | 10 | 0.1/omega and 0.05/omega | 2 and 0.5 | 0 | {0, pi} | dt=0.1,eta=2 poor; dt=0.05,eta=0.5 good vs Eq. 4.10 |
| 6 | rigid + Fixman | 10 | 0.05/omega | 0.5 | 0 | {0, pi} | curves coalesce vs Eq. 4.11 |
| 7 | rigid, no Fixman (effective) | 10 | 0.05/omega | 2 | 0 | {0, pi} | effective < transition-state |
| 8 | rigid + Fixman (effective) | 10 | 0.1/omega | 2 | 0 | {0, pi} | splitting disappears |
| 9 | flexible (effective) | 10 | 0.05/omega | 2 | 0 | {0, pi} | rate independent of phi_b |
| 10 | rigid + Fixman, stress | 1 | 0.05/omega | 0.25 | 0.5 | {0, pi} | lower curve phi_b=pi, upper phi_b=0 (Eq. 5.9) |
| 11 | rigid + Fixman, stress | 1 | 0.025/omega (2 lowest-T pts) | 0.25 | 0.25 | {0, pi} | reduced splitting |

## Qualitative invariants (assertions)

- Transition-state rate overestimates effective rate; gap grows with viscosity `eta`.
- Flexible-model equilibrium dihedral distribution is uniform even as spring constants
  `k_l, k_theta -> infinity` (limit and thermal average do not commute).
- Arrhenius plots under applied stress are CURVED, with opposite curvature direction for
  `phi_b=0` vs `phi_b=pi`; splitting increases at low temperature and with larger
  `sigma*l/E_b` (0.5 splits more than 0.25).
