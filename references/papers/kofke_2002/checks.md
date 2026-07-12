# Checks - Kofke 2002

Numeric fixtures an implementation can be tested against. All simulation values
are in reduced Lennard-Jones units ($\sigma = \epsilon = 1$).

## Simulation setup (LJ liquid, replica-exchange MC)

- Given: LJ liquid, number density `rho = 0.80`; system sizes
  `N in {32, 64, 108}` particles.
- Temperatures scanned: `T` from `1.0` to `2.0` in increments of `0.2`.
- Only one replica pair per simulation. Temperature pairs studied
  (`T0 <-> T1`):
  `1.0<->1.2, 1.0<->1.4, 1.0<->1.6, 1.0<->1.8, 1.0<->2.0, 1.2<->1.6,
   1.4<->1.6, 1.6<->2.0, 1.2<->2.0, 1.4<->2.0, 1.8<->2.0`.
- Trial mix: each MC trial is either an atom displacement or a replica
  exchange, chosen at random with displacements attempted 100x more often than
  exchanges.
- Trials per simulation: between `5e5` and `2e7`.
- Measured heat capacity at `T = 1.5`: `C_V/(N k) = 0.78`
  (i.e. per-molecule `c_V/k = 0.78`; extensive `C = C_V/k = 0.78 N`).

## Formula self-consistency fixtures

Let `B = beta1/beta0 = T0/T1` and `C = C_V/k` (extensive, `= 0.78 * N` at the
studied state).

1. Entropy-model / constant-C consistency (eq:11):
   given `B`, `C`, expect `DeltaS/k = -C ln(B)` (positive since `B < 1`).

2. Entropy-model estimate (eq:3, from eq:2 + eq:11):
   given `T0 < T1`, `N`, `c_V/k`, expect
   `p_acc ~ (T0/T1)^(N c_V/k) = exp(-DeltaS/k)`.
   Example: `N = 108`, `c_V/k = 0.78`, `T0=1.0`, `T1=2.0`
   -> exponent `N c_V/k = 84.24`, `(0.5)^84.24 ~ 2.9e-26` (entropy-only bound;
   the true acceptance is much larger because of the eq:12 attenuation term).

3. Exact constant-C acceptance (eq:10): given `C`, numerically integrate
   `p_acc = 4 * Gamma(2C)/Gamma(C)^2 * (2C+1)/C * INT_0^B kappa^C/(1+kappa)^{2(C+1)} dkappa`.
   Limits: as `B -> 1` (equal temperatures), `p_acc -> 1`; as `B -> 0`
   (distributions non-overlapping), `p_acc -> 0`.

4. Asymptotic form (eq:12, erratum-corrected): given `B`, `C -> inf`, expect
   `p_acc ~ exp(-DeltaS/k)/(pi C)^{1/2} * [4/(1+B)^2]^{C+1} * (1+B)/(1-B)`.
   Sanity: must approach eq:10 for large `C` and `B` not too close to 1.
   Regression guard against the ORIGINAL (wrong) 2002 text: the denominator of
   the last ratio is `(1-B)`, NOT `(1-B)^2`.

5. Gaussian model (eq:13, erratum-corrected): given `B`, `C`, expect
   `p_acc = erfc[ (C/2)^{1/2} * (1-B)/(1+B^2)^{1/2} ]`.
   Limits: `B -> 1` gives `erfc(0) = 1`; large argument -> 0.
   Regression guard against the ORIGINAL (wrong) 2002 text: the argument has the
   `(1+B^2)^{1/2}` denominator; the uncorrected form was `erfc[(C/2)^{1/2}(1-B)]`.

## Qualitative validation targets (from Fig. 2)

- For fixed `N`, `p_acc` correlates strongly with `DeltaS` (decreasing).
- Acceptance decays with `N`, but slower than the pure `exp(-DeltaS/k)`
  (eq:2) prediction; the eq:12 bracketed term with exponent `C+1` (> 1)
  attenuates the decay.
- eq:10 (integral, N=108) and eq:12 (asymptotic) coincide almost everywhere;
  the asymptotic form breaks down only as `DeltaS -> 0` (`B -> 1`).

## Formula selection guidance (temperature-ladder design)

Which acceptance formula to use when spacing a replica-exchange temperature
ladder, as a function of neighbor spacing `B = T0/T1`:

- Wide spacing (`B` well below 1, large `DeltaS`): eq:12 asymptotic is accurate
  and cheapest.
- Tight spacing (`B -> 1`, small `DeltaS`): eq:12 breaks down. Use the exact
  integral eq:10, or the erratum-corrected Gaussian eq:13, which Fig. 2 shows
  tracks the simulation data down to small `DeltaS`.
- Do NOT use eq:12 to set closely-spaced neighbors near unity; it will
  mis-predict the acceptance there.
