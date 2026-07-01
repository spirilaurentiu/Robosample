# Checks / Fixtures - Sohl-Dickstein 2014, LAHMC

## Shared experimental hyperparameters (Section 5)
- Leapfrog step length: ε = 1
- Leapfrog steps per L application: M = 10
- LAHMC max leapfrog applications: K = 4
- β ∈ {1, 0.1} (as stated per condition)
- Ill-conditioned Gaussians (Fig 2): covariance eigenvalues log-linearly spaced
  in [1, 10^6]; grid-search figure (Fig 3) uses eigenvalues {1, 10^5}.
- Rough-well energy (Eq 34): σ1 = 100, σ2 = 2.
- HMC optimal acceptance rate reference: ~65% (Neal 2010); HMC here is close to it.

## Claim
- LAHMC outperforms standard HMC for all explored hyperparameters, often by
  more than a factor of 2 in mixing time (gradient/function evaluations to reach
  autocorrelation 0.5).

## Table 1 - fraction of transitions to each target state
Columns: F ζ, L ζ, L^2 ζ, L^3 ζ, L^4 ζ. Each row is one (distribution, sampler)
condition. `given <distribution, sampler, β>, expect <fractions>`.

| Distribution | Sampler | β | Fζ | Lζ | L^2ζ | L^3ζ | L^4ζ |
|---|---|---|---|---|---|---|---|
| 2d Gaussian | HMC | 1   | 0.079 | 0.921 | 0     | 0     | 0     |
| 2d Gaussian | LAHMC | 1 | 0.000 | 0.921 | 0.035 | 0.044 | 0.000 |
| 2d Gaussian | HMC | 0.1 | 0.080 | 0.920 | 0     | 0     | 0     |
| 2d Gaussian | LAHMC | 0.1 | 0.000 | 0.921 | 0.035 | 0.044 | 0.000 |
| 100d Gaussian | HMC | 1   | 0.147 | 0.853 | 0     | 0     | 0     |
| 100d Gaussian | LAHMC | 1 | 0.047 | 0.852 | 0.059 | 0.035 | 0.006 |
| 100d Gaussian | HMC | 0.1 | 0.147 | 0.853 | 0     | 0     | 0     |
| 100d Gaussian | LAHMC | 0.1 | 0.047 | 0.852 | 0.059 | 0.035 | 0.006 |
| 2d Rough Well | HMC | 1   | 0.446 | 0.554 | 0     | 0     | 0     |
| 2d Rough Well | LAHMC | 1 | 0.292 | 0.554 | 0.099 | 0.036 | 0.019 |
| 2d Rough Well | HMC | 0.1 | 0.446 | 0.554 | 0     | 0     | 0     |
| 2d Rough Well | LAHMC | 0.1 | 0.292 | 0.554 | 0.100 | 0.036 | 0.019 |

Notes for regression tests:
- Standard HMC only ever transitions to Fζ or Lζ (all L^{>=2} columns are 0).
- LAHMC drives the Fζ (momentum-flip) fraction toward 0: exactly 0.000 for both
  2d Gaussian conditions; substantially reduced (0.147 -> 0.047, 0.446 -> 0.292)
  for 100d Gaussian and Rough Well.
- The Lζ fraction is essentially unchanged between HMC and LAHMC for a given
  distribution (e.g. 0.921, 0.853/0.852, 0.554), i.e. LAHMC redistributes the
  former flip mass onto L^2..L^4, not onto L.
- Row sums equal 1 (up to rounding).

## Unit / invariant checks for an implementation
- H(Fζ) = H(ζ) exactly (momentum flip conserves energy).
- L^{-1}ζ = FLFζ exactly (leapfrog reversibility) up to float error.
- |det ∂(Fζ)/∂ζ| = |det ∂(Lζ)/∂ζ| = 1 (volume preservation).
- Generalized detailed balance: p(ζ) π_{L^a}(ζ) = p(FL^a ζ) π_{L^a}(FL^a ζ).
- Σ_a π_{L^a}(ζ) + π_F(ζ) = 1 (total outgoing probability normalized).
- p(ζ')/p(ζ) = exp(H(ζ) - H(ζ')).
