# Checks - Understanding the Hastings Algorithm (Minh & Minh, 2014)

This is a theory/unification paper: no numeric benchmark tables. The concrete,
testable content is a set of exact ALGEBRAIC IDENTITIES between acceptance rules.
These are excellent regression fixtures for an MCMC accept/reject implementation.
Use any strictly-positive `p(x), p(y), gamma(x|y), gamma(y|x)` unless noted.

## Reduction / special-case fixtures

- Symmetric proposal `gamma(x|y) == gamma(y|x)`:
  given this, expect `alpha_MH(x,y) == min{p(y)/p(x), 1}` (eq:3 collapses to Metropolis eq:MT).

- MH == Algorithm M middle branch:
  given `L(x,y) < k(x,y) < H(x,y)`, expect `alpha_M(x,y) == alpha_MH(x,y) == min{ [gamma(x|y) p(y)] / [p(x) gamma(y|x)], 1 }` (eq:9, eq:10 middle).

- MAR with `C(x,y)=1`:
  given `C=1`, expect `alpha_MA(x,y) == alpha_MH(x,y)` (eq:20 => eq:3).

- MIR with `C(x,y)=1`:
  given `C=1`, expect `alpha_MI(x,y) == alpha_MH(x,y)` (eq:27 => eq:3). MAR and MIR give the SAME acceptance value for the same `C`.

- Barker from Algorithm M (majorizing side):
  given `k(x,y) = p(x)/gamma(x|y) + p(y)/gamma(y|x)` (which is `>= H`),
  expect `alpha_M(x,y) == alpha_BK(x,y) == (1 + [p(x) gamma(y|x)] / [gamma(x|y) p(y)])^{-1}`.

- Barker from Algorithm M (minorizing side):
  given `k(x,y) = ( gamma(x|y)/p(x) + gamma(y|x)/p(y) )^{-1}` (which is `<= L`),
  expect `alpha_M(x,y) == alpha_BK(x,y)` (same Barker value).

- Barker from MAR:
  given `C(x,y) = ( p(y)/gamma(y|x) + p(x)/gamma(x|y) ) / max{p(x)/gamma(x|y), p(y)/gamma(y|x)}` (`> 1`),
  expect `alpha_MA(x,y) == alpha_BK(x,y)` (eq:21).

## Ordering / bound fixtures

- Optimality (Peskun): for the same proposal, expect `alpha_MH(x,y) >= alpha_V(x,y)` for every Hastings variant V (max over all `C>=1`, achieved at `C=1`). In particular `alpha_MH >= alpha_BK`.
- Range: every `alpha_*(x,y)` in `[0,1]`.
- Bracket ordering: always `L(x,y) <= H(x,y)`, with `L=H` iff `p(x)/gamma(x|y) == p(y)/gamma(y|x)`.

## Detailed-balance invariant (the master test)

For every variant's transition kernel `P(y|x)` and all `x != y`:
- expect `p(x) P(y|x) == p(y) P(x|y)` (reversibility w.r.t. `p`, hence `pi` invariant).

Numeric example fixture (pick concrete values to unit-test the MH rule eq:3):
- given `p(x)=1, p(y)=2, gamma(y|x)=0.5, gamma(x|y)=0.5`,
  expect `alpha_MH = min{ (0.5*2)/(1*0.5), 1 } = min{2,1} = 1`.
- given `p(x)=2, p(y)=1, gamma(y|x)=0.5, gamma(x|y)=0.5`,
  expect `alpha_MH = min{ (0.5*1)/(2*0.5), 1 } = min{0.5,1} = 0.5`.
- same second case, Barker: `alpha_BK = (1 + (2*0.5)/(0.5*1))^{-1} = (1+2)^{-1} = 1/3 ~= 0.3333`
  (confirms `alpha_MH = 0.5 > alpha_BK = 0.333`).

## Duplication-probability fixture (Algorithm L, two-stage)

- given any `k(x,y)>0`, expect
  `P(x duplicated) = 1 - min{p(y)/(k*gamma(y|x)),1} * min{k*gamma(x|y)/p(x),1} = 1 - alpha_M(x,y)`.
  (Splitting the rejection into type-x + type-y must sum to the same `1 - alpha_M`.)
