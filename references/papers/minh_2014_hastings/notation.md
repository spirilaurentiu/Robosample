# Notation - Understanding the Hastings Algorithm (Minh & Minh, 2014)

Reduced/abstract sampling notation; no physical units. State space `(E, E)`.

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| `pi(.)` | normalized target distribution (invariant dist. of the chain) | probability density | `pi = p/P` |
| `p(.)` | UN-normalized target density | non-negative real | `pi(x) = p(x)/P` |
| `P` | normalizing constant, `P = integral_E p(x) dx` | positive real | usually unknown |
| `x` | current chain state `X_n` | point in `E` | `x ~ pi` |
| `y` | proposed next state | point in `E` | drawn from `gamma(.|x)` |
| `z` | dummy integration variable over `E` | point in `E` | - |
| `gamma(y|x)` | proposal density of `y` given current `x` | density | may or may not depend on `x`; `gamma(.)` if independent |
| `s(x,y)` | free symmetric function parameterizing Hastings acceptance | positive real | `s(x,y)=s(y,x)`, constrained by eq:1 to keep alpha in [0,1] |
| `k(x,y)` | free symmetric function parameterizing Algorithm M | positive real | `k(x,y)=k(y,x) > 0`; constant case excluded |
| `delta(x,y)` | symmetric function parameterizing the Stein algorithm | positive real | `delta(x,y)=delta(y,x)` |
| `C(x,y)` | symmetric multiplier for majorizing/minorizing coeff | real `>= 1` | `C=1` => Metropolis-Hastings (optimal) |
| `L(x,y)` | lower bracket `min{p(x)/gamma(x|y), p(y)/gamma(y|x)}` | positive real | symmetric |
| `H(x,y)` | upper bracket `max{p(x)/gamma(x|y), p(y)/gamma(y|x)}` | positive real | symmetric |
| `M` | absolute majorizing coefficient (AR/IMAR) | positive real | `M gamma(z) >= p(z)` for all `z`; "deficient" if only for some `z` |
| `M(x,y)` | RELATIVE majorizing coefficient (MAR) | positive real | symmetric, `>= H(x,y)` |
| `m` | absolute minorizing coefficient (IMIR) | positive real | `m gamma(z) <= p(z)` for all `z` |
| `m(x,y)` | RELATIVE minorizing coefficient (MIR) | positive real | symmetric, `<= L(x,y)` |
| `alpha_*(x,y)` | acceptance probability of variant `*` | in `[0,1]` | accept `y` if `r <= alpha` |
| `r, r1, r2` | uniform random deviates | `~ U(0,1)` | independent per draw |
| `I(a)` | indicator, 1 if `a` true else 0 | {0,1} | used in transition kernel self-loop |
| `X_n` | n-th chain sample | point in `E` | `n = 1,2,3,...` |

## Algorithm-variant acronyms

| acronym | name | proposal | coefficient regime |
|---|---|---|---|
| HA | Hastings algorithm | `gamma(.|x)` | general, via `s` |
| MH | Metropolis-Hastings | `gamma(.|x)` | `C=1` (optimal, max acceptance) |
| MT | Metropolis (1953) | symmetric `gamma` | - |
| BK | Barker (1965) | `gamma(.|x)` | `s=1` / `C>1` |
| ST | Stein | `gamma(.|x)` | via `delta` |
| M | Algorithm M | `gamma(.|x)` | any `k(x,y)>0` |
| AR | Acceptance-Rejection (i.i.d.) | `gamma(.)` | absolute `M` |
| IMAR | Independence Markovian Accept-Reject | `gamma(.)` | absolute `M` |
| MAR | Markovian Accept-Reject | `gamma(.|x)` | `k = M(x,y) >= H` (type-y duplications) |
| IMIR | Independence Markovian Minorizing | `gamma(.)` | absolute `m` |
| MIR | Markovian Minorizing | `gamma(.|x)` | `k = m(x,y) <= L` (type-x duplications) |
