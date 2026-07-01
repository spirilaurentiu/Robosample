# Understanding the Hastings Algorithm

David D. L. Minh (Illinois Institute of Technology) and Do Le (Paul) Minh
(California State University, Fullerton). arXiv:1408.4438v1, 2014.

## Abstract

The Hastings algorithm is a key tool in computational science. While mathematically
justified by detailed balance, it can be conceptually difficult to grasp. This paper
presents two complementary and intuitive ways to derive and understand the algorithm.
In this framework it is straightforward to see that the Metropolis-Hastings algorithm
has the highest acceptance probability of all Hastings algorithms.

Keywords: Hastings algorithm; Metropolis-Hastings; Markov Chain Monte Carlo; simulation.
MSC: Primary 65C05; Secondary 78M31, 80M31.

## 1 Introduction

### 1.1 The Hastings algorithm (HA)

The Hastings algorithm (HA) is an MCMC method that, rather than generating i.i.d.
samples from a target `pi(.)` on state space `(E, E)`, generates a Markov chain
`{X_n}` with `pi(.)` as its invariant distribution. Chain variates are dependent but
may still be used to estimate expectations under `pi(.)`. (The same symbol denotes both
a measure and its density.)

Often `pi(.) = p(.)/P` where the normalizing constant `P = integral_E p(x) dx` is
unknown; `p(.)` is the un-normalized target. Writing `x ~ pi(.)` and `x ~ p(.)`
interchangeably. `U(0,1)` is the uniform distribution on `(0,1)`. A proposal density
`gamma(.|x)` (possibly `x`-dependent) whose variates can be generated is required.

**Algorithm HA (Hastings).** Given `X_n = x ~ pi(.)`, generate `X_{n+1} ~ pi(.)` by:
1. HA1. generate `y ~ gamma(.|x)` and `r ~ U(0,1)`.
2. HA2. if `r <= alpha_HA(x,y)`, output `X_{n+1} = y`.
3. HA3. else output `X_{n+1} = x`.

Here `alpha_HA` (see eq:1 in equations.md) is defined via a symmetric function
`s(x,y)` constrained so `alpha_HA in [0,1]`. (Hastings 1970 stated this condition in
terms of the normalized `pi(.)` rather than `p(.)`.)

### 1.2 Special forms of HA

- **Metropolis-Hastings (MH):** with symmetric proposals, the Metropolis (1953)
  acceptance is `alpha_MT = min{p(y)/p(x), 1}` (eq:MT). Hastings generalized this by
  choosing `s = s_MH` (eq:2), giving the MH acceptance `alpha_MH` (eq:3).
- **Barker (BK):** Barker (1965) used `alpha_BK^(s) = (1 + p(x)/p(y))^{-1}` for
  symmetric proposals; Hastings generalized it by `s(x,y)=1`, giving `alpha_BK` (eq:4).
- **Another special form:** choosing `s` as in eq:5 gives the product-of-mins
  acceptance eq:6.

### 1.3 Detailed balance

To show a chain has invariant `pi(.)` it suffices that its transition kernel `P(.|.)`
satisfies detailed balance (reversibility) w.r.t. `p(.) = P pi(.)`:
`p(x) P(y|x) = P(x|y) p(y)` for all `x,y`. All kernels here split as
`P(y|x) = r1(y|x) + I(x=y) r2(y|x)`; since the self-loop part balances trivially,
only `x != y` is treated.

### Derivation (not implemented): HA satisfies detailed balance

For HA the transition kernel is (eq:7)
`P_HA(y|x) = alpha_HA(x,y) gamma(y|x) + I(x=y) integral_E (1 - alpha_HA(x,z)) gamma(z|x) dz`.
Substituting eq:1 and using the symmetry of `s`,
`p(x) P_HA(y|x) = p(x) s(x,y) [p(y) gamma(x|y)] / [p(x) gamma(y|x) + p(y) gamma(x|y)] gamma(y|x) = P_HA(x|y) p(y)`,
so detailed balance holds. Verifying HA works is simple; understanding *why* it works
motivates the two intuitive derivations below (Algorithm M; and the
majorizing/minorizing view MAR/MIR).

## 2 Algorithm M

### 2.1 Algorithm M

**Algorithm M.** Given `X_n = x`, generate `y ~ gamma(.|x)`, `r ~ U(0,1)`; accept
`y` if `r <= alpha_M(x,y)` (eq:8), where `k(x,y) > 0` is a free symmetric function.
Its kernel mirrors eq:7 and satisfies detailed balance for any symmetric `k`.
If `k` is a positive constant `k`, then `p(.)/k` is just another un-normalized target
and eq:8 reduces to eq:6; that constant case is excluded hereafter.

Defining brackets `L(x,y)`, `H(x,y)` (see eq:L-H), Algorithm M's acceptance is the
piecewise eq:10. MH is the middle branch (`L < k < H`, eq:9). Barker arises two ways:
setting `k = p(x)/gamma(x|y) + p(y)/gamma(y|x) >= H`, or
`k = (gamma(x|y)/p(x) + gamma(y|x)/p(y))^{-1} <= L`; both give `alpha_BK` (eq:4).

### 2.2 Algorithm M and HA are equivalent

### Derivation (not implemented): equivalence

HA is a special case of M: for any Hastings `s`, define `M_s >= H` (eq:11) or
`m_s <= L` (eq:13); setting `k = M_s` (eq:12) or `k = m_s` (eq:14) in eq:10 reproduces
`alpha_HA`. Conversely M is a special case of HA, by cases on `k` vs `L,H` (using
`s = sMH` in the middle case). The map `s <-> k` is not one-to-one: the set of valid
`k>0` is larger than the set of `s` (which must also satisfy eq:1); every `s` yields at
least two `k` (namely `M_s` and `m_s`), and all `k in (L,H)` map to `sMH`.

### 2.3 Algorithm M and the Stein algorithm

Stein (in Liu 2001, p.112) proposed acceptance `alpha_ST = delta(x,y)/[p(x) gamma(y|x)]`
with symmetric `delta` (eq:15). By the same logic Algorithm M, Stein, and HA are all
equivalent (via `M_delta = p(x)p(y)/delta >= H` or `m_delta = delta/[gamma(x|y)gamma(y|x)] <= L`).
Previously the Stein-HA relationship was unclear; here they are unified. The point of
introducing Algorithm M is that it can be built up *intuitively* from
Acceptance-Rejection, rather than merely accepted because it satisfies detailed balance.

## 3 Markovian Acceptance-Rejection (MAR)

### 3.1 Acceptance-Rejection (AR)

AR (von Neumann 1951) draws i.i.d. samples from `p(.)` using proposal `gamma(.)` and a
majorizing coefficient `M` with `M gamma(z) >= p(z)` for all `z`. Loop: draw
`y ~ gamma(.)`, `r ~ U(0,1)`; accept (`reject=0`, output `y`) if `r <= p(y)/[M gamma(y)]`
(eq:AR), else repeat. The pairs `(y, r M gamma(y))` are uniform under `M gamma(.)`;
accepted ones are uniform under `p(.)`, hence `~ pi(.)`.

### 3.2 Independence Markovian Acceptance-Rejection (IMAR)

IMAR modifies AR: when `y` is rejected, the current `x` is REPEATED (rather than
re-drawing). Accept `y` if `r <= alpha_IMA = p(y)/[M gamma(y)]` (eq:IMAR), else output
`x`. This makes samples dependent (a Markov chain) but preserves detailed balance:
`p(x) P_IMA(y|x) = p(x) [p(y)/(M gamma(y))] gamma(y) = P_IMA(x|y) p(y)`.
Expected duplications per delivered `z` is `M-1`, so total deliveries proportional to
`M p(z)`, i.e. to `pi(z)`. With a *deficient* `M` (only `M gamma(z) >= p(z)` for some
`z`), AR/IMAR sample `min{p(.), M gamma(.)}`, giving `alpha_D = min{p(y)/(M gamma(y)), 1}`;
they still deliver `y ~ pi(.)` within `{z : M gamma(z) >= p(z)}`.

### 3.3 Markovian Acceptance-Rejection (MAR)

When `gamma(.|x)` depends on `x`, an absolute `M` is too restrictive (a pair with
`gamma(eta|xi)=0`, `p(eta)>0` forces `M = infinity`). Because MAR delivers either `x`
or `y` each step, only a symmetric RELATIVE majorizing coefficient `M(x,y)` is needed
with `M(x,y) gamma(x|y) >= p(x)` and `M(x,y) gamma(y|x) >= p(y)`, i.e. `M(x,y) >= H(x,y)`
(eq:16). Write `M(x,y) = C(x,y) H(x,y)` with symmetric `C >= 1` (eq:17).

**Algorithm MAR.** Accept `y` if `r <= alpha_MA(x,y)` (eq:18 = eq:19 = eq:20). It
satisfies detailed balance for symmetric `M`. Both `x,y` lie in the region
`{z : M(x,z) gamma(z|x) >= p(z)}`, so a relative `M` deficient as an absolute
coefficient still suffices.

### 3.4 BK in MAR; 3.5 MH in MAR

Choosing `C` as in eq:21 (`> 1`) makes MAR reproduce Barker (eq:4). Setting `C = 1`
makes `alpha_MA = alpha_MH` (eq:3), so MAR becomes MH. Peskun (1973) proved (via
partial ordering of kernels) that `alpha_MH` is optimal - minimal asymptotic variance
of sample-path averages. In the MAR view this is immediate: eq:20 is maximized at
`C = 1`. Thus MH has the highest acceptance probability of all Hastings algorithms.
Intuitively, the most efficient majorizer "touches" `p(.)` at one point; `C=1` makes
`M(x,y) gamma(x|y) = p(x)` or `M(x,y) gamma(y|x) = p(y)`. Any `C>1` (as in BK)
needlessly rejects proposals.

### 3.6 MAR and Algorithm M

### Derivation (not implemented)

MAR is Algorithm M with `k = M(x,y) >= H`. Conversely, for any `k>0`,
`M_k(x,y) = k max{p(x)/(k gamma(x|y)),1} max{p(y)/(k gamma(y|x)),1}` is a relative
majorizing coefficient (`>= H`, checked by three cases on `k` vs `L,H`), and yields
`alpha_MA = alpha_M`. Hence M and MAR are equivalent.

### 3.7 MAR and HA

Since MAR = M = HA. Directly: for `M(x,y) >= H`, set `s` as in eq:24 to recover
`alpha_MA` in HA; conversely set `M = M_s >= H`. There is a one-to-one map between
symmetric `s` satisfying eq:1 and symmetric `M` of the form eq:17 - but `M(.,.)` has an
intuitive meaning (relative majorizing coefficient) whereas `s(.,.)` is opaque.

## 4 Markovian Minorizing (MIR)

### 4.1 Independence Markovian Minorizing (IMIR)

Dual to IMAR. With state-independent `gamma(.)` whose support contains that of `p(.)`,
and an absolute minorizing coefficient `m` with `m gamma(z) <= p(z)` for all `z`:
accept `y` if `r <= alpha_IMI = m gamma(x)/p(x)` (eq:IMIR; note it depends on the
CURRENT `x`), else output `x`. Kernel `P_IMI(y|x) = alpha_IMI(x,y) gamma(y)` satisfies
detailed balance. With a deficient `m`, IMIR samples `max{p(.), m gamma(.)}` giving
`alpha_d = min{m gamma(x)/p(x), 1}` and still delivers `y ~ pi(.)` within
`{z : m gamma(z) <= p(z)}`.

An equivalent, more intuitive rewrite (Algorithm IMJ): draw `r`; if
`r <= m gamma(x)/p(x)` then draw `y ~ gamma(.)` and output it, else output `x`
(no need to draw `y` on rejection). Each delivered `x` is duplicated until first
success of Bernoulli trials with success prob `m gamma(x)/p(x)`; expected duplications
`p(x)/[m gamma(x)]`, giving total deliveries proportional to `p(x)`. (Minh et al. 2012
used `m` to make MCMC regenerative.)

### 4.2 Markovian Minorizing (MIR)

For `x`-dependent proposals, an absolute `m` is too restrictive; use a symmetric
RELATIVE minorizing coefficient `m(x,y)` with `m(x,y) gamma(x|y) <= p(x)`, i.e.
`m(x,y) <= L(x,y)` (eq:25). Write `m(x,y) = L(x,y)/C(x,y)` with symmetric `C >= 1`
(eq:26). **Algorithm MIR:** accept `y` if `r <= alpha_MI(x,y)` (eq:27), else output `x`.
Kernel `P_MI(y|x) = alpha_MI(x,y) gamma(y|x)` satisfies detailed balance.

### 4.3 / 4.4 MIR and HA / Algorithm M

MIR's acceptance (eq:27) equals MAR's (eq:20), so MIR is equivalent to HA and to
Algorithm M (with `k = m(x,y) <= L`, type-x duplications only). The converse uses
`m_k(x,y) = k min{p(x)/(k gamma(x|y)),1} min{p(y)/(k gamma(y|x)),1} <= L` (three cases),
giving `alpha_MI = alpha_M`. There is a one-to-one map between symmetric `s` (eq:1) and
symmetric `m` (eq:25), with `m` interpretable as a relative minorizing coefficient.

## 5 Summary

Algorithm M can be written in two stages (Algorithm L): reject as "type-x"
(`r1 > min{k gamma(x|y)/p(x), 1}`) or, failing that, as "type-y"
(`r2 > min{p(y)/(k gamma(y|x)), 1}`); otherwise accept `y`. The total duplication
probability equals `1 - alpha_M` (eq:L-twostage).

As `k(x,y)` sweeps from high to low:
- **Case 1, `k >= H`:** `k` is a relative majorizing coefficient (MAR, type-y
  duplications only). As `k` decreases to `H`, `alpha_M = alpha_MA` rises to its maximum
  `alpha_MH`; the paired minorizing coefficient `m_k` rises to `L`.
- **Case 2, `L < k < H`:** `k` is too deficient for MAR alone; both type-x and type-y
  duplications occur, but `alpha_M` stays pinned at `alpha_MH`. Here `M_k = H`, `m_k = L`.
- **Case 3, `k <= L`:** `k` is a relative minorizing coefficient (MIR, type-x
  duplications only). As `k` decreases below `L`, `alpha_M = alpha_MI` falls below
  `alpha_MH`.

Algorithm M unifies MAR (=HA), MIR (=HA) and MH (the optimal HA). It is not more
general than HA, but is easier to understand intuitively.

## Key references (harvested)

Hastings 1970; Metropolis, Rosenbluth, Rosenbluth, Teller & Teller 1953; Barker 1965;
Peskun 1973; Tierney 1994; Chib & Greenberg 1995; Billera & Diaconis 2001; Liu 2001
(Monte Carlo Strategies; contains Stein's algorithm); Minh, Minh & Nguyen 2012
(regenerative MCMC); von Neumann 1951 (Acceptance-Rejection); Minh 2001; Dongarra &
Sullivan 2000.
