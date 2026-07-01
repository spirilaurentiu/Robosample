# Equations - Understanding the Hastings Algorithm (Minh & Minh, 2014)

All algorithms sample a target `p(.)` (un-normalized, `pi = p/P`) via a Markov chain.
Every acceptance rule below is a drop-in for the accept/reject step of a
Metropolis-Hastings-type MCMC move. `gamma(y|x)` is the proposal density of `y`
given current state `x`.

<!-- eq:1 -->
$$ 0 \le \alpha_{HA}(x,y) = s(x,y) \left( 1 + \frac{p(x)}{\gamma(x|y)} \frac{\gamma(y|x)}{p(y)} \right)^{-1} \le 1 $$
- **what:** General Hastings acceptance probability, parameterized by a symmetric function `s(x,y)` that must be chosen so the whole expression stays in `[0,1]` (the "Hastings condition").
- **symbols:** `p(x)` - un-normalized target density at `x`; `gamma(y|x)` - proposal density of `y` given `x`; `s(x,y)=s(y,x)` - free symmetric function; `alpha_HA` - accept probability.

<!-- eq:MT (Metropolis, symmetric proposal) -->
$$ \alpha_{MT}(x,y) = \min\left\{ \frac{p(y)}{p(x)},\, 1 \right\} $$
- **what:** Original Metropolis (1953) acceptance for a SYMMETRIC proposal `gamma(x|y)=gamma(y|x)`.
- **symbols:** as above; requires symmetric proposal so proposal ratio cancels.

<!-- eq:2 -->
$$ s_{MH}(x,y) = \left(1 + \frac{p(x)}{\gamma(x|y)} \frac{\gamma(y|x)}{p(y)}\right) \min\left\{\frac{\gamma(x|y)}{p(x)} \frac{p(y)}{\gamma(y|x)},\, 1\right\} $$
- **what:** The symmetric function `s` that turns the general Hastings rule (eq:1) into Metropolis-Hastings.
- **symbols:** as eq:1; `s_MH` symmetric in `(x,y)`.

<!-- eq:3 -->
$$ \alpha_{MH}(x,y) = \min\left\{ \frac{\gamma(x|y)}{p(x)} \frac{p(y)}{\gamma(y|x)},\, 1 \right\} $$
- **what:** Metropolis-Hastings acceptance probability (the standard rule). Proven optimal (Peskun 1973): highest acceptance of all Hastings algorithms => minimal asymptotic variance.
- **symbols:** ratio `= [p(y) gamma(x|y)] / [p(x) gamma(y|x)]`; capped at 1.

<!-- eq:BK-symmetric (Barker, symmetric proposal) -->
$$ \alpha_{BK}^{(s)}(x,y) = \left(1 + \frac{p(x)}{p(y)}\right)^{-1} $$
- **what:** Barker (1965) acceptance for a SYMMETRIC proposal.
- **symbols:** as above.

<!-- eq:4 -->
$$ \alpha_{BK}(x,y) = \left(1 + \frac{p(x)}{\gamma(x|y)} \frac{\gamma(y|x)}{p(y)}\right)^{-1} $$
- **what:** Generalized Barker acceptance (Hastings rule eq:1 with `s(x,y)=1`). Always < the MH value.
- **symbols:** as eq:1.

<!-- eq:6 -->
$$ \alpha(x,y) = \min\left(\frac{\gamma(x|y)}{p(x)}, 1\right) \min\left(\frac{p(y)}{\gamma(y|x)}, 1\right) \le 1 $$
- **what:** Another special Hastings acceptance (product-of-mins form), from `s` in eq:5. Equals `alpha_M` (eq:8) when `k` is a positive constant.
- **symbols:** as eq:1.

<!-- eq:8 -->
$$ \alpha_M(x,y) = \min\left\{\frac{k(x,y)\gamma(x|y)}{p(x)}, 1\right\} \min\left\{\frac{p(y)}{k(x,y)\gamma(y|x)}, 1\right\} \le 1 $$
- **what:** Algorithm M acceptance probability - unifying form parameterized by a free positive symmetric function `k(x,y)`. Equivalent to the full Hastings family.
- **symbols:** `k(x,y)=k(y,x) > 0` - free symmetric "tuning" function; excluding the constant case `k(x,y)=k`.

<!-- eq:10 -->
$$ \alpha_M(x,y) = \begin{cases} \dfrac{p(y)}{k(x,y)\gamma(y|x)} & \text{if } k(x,y) \ge H(x,y) \\[2mm] \min\left\{\dfrac{\gamma(x|y)}{p(x)} \dfrac{p(y)}{\gamma(y|x)}, 1\right\} & \text{if } L(x,y) < k(x,y) < H(x,y) \\[2mm] \dfrac{k(x,y)\gamma(x|y)}{p(x)} & \text{if } k(x,y) \le L(x,y) \end{cases} $$
- **what:** Piecewise form of Algorithm M in terms of the bracketing thresholds `L,H`. Middle branch = MH (eq:3). This is the practical decision tree.
- **symbols:** `L(x,y)=min{p(x)/gamma(x|y), p(y)/gamma(y|x)}`; `H(x,y)=max{...}` (both symmetric).

<!-- eq:L-H (definitions) -->
$$ L(x,y) = \min\left\{\frac{p(x)}{\gamma(x|y)}, \frac{p(y)}{\gamma(y|x)}\right\}, \qquad H(x,y) = \max\left\{\frac{p(x)}{\gamma(x|y)}, \frac{p(y)}{\gamma(y|x)}\right\} $$
- **what:** Lower/upper bracket functions used throughout; both symmetric in `(x,y)`.
- **symbols:** as above.

<!-- eq:11 -->
$$ M_s(x,y) = \frac{1}{s(x,y)} \left( \frac{p(x)}{\gamma(x|y)} + \frac{p(y)}{\gamma(y|x)} \right) \ge H(x,y) $$
- **what:** Map from a Hastings `s` to an Algorithm-M coefficient `k=M_s` in the majorizing regime (`>= H`). Recovers `alpha_HA` via eq:10 top branch.
- **symbols:** as eq:1; `M_s` symmetric, `>= H(x,y)`.

<!-- eq:13 -->
$$ m_s(x,y) = s(x,y) \left( \frac{\gamma(x|y)}{p(x)} + \frac{\gamma(y|x)}{p(y)} \right)^{-1} \le L(x,y) $$
- **what:** Map from a Hastings `s` to an Algorithm-M coefficient `k=m_s` in the minorizing regime (`<= L`). Recovers `alpha_HA` via eq:10 bottom branch.
- **symbols:** as eq:1; `m_s` symmetric, `<= L(x,y)`.

<!-- eq:15 -->
$$ 0 \le \alpha_{ST}(x,y) = \frac{\delta(x,y)}{p(x)\gamma(y|x)} \le 1 $$
- **what:** Stein algorithm acceptance, parameterized by symmetric `delta(x,y)`. Shown equivalent to the Hastings family.
- **symbols:** `delta(x,y)=delta(y,x)` - symmetric function chosen so ratio in `[0,1]`.

<!-- eq:AR (Acceptance-Rejection accept test) -->
$$ r \le \frac{p(y)}{M\gamma(y)} \le 1 $$
- **what:** von Neumann accept/reject test for i.i.d. sampling. `M` is an absolute majorizing coefficient s.t. `M gamma(z) >= p(z)` for all `z`. Rejected draws are re-drawn (independent samples).
- **symbols:** `r ~ U(0,1)`; `gamma(.)` state-independent proposal; `M` absolute majorizing coeff.

<!-- eq:IMAR (Independence Markovian Accept-Reject) -->
$$ \alpha_{IMA}(x,y) = \frac{p(y)}{M\gamma(y)} \le 1 $$
- **what:** IMAR acceptance. Same test as AR but on rejection the CURRENT `x` is repeated (Markov chain, not i.i.d.). Deficient-M variant: `alpha_D = min{p(y)/(M gamma(y)), 1}`.
- **symbols:** as eq:AR.

<!-- eq:16 -->
$$ M(x,y) \ge \max\left\{\frac{p(x)}{\gamma(x|y)}, \frac{p(y)}{\gamma(y|x)}\right\} = H(x,y) $$
- **what:** Requirement on the RELATIVE majorizing coefficient `M(x,y)` for x-dependent proposals (MAR). Must be symmetric.
- **symbols:** `M(x,y)=M(y,x)` relative majorizing coefficient.

<!-- eq:17 -->
$$ M(x,y) = C(x,y) \max\left\{ \frac{p(x)}{\gamma(x|y)}, \frac{p(y)}{\gamma(y|x)} \right\} $$
- **what:** Parameterization of the relative majorizing coefficient via a symmetric `C(x,y) >= 1`.
- **symbols:** `C(x,y)=C(y,x) >= 1` - free symmetric multiplier; `C=1` gives MH, `C>1` gives Barker-type.

<!-- eq:18 -->
$$ \alpha_{MA}(x,y) = \frac{p(y)}{M(x,y)\gamma(y|x)} $$
- **what:** MAR (Markovian Acceptance-Rejection) acceptance probability in terms of `M(x,y)`.
- **symbols:** as eq:16, eq:17.

<!-- eq:20 -->
$$ \alpha_{MA}(x,y) = \frac{1}{C(x,y)} \min\left\{ \frac{\gamma(x|y)}{p(x)} \frac{p(y)}{\gamma(y|x)}, 1 \right\} \le 1 $$
- **what:** MAR acceptance in terms of `C`. `C=1` => exactly `alpha_MH` (eq:3); this is why MH is the maximum-acceptance member.
- **symbols:** as eq:17.

<!-- eq:21 -->
$$ C(x,y) = \left(\frac{p(y)}{\gamma(y|x)} + \frac{p(x)}{\gamma(x|y)}\right) \left(\max\left\{\frac{p(x)}{\gamma(x|y)}, \frac{p(y)}{\gamma(y|x)}\right\}\right)^{-1} > 1 $$
- **what:** Choice of `C` (in eq:20) that makes MAR reproduce Barker's `alpha_BK` (eq:4).
- **symbols:** as eq:17; this `C_BK > 1`.

<!-- eq:25 -->
$$ m(x,y) \le \min\left\{\frac{p(x)}{\gamma(x|y)}, \frac{p(y)}{\gamma(y|x)}\right\} = L(x,y) $$
- **what:** Requirement on the RELATIVE minorizing coefficient `m(x,y)` for MIR. Must be symmetric.
- **symbols:** `m(x,y)=m(y,x)` relative minorizing coefficient.

<!-- eq:26 -->
$$ m(x,y) = \frac{1}{C(x,y)} \min\left\{ \frac{p(x)}{\gamma(x|y)}, \frac{p(y)}{\gamma(y|x)} \right\} $$
- **what:** Parameterization of the relative minorizing coefficient via symmetric `C(x,y) >= 1`.
- **symbols:** as eq:17.

<!-- eq:27 -->
$$ \alpha_{MI}(x,y) = \frac{m(x,y)\gamma(x|y)}{p(x)} = \frac{1}{C(x,y)} \min\left\{\frac{\gamma(x|y)}{p(x)} \frac{p(y)}{\gamma(y|x)}, 1\right\} \le 1 $$
- **what:** MIR (Markovian Minorizing) acceptance probability. Identical in value to MAR (eq:20) => also equivalent to HA. `C=1` => MH.
- **symbols:** as eq:25, eq:26.

<!-- eq:IMIR (Independence Markovian Minorizing accept test) -->
$$ \alpha_{IMI}(x,y) = \frac{m\gamma(x)}{p(x)} \le 1 $$
- **what:** IMIR acceptance (state-independent proposal). `m` absolute minorizing coeff, `m gamma(z) <= p(z)`. Note: test depends on CURRENT state `x`, not proposal `y`. Deficient-m variant: `alpha_d = min{m gamma(x)/p(x), 1}`.
- **symbols:** `m` absolute minorizing coefficient; `gamma(.)` state-independent proposal.

<!-- eq:L-twostage (Algorithm L, duplication probability) -->
$$ P(\text{x duplicated}) = 1 - \min\left\{\frac{p(y)}{k(x,y)\gamma(y|x)}, 1\right\} \min\left\{\frac{k(x,y)\gamma(x|y)}{p(x)}, 1\right\} $$
- **what:** Total probability that state `x` is repeated in Algorithm M written as a two-stage (type-x then type-y) rejection. `= 1 - alpha_M`. Type-x rejection prob `= 1 - min{k*gamma(x|y)/p(x),1}`; type-y `= min{k*gamma(x|y)/p(x),1}(1 - min{p(y)/(k*gamma(y|x)),1})`.
- **symbols:** `r1, r2 ~ U(0,1)` two independent uniforms; `k(x,y)` free symmetric coefficient.

## Derivations (not implemented)

The bulk of the paper proves equivalences (HA <=> Algorithm M <=> Stein <=> MAR <=> MIR)
and detailed balance `p(x)P(y|x) = p(y)P(x|y)` for each transition kernel. These are
proof-algebra (substituting the `M_s`, `m_s`, `M_k`, `m_k` maps into eq:1/eq:10/eq:15/eq:18/eq:27)
and are not separately implementable beyond the acceptance formulas above. The general
transition kernel form is
$$ P(y|x) = \alpha(x,y)\gamma(y|x) + I(x=y)\int_E (1-\alpha(x,z))\gamma(z|x)\,dz $$
(first term = accepted move to `y`; second = self-loop from rejected proposals).
