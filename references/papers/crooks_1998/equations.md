# Crooks (1998) - Implementable equations

<!-- eq:1 -->
$$ \overline{e^{-\beta W}} = e^{-\beta \Delta F} $$
- **what:** Jarzynski/Crooks work equality. The path-averaged exponential of the negative reduced switching work equals the exponential of the negative reduced free energy difference. Estimator for $\Delta F$ from nonequilibrium switching runs.
- **symbols:** $W$ - total switching work along a path (energy); $\Delta F = F(\beta,\lambda_\tau) - F(\beta,\lambda_0)$ - Helmholtz free energy difference (energy); $\beta = 1/k_B T$ (1/energy); overbar - average over forward paths starting from canonical equilibrium at $\lambda_0$.

<!-- eq:2 -->
$$ \langle e^{-\beta W} \rangle_0 = e^{-\beta \Delta F} $$
- **what:** Infinitely-fast-switching (Zwanzig thermodynamic-perturbation) limit of Eq. (1). $W$ is the instantaneous work to change $\lambda$ while holding the configuration fixed; average is over the equilibrium ensemble at the initial $\lambda$.
- **symbols:** $\langle\cdot\rangle_0$ - equilibrium average with $\lambda$ fixed at initial value; $W$ - instantaneous work $E(i,\lambda_1)-E(i,\lambda_0)$ at fixed configuration $i$.

<!-- eq:3 -->
$$ P(A \mid \lambda) = \frac{e^{-\beta E(A, \lambda)}}{\sum_{i} e^{-\beta E(i, \lambda)}} = \exp[\,\beta F(\beta, \lambda) - \beta E(A, \lambda)\,] $$
- **what:** Canonical equilibrium probability of internal state $A$ at fixed control parameter $\lambda$.
- **symbols:** $A$ - internal state; $E(A,\lambda)$ - state energy (energy); $F(\beta,\lambda) = -\beta^{-1}\ln\sum_i e^{-\beta E(i,\lambda)}$ - Helmholtz free energy; sum over all states.

<!-- eq:5 -->
$$ W = \sum_{t=0}^{\tau-1} [\,E(i_t, \lambda_{t+1}) - E(i_t, \lambda_t)\,] $$
- **what:** Total work performed on the system over the protocol: sum of energy changes from each control-parameter update at fixed internal state. This is the quantity accumulated per switching step and fed into Eqs. (1)/(2).
- **symbols:** $i_t$ - internal state at step $t$; $\lambda_t$ - control parameter at step $t$; $\tau$ - number of time steps (protocol length).

<!-- eq:6 -->
$$ Q = \sum_{t=1}^{\tau} [\,E(i_t, \lambda_t) - E(i_{t-1}, \lambda_t)\,] $$
- **what:** Total heat exchanged with the reservoir: sum of energy changes from each internal-state evolution at fixed control parameter.
- **symbols:** $Q$ - heat absorbed by the system from the bath (energy); other symbols as in Eq. (5).

<!-- eq:7 -->
$$ \Delta E = Q + W = E(i_\tau, \lambda_\tau) - E(i_0, \lambda_0) $$
- **what:** First law: total energy change equals heat plus work, and equals the endpoint energy difference. Sign convention: $Q$ and $W$ both add energy to the system.
- **symbols:** $\Delta E$ - total energy change (energy); $i_0,\lambda_0$ - initial state/parameter; $i_\tau,\lambda_\tau$ - final state/parameter.

<!-- eq:8 -->
$$ \frac{P(A \xrightarrow{\lambda} B)}{P(A \xleftarrow{\lambda} B)} = \frac{P(B \mid \lambda)}{P(A \mid \lambda)} = \frac{e^{-\beta E(B, \lambda)}}{e^{-\beta E(A, \lambda)}} $$
- **what:** Single-step detailed balance (microscopic reversibility) at fixed control parameter. $P(A\xleftarrow{\lambda}B)$ denotes the reverse transition $B\to A$ at the same $\lambda$.
- **symbols:** $P(A \xrightarrow{\lambda} B)$ - transition probability $A\to B$ at fixed $\lambda$; $E(A,\lambda),E(B,\lambda)$ - state energies.

<!-- eq:9 -->
$$ \frac{P(i_0 \xrightarrow{\lambda_1} i_1 \cdots \xrightarrow{\lambda_\tau} i_\tau)}{P(i_0 \xleftarrow{\lambda_1} i_1 \cdots \xleftarrow{\lambda_\tau} i_\tau)} = e^{-\beta Q} $$
- **what:** Multi-step (generalized) detailed balance: the ratio of forward-path to reverse-path transition probability equals $e^{-\beta Q}$, with $Q$ the heat exchanged along the forward path. Holds regardless of the work performed.
- **symbols:** numerator - forward path transition probability (product of single steps, Eq. path-factor); denominator - reverse path transition probability; $Q$ - total heat (Eq. 6).

<!-- eq:10 -->
$$ \frac{P(i_0 \mid \lambda_0)\, P(i_0 \xrightarrow{\lambda_1} i_1 \cdots \xrightarrow{\lambda_\tau} i_\tau)}{P(i_\tau \mid \lambda_\tau)\, P(i_0 \xleftarrow{\lambda_1} i_1 \cdots \xleftarrow{\lambda_\tau} i_\tau)} = e^{\beta \Delta E - \beta \Delta F}\, e^{-\beta Q} = e^{\beta W_d} $$
- **what:** Crooks path-ensemble ratio. With equilibrium endpoints, the ratio of forward to reverse joint path probability equals $e^{\beta W_d}$, the exponential of the reduced dissipative work. This is the core identity used to convert forward-path averages to reverse-path averages.
- **symbols:** $P(i_0\mid\lambda_0)$ - initial equilibrium weight; $P(i_\tau\mid\lambda_\tau)$ - final equilibrium weight; $\Delta E = Q+W$; $W_d = W - \Delta F$ - dissipative work (energy); $\Delta F$ - reversible work.

<!-- eq:13 -->
$$ \frac{P(A \xrightarrow{\lambda} B)}{P(A \xleftarrow{\lambda} B)} = \frac{e^{-\beta E(B, \lambda) - \beta p V(B, \lambda)}}{e^{-\beta E(A, \lambda) - \beta p V(A, \lambda)}} = e^{-\beta Q - \beta p \Delta V} $$
- **what:** Isothermal-isobaric detailed balance. Same structure as Eq. (8) with the enthalpic $pV$ term added; yields the NPT generalization $\overline{e^{-\beta W}} = e^{-\beta \Delta G}$.
- **symbols:** $p$ - pressure; $V(A,\lambda)$ - volume of state $A$; $\Delta V$ - volume change ($V(B)-V(A)$); $\Delta G$ - Gibbs free energy difference. <!-- CHECK: denominator pV corrected to V(A,lambda) (raw OCR read V(B,lambda)) so the ratio gives e^{-beta p Delta V}. -->

## Derivations note

Equations (11) and (12) (Section 3) are the proof of Eq. (1): substitute Eq. (10)
to rewrite the forward-path average $\overline{e^{-\beta W}}$ as a reverse-path
average, then use path independence of $\Delta F$ and normalization of the reverse
path measure. They are proof-algebra, not standalone implementable formulas; kept
in `paper.md` under the Derivation marker.
