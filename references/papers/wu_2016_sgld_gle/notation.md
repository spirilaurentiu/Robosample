# Notation - Wu, Brooks, Vanden-Eijnden 2016 (SGLD-GLE)

Standard molecular-dynamics units (CHARMM convention in the examples). All
per-particle vectors are in $\mathbb{R}^3$; index $i$ runs over particles.
Detailed-balance / canonical-sampling convention: temperature enters only
through $kT$ (thermal energy).

| symbol | meaning | units / dtype / shape | convention |
|---|---|---|---|
| $\mathbf{r}_i$ | position of particle $i$ | length (Å), $\mathbb{R}^3$ | |
| $\mathbf{p}_i$ | momentum of particle $i$ | mass·velocity, $\mathbb{R}^3$ | $\mathbf{p}_i = m_i\dot{\mathbf r}_i$ |
| $\dot{\mathbf{r}}_i$ | velocity of particle $i$ | length/time, $\mathbb{R}^3$ | |
| $\dot{\mathbf{p}}_i$ | force (momentum time-derivative) | force, $\mathbb{R}^3$ | |
| $\mathbf{f}_i$ | conservative interaction force $-\nabla_{\mathbf r_i} E_p$ | force, $\mathbb{R}^3$ | from force field (CHARMM in examples) |
| $\mathbf{g}_i$ | SGLD guiding force | force, $\mathbb{R}^3$ | eq. (3) |
| $\mathbf{g}_i^{(GLE)}$ | SGLD-GLE guiding force (auxiliary variable) | force, $\mathbb{R}^3$ | eq. (15); evolves by eq. (16) |
| $\mathbf{R}_i$ | white Gaussian random force | force, $\mathbb{R}^3$ | zero-mean, covariance eq. (2) |
| $\eta_i$ | colored GLE noise | force, $\mathbb{R}^3$ | zero-mean, covariance eq. (8) |
| $\widetilde{\;\cdot\;}$ | low-frequency local time-average operator | same as operand | exp filter eq. (4), time const $t_L$ |
| $\widetilde{\mathbf{p}}_i$ | low-frequency (averaged) momentum | mass·velocity, $\mathbb{R}^3$ | |
| $\widetilde{\mathbf{R}}_i$ | low-frequency (averaged) random force | force, $\mathbb{R}^3$ | |
| $m_i$ | mass of particle $i$ | mass (amu) | |
| $\gamma$ | collision / friction frequency | 1/time (ps$^{-1}$) | same as LD friction |
| $k$ | Boltzmann constant | energy/temperature | $kT$ = thermal energy |
| $T$ | simulation temperature | K | |
| $t_L$ | local averaging time | time (ps) | sets which modes are enhanced |
| $\delta t$ | integration timestep | time (fs/ps) | |
| $\lambda$ | guiding factor (strength of enhancement) | dimensionless | $\lambda \in [0,1]$ in examples; $\lambda=0$ = plain LD |
| $\mu$ | GLE guiding parameter | dimensionless | $\mu \ge 0$, $\lambda = \mu(2-\mu)$, $\mu = 1-\sqrt{1-\lambda}$ (small root used) |
| $\xi$ | energy-conservation factor (SGLD only) | dimensionless | eq. (6); absent in SGLD-GLE |
| $K(t)$ | memory kernel of GLE | 1/time | eq. (9) |
| $\rho$ | phase-space probability density | | extended eq. (18) |
| $E_p$ | potential energy | energy | |
| $\omega$ | angular frequency (spectra) | 1/time | |
| $C(t)$ | velocity autocorrelation function | (velocity)$^2$ | eq. (20) |
| $\delta(t)$ | Dirac delta | 1/time | convention eq. (10): half-sided = 1/2 |
| $\delta_{ij}$ | Kronecker delta | dimensionless | |

## Sign / convention notes

- Guiding force does **no net work**: $\sum_i \mathbf{g}_i\cdot\dot{\mathbf r}_i = 0$ (SGLD, via $\xi$). In SGLD-GLE detailed balance is restored differently (random averaging term), so the $\xi$ term is dropped.
- $\lambda=0 \Rightarrow \mu=0 \Rightarrow \mathbf{g}^{(GLE)}=0$: SGLD-GLE reduces to plain Langevin dynamics.
- Enhancement scales with $\gamma$: guiding force $\propto \gamma$, so gain over LD grows with friction; $\gamma\to 0$ recovers LD.
- The exponential average (eq. 4) is normalized by $1/t_L$ so that $\widetilde{P}$ has the same units/magnitude scale as $P$ for slowly varying $P$.
