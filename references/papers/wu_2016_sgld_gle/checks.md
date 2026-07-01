# Checks / fixtures - Wu, Brooks, Vanden-Eijnden 2016 (SGLD-GLE)

Concrete numbers from the paper for regression testing an implementation.

## Parameter relation (unit test, exact)

- Given guiding factor $\lambda$, expect GLE parameter roots $\mu = 1 \pm \sqrt{1-\lambda}$ from $\lambda = \mu(2-\mu)$ (eq. 13).
  - $\lambda = 0 \Rightarrow \mu \in \{0, 2\}$ (physical/used root $\mu=0$; guiding force vanishes).
  - $\lambda = 1 \Rightarrow \mu = 1$ (double root).
  - $\lambda = 0.75 \Rightarrow \mu \in \{0.5, 1.5\}$.
- Given $\mu$, expect $\lambda = \mu(2-\mu)$. E.g. $\mu=0.5 \Rightarrow \lambda=0.75$.

## Exponential-average incremental update (unit test)

- Given $\widetilde{P}(t-\delta t)$ and $P(t)$, expect
  $\widetilde{P}(t) = (1-\delta t/t_L)\,\widetilde{P}(t-\delta t) + (\delta t/t_L)\,P(t)$ (eq. 4).
  - With $t_L = 0.2$ ps, $\delta t = 0.002$ ps (2 fs): mixing weight $\delta t/t_L = 0.01$.
  - With $t_L = 0.2$ ps, $\delta t = 0.001$ ps (1 fs, argon): weight $= 0.005$.

## Canonical-sampling correctness (statistical / integration test)

- Given SGLD-GLE run at fixed $T$, expect the marginal position-momentum
  distribution (eq. 19) to be the exact canonical (NVT / NPT) Maxwell-Boltzmann,
  independent of $\lambda$. Operationally: potential-energy distribution and
  $\psi$-angle distribution must overlap the plain-LD result (unlike SGLD or
  high-T LD, which shift).

## Alanine dipeptide fixture

- Force field: CHARMM all-atom; distance-dependent dielectric $= 4r$; no non-bonded cutoff.
- Timestep: 2 fs; SHAKE on all bonds; length 20 ns; save every 2 ps.
- Temperature: 300 K (high-T LD run at 350 K); collision frequency $\gamma = 10$/ps.
- SGLD / SGLD-GLE local averaging time: $t_L = 0.2$ ps.
- $\phi$-$\psi$ metastable regions: I around $(-90°, -70°)$, II around $(-90°, 160°)$; small barrier between them.
- Transition counts (I <-> II) for SGLD-GLE: $\lambda=0$ gives **353** transitions, $\lambda=1$ gives **515** transitions, with average potential energy nearly constant across $\lambda$.
- Expectation: SGLD-GLE at $\lambda=1$ has more transitions than LD and than high-T LD, while its $\psi$ and potential-energy distributions match LD (canonical).

## Liquid argon fixture (NPT)

- Potential: Lennard-Jones 6-12, $\varepsilon = 119.8$ K (i.e. $\varepsilon/k = 119.8$ K), $\sigma = 3.405$ Å.
- System: 500 argon atoms in cubic periodic box $28.53 \times 28.53 \times 28.53$ Å$^3$.
- Timestep: 1 fs; collision frequency $\gamma = 10$ ps$^{-1}$; $T = 100$ K (high-T LD at 150 K); target pressure $P = 1$ atm.
- Length: 10 ns; store coordinates + velocities every 0.05 ps.
- Guiding factor $\lambda = 1$ for SGLD and SGLD-GLE runs.
- Expectation: SGLD-GLE average energy and volume vs diffusion constant nearly identical to LD (canonical), while diffusion constant is significantly larger; high-T LD and SGLD shift energy/volume upward.

## Random-force covariance (property test)

- Given SGLD white noise, expect $\langle \mathbf{R}_j(0)\mathbf{R}_i(t)\rangle = 2 m_i kT \gamma\,\delta(t)\delta_{ij}$ (eq. 2).
- Given colored GLE noise, expect $\langle \eta_i(t)\eta_j(t')\rangle = \delta_{ij} m_i kT \gamma\, K(t-t')$ with $K(t)=2\delta(t)-(\lambda/t_L)e^{-t/t_L}$ (eqs. 8, 9).
