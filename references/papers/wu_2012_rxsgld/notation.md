# Notation - Wu 2012 RXSGLD

Conventions: `β = 1/(kT)`. Bold = per-particle 3-vector. A tilde `~` over any
quantity denotes its low-frequency (locally time-averaged, eq:6) component;
`P - P_tilde` is the high-frequency component. Superscript `(m)` = stage m; `[i]`
or `(i)` = replica/conformation index. Reduced energy units of `kT_0` are used in
the double-well test (`T_0 = 50 K`).

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| `i` | particle index | int | sum over all particles |
| `m, n` | stage indices | int | m alternately odd/even for exchange sweeps |
| `p_i` | momentum of particle i | vector R^3, mass*length/time | |
| `p_tilde_i` | low-frequency momentum | vector R^3 | eq:6 applied to p |
| `p_dot_i` | time derivative of momentum | vector R^3 | eq:1 |
| `f_i` | interaction (force-field) force | vector R^3 | |
| `g_i` | guiding force | vector R^3 | eq:3 |
| `R_i` | random (Langevin) force | vector R^3 | eq:2 |
| `r_dot_i` | velocity of particle i | vector R^3 | |
| `m_i` | mass of particle i | mass | |
| `γ_i` | collision (friction) frequency | 1/time | 100/ps (double well); 1/ps (peptide) |
| `λ_i` | guiding factor (per particle) | dimensionless | λ_i=0 -> plain Langevin |
| `ξ` | energy conservation factor | dimensionless | recomputed each step, eq:5 |
| `k` | Boltzmann constant | energy/temperature | |
| `T` | simulation temperature | K | |
| `β` | inverse temperature 1/(kT) | 1/energy | |
| `t_L` | local averaging time | time | 0.2 ps in all runs |
| `δt` | time step | time | 1 fs (double well) |
| `t_est` | estimation time for evolving averages | time | typically 10 t_L |
| `P, P_tilde` | any property and its low-freq part | | eq:6 |
| `E_p` | potential energy | energy | |
| `E_tilde_p` | low-frequency potential energy | energy | |
| `E_p - E_tilde_p` | high-frequency potential energy | energy | |
| `λ_lf, λ_hf` | low/high frequency energy factors | dimensionless | =1 for LD; eq:8, eq:9 |
| `χ_lf, χ_hf` | low/high frequency collision factors | dimensionless | =1 for LD; eq:10, eq:11 |
| `T_tilde` | low-frequency temperature | K | eq:12 |
| `T_tilde_0` | reference low-freq temperature (zero guiding) | K | |
| `T_SG` | self-guiding temperature | K | eq:16; search ability metric |
| `T_SG^0` | target self-guiding temperature | K | tuning target for {λ_i} |
| `Θ_SGLD` | SGLD configurational partition function | dimensionless | eq:7 |
| `Θ_LD` | canonical (LD) partition function | dimensionless | eq:13 |
| `w_SGLD` | SGLD reweighting factor | dimensionless | eq:14 |
| `N_DF` | number of degrees of freedom | int | |
| `k` (stage count) | number of stages above base | int | total stages = k+1 = 8 in this work |
| `μ_tilde_m` | low-frequency exchange coefficient | 1/energy | eq:19a; =0 for TRXLD |
| `μ_m` | high-frequency exchange coefficient | 1/energy | eq:19b; =β_m for TRXLD |
| `ρ_SGLD` | SGLD distribution probability | | eq:18 |
| `π_RX` | replica exchange probability | dimensionless | eq:20; accept min{1,π} |
| `s_mn` | momentum temperature-scaling factor | dimensionless | sqrt(T_n/T_m), eq:24 |
| `s_tilde_mn` | low-freq momentum scaling factor | dimensionless | sqrt(T_tilde_n/T_tilde_m), eq:26 |
| `X_m^[i]` | conformation i at stage m | configuration | |
| `ε_p(x,y,z)` | double-well test potential | energy | eq:27 |
| `a` | x,z stiffness parameter | energy | 20000 kT_0 |
| `b` | double-well depth parameter | energy | 160 kT_0 |
| `w` | well separation | length | 2 Å |
| `s` | skew parameter | energy | 0, kT_0, or 2kT_0 |
| `T_0` | base temperature (double well) | K | 50 K |
| `r_xz` | radial coordinate sqrt(x^2+z^2) | length | |
| `x_1` | fraction of distribution in well near y=0 | dimensionless | Table I observable |
| `s_i` | clustering subset variable | | SIC method |
| `R_i(j)` | region j of subset i | | |
| `I_i` | region index of subset i | int in 1..k_i | eq:31 |
| `k_i` | number of regions of subset i | int | |
| `N_c` | total possible clusters | int | prod_i k_i |
| `CSR` | conformational searching relevancy | fraction | base-stage / all-stage clusters |
