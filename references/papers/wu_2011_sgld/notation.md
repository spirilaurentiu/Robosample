# Notation — Wu & Brooks 2011, SGLD

Units: CHARMM convention. Energies in kcal/mol, temperature in K, length in Å, time in ps, mass in
amu. `p_i = m_i ṙ_i`. Tilde (`~`) always denotes a local/evolving average (low-frequency property);
`x - x̃` is the high-frequency counterpart. `⟨·⟩_L` also means the evolving average in this paper.

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| `p_i`, `r_i`, `f_i` | momentum, position, interaction force of atom i | vec3 per atom | `f_i` must include constraint force |
| `p̃_i`, `f̃_i`, `g̃_i` | low-frequency (evolving-averaged) momentum/force/guiding force | vec3 per atom | evolving average, eq:1 |
| `g_i`, `g'_i` | guiding force; uncorrected guiding force `= λ_i γ_i p̃_i` | vec3 per atom | eq:10, eq:A7 |
| `R_i` | random (stochastic) force | vec3 per atom | Gaussian, zero mean, eq:8/A1 |
| `γ_i` | collision (friction) frequency | 1/ps | LD damping; SGLD uses same value |
| `λ_i`, `λ` | guiding factor (input parameter) | dimensionless | `λ=0` recovers LD; `λ>0` accelerates search |
| `ξ` | energy-conservation factor | dimensionless | makes total guiding-force power zero, eq:12/A6 |
| `t_L` | local average time `= L·δt` | ps | sets low/high-freq split; typ. 0.2 ps |
| L | local averaging size (# points) | int | |
| `δt` | time step | ps | 0.001 (1 fs) or 0.002 (2 fs) in examples |
| ϖ | signal frequency (test function) | 1/time | `ϖ_L = 1/t_L` is the split frequency |
| `E_p`, `Ẽ_p` | potential energy; low-frequency potential energy | kcal/mol | |
| `Ē_p` | trajectory-average potential energy | kcal/mol | subtracted in eq:A9 to avoid overflow |
| `E_k`, `Ẽ_k` | kinetic / low-frequency kinetic energy | kcal/mol | eq:5 |
| T | target (bath) temperature | K | |
| `T̃` | low-frequency temperature | K | eq:6, computed live |
| `T̃_0` | reference low-frequency temperature (`T̃` at `λ=0`) | K | eq:22: `T̃_0 = T̃ χ_lf` |
| `T_lf`, `T_hf` | effective temp in low/high-frequency space | K | `T_lf = T/χ_lf` |
| `T_sg` | self-guiding temperature | K | search-ability measure, eq:27 |
| `λ_lf`, `λ_hf` | low/high-frequency energy factors | dimensionless | eq:15/17, eq:A8; `=1` at `λ=0` |
| `χ_lf` | low-frequency collision factor `= T̃_0/T̃` | dimensionless | eq:20, eq:A8 |
| `χ_i` | leap-frog velocity scaling parameter | dimensionless | eq:A11 |
| `N_DF` | number of degrees of freedom | int | |
| `k` | Boltzmann constant | kcal/(mol·K) | |
| Ω | density of states | | in partition function |
| FLF, FHF, GLF, GHF, PPLF, GPLF | running accumulators for factors | scalar | eq:A8 |
| `f_i^CON`, `r_i^CON` | constraint force / constrained position | vec3 per atom | from SHAKE, eq:A13 |
