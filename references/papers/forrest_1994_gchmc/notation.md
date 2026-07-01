# Notation - Forrest & Suter 1994

Conventions: reduced units with bead mass `mb = 1` unless noted; `β = 1/(kB T)`.
A **tilde** marks fictitious/computational quantities (momenta, moments of
inertia, kinetic energy, time). Angles `φ_k`: indices 1-3 are Euler angles,
`4 ≤ k ≤ Nb` are torsional (dihedral) angles; `φ4` innermost torsion, `φ_Nb`
outermost.

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| `Nc` | number of chains | int | 20 in experiments |
| `Nb` | beads per chain | int | 24 (C24) in experiments |
| `c, c'` | chain index | int | 1..Nc |
| `k` | angle index of a chain | int | 1..Nb (1-3 Euler, 4-Nb torsion) |
| `i, j` | bead/DOF index | int | |
| `r_i^(c)` | position of bead i of chain c | Å, R^3 | Cartesian, lab frame |
| `r_1^(c)` | chain-origin (first bead) position | Å, R^3 | 3 explicit Cartesian coords per chain |
| `φ_k^(c)` | k-th generalized angle of chain c | rad | Euler (k≤3) or torsion (k≥4) |
| `φ2^(c)` | second Euler angle | rad | appears in Jacobian sin factor |
| `{x}` | full Cartesian coordinate set | 3 Nc Nb | unconstrained flexible model |
| `{q}` | generalized coordinate set | Nc(Nb+3) | rigid model: `{r_1^(c), φ_k^(c)}` |
| `E` (`𝓔`) | total potential energy | kJ/mol | eq 1 |
| `V_φ` | torsional (Ryckaert-Bellemans) potential | kJ/mol | eq 2 |
| `V_LJ` | 12-6 Lennard-Jones potential | J/mol | eq 3 |
| `C` | RB torsional prefactor | 9.0 kJ/mol | eq 2 (also used as a normalization const in eq 4) |
| `a0..a5` | RB coefficients | dimensionless | 1, 1.31, -1.414, -0.3297, 2.828, -3.3943 |
| `ε` | LJ well depth | 410 J/mol | all beads |
| `σ` | LJ diameter | 3.94 Å | all beads |
| `lb` | fixed bond length | 1.53 Å | rigid constraint |
| `θ` | fixed bond angle | 112° | rigid constraint |
| `mb` | bead (CH2) mass | 1 (reduced) or 14 g/mol | reduced units unless stated |
| `mc` | mass of chain c | reduced | = Nb·mb |
| `T` | temperature | 480 K | melt; high-T start ≈ 1000 K |
| `kB` | Boltzmann constant | | |
| `β` | inverse temperature | 1/(kB T) | |
| `p_1^(c)` | Cartesian momentum of chain-origin | R^3 | conjugate to r_1^(c) |
| `π̃_k^(c)` | fictitious momentum conjugate to φ_k^(c) | | drawn N(0, Ĩ_k/β) |
| `Ĩ_k^(c)` | effective (fictitious) moment of inertia | reduced | free constant; best = ⟨I_k⟩ |
| `I_k` (`𝓘_k`) | instantaneous moment of inertia of angle φ_k | reduced | eq 14, conformation-dependent |
| `⟨I_k⟩` | equilibrium mean of I_k | reduced | time-independent; recommended Ĩ_k |
| `I_ij` | true inertia tensor (generalized coords) | | eq 12, coupled/time-dependent |
| `ω_k` | actual angular velocity of angle k | rad/time | true dynamics (contrast) |
| `e_k` | unit rotation axis of angle φ_k | R^3, unit | eq 14 |
| `P_k` | point the rotation axis passes through | Å, R^3 | eq 14 |
| `M_k` | set of beads moved by angle φ_k | index set | eq 14 |
| `u_i` | unit vector along i-th bond | R^3, unit | eq 17 |
| `K̃` | fictitious kinetic energy | | eq 8 |
| `H̃` | total (fictitious) Hamiltonian = E + K̃ | | |
| `ΔH̃` | discretization error of H̃ over a trajectory | | eq 11 |
| `ΔH̃*` | ΔH̃ - kB T ln(Πc) | | recovers standard identity, eq 18 |
| `D({q})` | Cartesian→GC Jacobian × Fixman weight | dimensionless | eq 6; Fixman part ignored, only Euler-angle sin retained |
| `PS[a→b]` | proposal probability a→b | | |
| `Pacc` / `PA` | acceptance probability | | eq 7, 11, A3, A5 |
| `t̃` | fictitious (computer) time | | |
| `δt̃_MD` | MD integration time-step | reduced (e.g. 1.5e-3) | leap-frog |
| `N_MD` | MD steps per MC step (trajectory length) | int | optimum ≈ 75 |
| `Δt̃` | trajectory length = N_MD δt̃_MD | | |
| `Δt` | real-time interval per MC step | fs | ≈150 fs; per MD step ≈15 fs |
| `⟨s^2⟩` | mean-square radius of gyration | ≈40 Å^2 | |
| `f_EEV` | end-to-end vector autocorrelation | dimensionless | eq 16 |
| `f_BCF;1`, `f_BCF;2` | bond orientational autocorrelations | dimensionless | eq 17; BCF;2 is P2 Legendre |
| `τ_s` | time for c.o.m. to diffuse one radius of gyration | MD steps | efficiency metric |
| `τ_0.75` | time for f_EEV to decay to 75% | MD steps | efficiency metric |

## Sign / convention notes

- Fixman correction is deliberately **omitted**; only the Euler-angle Jacobian
  contribution (product of `sin φ2` over chains) is kept in acceptance.
- The internal energy in the equations of motion equals the actual internal
  energy here (not fictitious), but the kinetic energy is fictitious.
- Momenta are refreshed (redrawn) **every** MC step regardless of accept/reject -
  required for detailed balance.
- Leap-frog is mandatory: time-reversibility + phase-space-volume conservation
  make the acceptance factor (11)/(A3) valid; other integrators need an extra
  Jacobian.
