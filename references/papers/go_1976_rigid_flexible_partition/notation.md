# Notation - Gō & Scheraga 1976

Reduced-unit / sign conventions:
- $\beta = 1/kT$ throughout (statistical-mechanics convention).
- Zero-point energy of a mode = $\pi\hbar\nu_i = \tfrac12 h\nu_i$ (paper writes $2\pi\hbar\nu_i$ where a full $h\nu_i$ is meant, since $h=2\pi\hbar$).
- Dimensionless mode variable $x_i = 2\pi\hbar\nu_i/kT = h\nu_i/kT$.
- Superscript $+$ (or $\dagger$) denotes matrix/vector transpose.
- "Hard" vs "soft" partition is central: soft = dihedral + external; hard = bond lengths + bond angles.

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $Z$ | configurational partition function | scalar | eq 1 = flexible, eq 2 = rigid |
| $Z_f$ | flexible-model partition function | scalar | eqs 8, 16 |
| $Z_r$ | rigid-model partition function | scalar | eqs 12, 13 |
| $Z_{\rm QM}$ | quantum-correct partition function | scalar | eq 17 |
| $Q$ | soft variables (dihedrals + 6 external) | vector, length $m$ | integration variables in final form |
| $Q'$ | hard variables (bond lengths + bond angles) | vector, length $l$ | integrated out / fixed |
| $q_i, q_i'$ | $i$th soft / hard variable | scalar | |
| $q_{i0}'$ | strain-free (min-energy) value of hard var $i$ | scalar | approximated Q-independent |
| $m$ | number of soft variables | integer | includes 6 external |
| $l$ | number of hard variables | integer | $3n = l + m$ |
| $n$ | number of atoms | integer | |
| $x_{k\alpha}$ | Cartesian coordinate, atom $k$, component $\alpha$ | length | $\alpha\in\{1,2,3\}$ |
| $m_k$ | mass of atom $k$ | mass | |
| $p_{k\alpha}$ | Cartesian momentum $= m_k\dot x_{k\alpha}$ | momentum | |
| $P_r$ | soft generalized momentum, rigid model $=\mathbf{H}^0\dot Q$ | vector, length $m$ | |
| $F(Q)$ | conformational (free) energy / potential of mean force | energy | $F=U(Q)+V(Q)$ (intramol + solvation) |
| $F(Q,Q')$ | full conf. energy incl. hard-variable strain | energy | eq 3 |
| $F_0(Q)$ | min conf. energy over hard vars for given $Q$ | energy | rigid-model energy |
| $f_{ij}''$ | hard-variable harmonic force constants | energy/(coord$^2$) | Q-independent (assumed) |
| $\mathbf{F}''$ | $l\times l$ force-constant matrix, $(i,j)=f_{ij}''$ | matrix | |
| $\mathcal{H}$ | full $3n\times3n$ mass-metric matrix (KE quadratic form) | matrix | eq 5 |
| $\mathbf{H}^0$ | soft-soft block of $\mathcal{H}$ | $m\times m$ | |
| $\mathbf{H}'$ | soft-hard coupling block | $m\times l$ | |
| $\mathbf{H}''$ | hard-hard block | $l\times l$ | |
| $\mathcal{G}$ | $=\mathcal{H}^{-1}$, inverse mass metric | $3n\times3n$ | eq 11 |
| $\mathbf{G}^0,\mathbf{G}',\mathbf{G}''$ | blocks of $\mathcal{G}$ | $m\times m$, $m\times l$, $l\times l$ | |
| $\mathbf{G}$ | soft inverse mass metric $=(\mathbf{H}^0)^{-1}=\mathbf{G}^0-\mathbf{G}'\mathbf{G}''^{-1}\mathbf{G}'^+$ | $m\times m$ | eq 10; the $(\det\mathbf{G})^{-1/2}$ weight |
| $g$ | pure geometric metric tensor (equal-mass case) | $m\times m$ | $\mathbf{H}^0=m_0 g$, eq C-1 |
| $D$ | Cartesian->internal Jacobian | scalar | eq 14; indep. of dihedrals |
| $D(Q_0')$ | Jacobian at strain-free hard values | scalar | |
| $\nu_i$ | vibrational frequency of hard mode $i$ | Hz (cm$^{-1}$ in text) | from $\det(\mathbf{F}''\mathbf{G}'')$, eq 15 |
| $x_i$ | $=2\pi\hbar\nu_i/kT = h\nu_i/kT$ | dimensionless | |
| $\Gamma_r,\Gamma_f$ | QM/classical vib. partition-function ratios | dimensionless | eqs 18, 19 |
| $g_r,g_f$ | per-mode log-sensitivities of $\Gamma_r,\Gamma_f$ | dimensionless | eqs 22, 23; both negative |
| $V$ | system volume | volume | from translational integration |
| $8\pi^2$ | rotational (Euler-angle) phase volume | constant | overall rotation |
| $\hbar$ | reduced Planck constant | action | $h=2\pi\hbar$ |
| $k, T, \beta$ | Boltzmann const, temperature, $1/kT$ | | |
