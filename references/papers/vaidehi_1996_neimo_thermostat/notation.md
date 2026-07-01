# Notation — Vaidehi, Jain, Goddard 1996 (NEIMO constant-temperature)

**Unit/sign conventions.** MD units: length = Å, time unit = 0.0488 ps, energy = kcal/mol. $k$ = Boltzmann constant. Primes denote *virtual* (Nosé) variables; unprimed denote *real* variables — this is **opposite the convention of Nosé (1984)**. Repeated Latin indices $i,j$ sum over $\mathcal{N}$ generalized (internal) coordinates; $\alpha$ sums over the $3N$ Cartesian coordinates. Reported temperature fluctuations use the normalized statistic $((\mathcal{N}/2)\langle\delta T_{\rm calc}^2\rangle)^{1/2}$, which equals $T_B$ at equilibrium.

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $N$ | number of atoms | scalar int | — |
| $\mathcal{N}$ | number of generalized (internal) degrees of freedom; includes 6 base-body dof per chain | scalar int | torsion-only: $\mathcal{N} = $ #hinges + 6 |
| $\boldsymbol{\theta}$ | generalized coordinates (torsion angles) | $\mathbb{R}^{\mathcal{N}}$ | — |
| $\dot{\boldsymbol{\theta}}, \ddot{\boldsymbol{\theta}}$ | generalized velocity, acceleration | $\mathbb{R}^{\mathcal{N}}$ | real time |
| $\theta_i', \dot{\theta}_i'$ | virtual (Nosé) coordinate/velocity | $\mathbb{R}$ | $\theta_i=\theta_i'$, $\dot{\theta}_i = s\dot{\theta}_i'$ |
| $\mathcal{M}, M_{ij}$ | mass matrix / moment-of-inertia tensor and its elements | $\mathbb{R}^{\mathcal{N}\times\mathcal{N}}$, symmetric $M_{ij}=M_{ji}$ | never explicitly inverted in NEIMO |
| $\mathcal{T}$ | generalized force (torque) vector | $\mathbb{R}^{\mathcal{N}}$ | subscript N/H = Nosé/Hoover form |
| $\mathbf{C}$ | Coriolis / velocity-dependent forces | $\mathbb{R}^{\mathcal{N}}$ | subscript N/H |
| $\mathbf{F}_N, \mathbf{F}_H$ | thermostat friction-like force | $\mathbb{R}^{\mathcal{N}}$ | folded into Coriolis term |
| $\Phi$ | total physical potential energy | kcal/mol | $\Phi = \Phi_{\rm int}+\Phi_{\rm NB}$ |
| $\Phi_{\rm int}$ | internal (torsional) potential | kcal/mol | gradient wrt $\theta$ |
| $\Phi_{\rm NB}$ | nonbond potential (Coulomb + vdW + external) | kcal/mol | gradient computed wrt Cartesian $X_\alpha$ |
| $X_\alpha$ | Cartesian coordinate $\alpha$ | Å | $\alpha = 1..3N$ |
| $s$ | Nosé bath (time-scaling) coordinate | dimensionless | $\langle s\rangle \approx 1$ at equilibrium |
| $p_s$ | conjugate momentum of $s$ | — | — |
| $\dot{s}, \ddot{s}$ | bath velocity, acceleration | — | — |
| $Q$ | Nosé thermostat mass | energy·time² | tunes equilibration, not the distribution |
| $\tau_s$ | bath relaxation time (encodes $Q$) | ps | $\tau_s^2 = Q/(\mathcal{N}kT_B)$ |
| $\tau_{\rm slong}$ | long-time $\langle s\rangle$ relaxation | ps | $=\sqrt{\mathcal{N}}\tau_s$ |
| $\zeta$ | Hoover friction coefficient | 1/time | $\zeta = (1/s)\,\mathrm{d}s/\mathrm{d}t_s$ |
| $g$ | thermostat dof count | scalar | $g=\mathcal{N}+1$ (Nosé), $g=\mathcal{N}$ (Hoover) |
| $k$ | Boltzmann constant | — | — |
| $T$ | instantaneous temperature | K | from KE via eq:17 |
| $T_B$ | heat-bath (target) temperature | K | 300 K in all runs |
| $T_{\rm calc}$ | time-averaged calculated temperature | K | should equal $T_B$ |
| $\delta$ | integration time step | fs (or ps) | 1–30 fs range tested |
| $t, t_s$ | real time, virtual (Nosé) time | ps | $\mathrm{d}t = \mathrm{d}t_s/s$ |
| $\mathcal{M}^{-1}$ | inverse mass matrix | — | factorized as $[\mathbf{I}-\mathbf{H}\psi\mathbf{K}]^T\mathcal{D}^{-1}[\mathbf{I}-\mathbf{H}\psi\mathbf{K}]$ |
| $\mathbf{H}$ | hinge matrix (relative-motion characteristics, $m$ dof/hinge) | — | $m=1$ for torsion-only |
| $\psi$ | lower-diagonal spatial transformation matrix | — | — |
| $\mathbf{K}$ | spatial operator (body inertia + hinge chars) | — | — |
| $\mathcal{D}$ | block-diagonal factor | $\mathcal{N}$ blocks $m\times m$ | $1\times1$ for torsion-only |
| $\mathbf{P}$ | articulated-body inertia | — | — |
| $\mathbf{a}, \mathbf{b}$ | Coriolis and gyroscopic forces | — | — |
| $\hat{\mathbf{f}}_c$ | Cartesian/external force contribution | — | — |
| cluster | rigid body (group of atoms moving as a unit) | — | 1 atom → whole domain |
| hinge | joint between clusters | 1–6 dof | 1 = torsion, 5 = fixed-distance, 6 = base |
| chain | union of hinge-connected clusters | — | has one base cluster (6-dof hinge) |
