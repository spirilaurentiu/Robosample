# Notation - Wagner 2013 (GNEIMO / CICMD)

Spatial Operator Algebra (SOA) conventions follow Jain, *Robot and Multibody Dynamics* (2010).
The molecule is a tree of rigid **clusters** (rigid bodies) connected by **hinges** (joints, mobilizers).
Torsional MD = 1-DOF internal hinges; the base cluster carries a full 6-DOF hinge.

| symbol | meaning | units/dtype/shape | convention |
|---|---|---|---|
| $\mathcal{N}$ | number of generalized (hinge) DOF | scalar int | system rigid-body DOF (6) subtracted for thermal target (eq:1) |
| $n$ | number of clusters; base-cluster index | scalar int | clusters numbered tips->base; base = n |
| $T$ | temperature | K | 300-310 K in tests, up to 1050 K in REXMD |
| $k$ | Boltzmann constant | energy/K | printed as bold $\mathbf{k}$ in source |
| $\Re_e$ | system kinetic energy | energy | fraktur-R sub e in source; target set by eq:1 |
| $\dot{\theta}$ | generalized (hinge) velocities | $\mathcal{R}^{\mathcal{N}}$ | stacked over all hinges |
| $\dot{\theta}(k),\,\theta(k)$ | hinge velocity / coordinate of cluster k | per-hinge vector | $\dot{\theta}(n)=\mathcal{V}(n)$ for 6-DOF base |
| $\mathcal{M}$ (script) | mass matrix (generalized) | $\mathcal{N}\times\mathcal{N}$, dense, config-dependent | eq:2 factorization; nonseparable Hamiltonian |
| $\mathcal{M}$ (plain, eq:11) | total system mass $\sum_i m_i$ | scalar (mass) | context distinguishes from mass matrix |
| $I$ | identity operator/matrix | - | $I_6$ = 6x6 identity |
| $H$ | joint/hinge map operator | block; maps hinge->spatial | $H^{*}(n)=I$ for 6-DOF base hinge |
| $\phi$ | SOA rigid-body transformation (shift) operator | block lower-triangular | $\phi(j,k)$ = 6x6 transform between clusters j,k |
| $\tilde{\phi}$ | $\phi - I$ | operator | strictly lower-triangular part |
| $\psi$ | SOA articulated transformation operator | block | distinct from $\phi$; appears in inverse transform eq:4 |
| $\mathcal{K}$ | SOA gain operator | block | from articulated-body forward dynamics |
| $\mathcal{D}$ | block-diagonal articulated hinge inertia | block-diagonal, SPD | $\mathcal{D}^{1/2}$ used to define modal coords |
| $v,\,v(k)$ | modal velocity coordinates | $\mathcal{R}^{\mathcal{N}}$ | independent; drawn from $\mathcal{N}(0,1)$ per Boltzmann |
| $M(k)$ | spatial inertia of cluster k | 6x6 SPD | per-cluster |
| $\mathcal{R}(k)$ | composite rigid-body spatial inertia of cluster k + all children | 6x6 SPD | accumulated tip->base |
| $M_S=\mathcal{R}(n)$ | whole-system spatial inertia at base | 6x6 SPD | referenced to base-cluster frame |
| $\mathcal{V}(k)$ | spatial velocity of cluster k | 6-vector (angular; linear) | body/spatial 6-velocity |
| $\mathcal{V}_{CM}$ | CM spatial velocity | 6-vector | solved from eq:6/eq:9 |
| $\mathfrak{h}_S$ | system spatial momentum at base | 6-vector | fraktur-h in source |
| $E$ | base pick-off operator | $6\times 6n$, $[0_6,\dots,0_6,I_6]$ | selects base cluster |
| $\delta_V$ | base-cluster spatial-velocity correction | 6-vector | $=-\mathcal{V}_{CM}$ to null momentum |
| $\zeta$ | Nose-Hoover thermostat friction | 1/time | bath coupling |
| $s_t$ | Nose bath variable | dimensionless | $\propto$ bath potential energy; $\ln s_t$ used |
| $KE_{CM}$ | CM kinetic energy | energy | diagnostic for flying-ice-cube |
| $c$ | fit constant | energy | $\propto$ initial CM kinetic energy |
| $B$ | crystallographic B-factor | Å$^2$ | $B=(8\pi^2/3)\mathrm{RMSF}^2$ |
| $\mathrm{RMSF}$ | root-mean-square fluctuation (per residue) | Å | vs time-average position |
| CRMSD | Cartesian backbone RMSD vs crystal | Å | folding/refinement metric |

Conventions / gotchas:
- `*` superscript = transpose/adjoint of an operator.
- The `-6` in eq:1 removes overall translation+rotation; the base cluster owns those 6 DOF.
- Modal coords v are the ONLY coordinates in which equipartition holds for CICMD; do not draw $\dot\theta$ directly from Boltzmann.
- Hamiltonian is **nonseparable** (KE depends on configuration) -> symplectic Verlet is not directly usable; use Lobatto IIIa-b / RK4 / adaptive CVODE.
- All O(N) sweeps: composite inertia $\mathcal{R}$ tip->base; modal-to-hinge transform (eq:4) base->tips scatter.
