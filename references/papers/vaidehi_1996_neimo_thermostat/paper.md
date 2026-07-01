# Constant Temperature Constrained Molecular Dynamics: The Newton–Euler Inverse Mass Operator (NEIMO) Method

Nagarajan Vaidehi, Abhinandan Jain, William A. Goddard III. *J. Phys. Chem.* 1996, 100, 10508–10517. DOI: 10.1021/jp953043o.

## Abstract

The Newton–Euler inverse mass operator (NEIMO) method for internal coordinate molecular dynamics (MD) of macromolecules leads to stable dynamics for time steps about 10 times larger than conventional Cartesian dynamics (e.g., 20–30 fs rather than 1–2 fs for systems containing hydrogens). NEIMO is practical for large systems since the computation time scales linearly with the number of degrees of freedom $\mathcal{N}$ (instead of the $\mathcal{N}^3$ scaling for conventional constrained MD). This paper generalizes the NEIMO formalism to the Nosé (and Hoover) thermostat to derive equations of motion for constrained canonical-ensemble MD, and examines the optimum thermostat mass $Q$ (equivalently the relaxation time $\tau_s$) for NEIMO–Hoover dynamics of polymers. NEIMO–Hoover simulations on amorphous poly(vinyl chloride) and poly(vinylidene fluoride) with time steps of 20–30 fs give stable dynamics.

## 1. Introduction

For studies of conformations and dynamics of polymers and proteins it is useful to constrain structural properties such as bond lengths and bond angles so the focus is on dihedral angles. Constraints reduce the number of degrees of freedom (dof) from $3N$ ($N$ = number of atoms) to $\mathcal{N}$. For polyethylene $C_pH_{2p+2}$, $3N = 9p + 6$ while the number of torsional dof is $\mathcal{N} = p - 1$. Constrained MD also permits larger time steps.

With constrained MD, Newton's equation of motion becomes eq:1, solved for the acceleration via eq:2. The problem is that $\mathcal{M}$ is an $\mathcal{N} \times \mathcal{N}$ matrix, so computing $\mathcal{M}^{-1}$ costs $O(\mathcal{N}^3)$. Since the rest of the calculation is $O(\mathcal{N})$, the inversion dominates for large systems (e.g. inverting a $332 \times 332$ matrix every step for 1001 atoms). NEIMO computes $\ddot{\boldsymbol{\theta}}$ directly without forming or inverting $\mathcal{M}$, using a spatial operator algebra formalism that scales linearly with $\mathcal{N}$.

Standard Newtonian dynamics conserves total energy, describing a *microcanonical* ensemble (E, V, N conserved). To simulate a *canonical* ensemble (T, V, N), the Nosé and Hoover thermostat extensions are applied to the constrained internal-coordinate framework here.

## 2. NEIMO Method in the Microcanonical Ensemble

NEIMO offers an $O(\mathcal{N})$ algorithm for solving eq:1 for constant-energy dynamics. The system is partitioned into rigid bodies (*clusters*) connected by *hinges*, each hinge having 1 to 6 dof. A cluster is a group of atoms moving as a rigid unit (single atom, methylene group, phenyl ring, helix, or an entire domain). Freezing all bonds and angles makes each hinge one-dimensional (a torsion). If only the distance between two clusters is fixed, the hinge is five-dimensional. The union of clusters connected by hinges is a *chain*, which may be open or cross-linked. Each chain has a *base cluster* connected to the reference frame by a full-6-dof hinge, giving each chain freedom to translate and orient.

Within a chain, adjacent clusters are *parent* and *child*. The base has no parent. Branching outward forms a topological tree; each child has exactly one parent. Clusters with no children are *tips*. In a serial chain (linear polymer) there is a single tip and each intermediate cluster has a unique parent and child.

For each cluster, spatial quantities (position, velocity, acceleration, momentum) are computed recursively from the parent cluster. The spatial operator corresponding to the mass matrix is factorized using innovations operator factorization into square, invertible factors, so NEIMO involves no explicit calculation or inversion of the mass matrix. There are three steps:

1. **Velocities.** A base-to-tips recursion: the spatial velocity of each cluster is computed from its hinge motion and its parent's motion.
2. **Forces.** A tips-to-base recursion: the effective force ($\mathcal{T} - \mathbf{C}$) on each cluster is derived from direct forces (Cartesian forces, hinge torques, Coriolis and other velocity-dependent forces) plus indirect forces from its children. Requires velocities (for Coriolis). Quantities related to $\mathcal{M}^{-1}$ are also computed here.
3. **Accelerations.** A final base-to-tips recursion to compute $\ddot{\boldsymbol{\theta}}$ for all clusters.

The microcanonical NEIMO algorithm (combined with POLYGRAF) allows time steps of 20–30 fs even with explicit hydrogens; 20 fs generally gives stable dynamics for amorphous polymers under periodic boundary conditions.

## 3. Canonical Ensemble NEIMO Method

Constant-kinetic-energy methods (velocity/momentum scaling) give incorrect kinetic-energy fluctuations and hence do not produce a canonical momentum distribution. Andersen used stochastic collisions. Nosé made a major advance by (i) extending the system with an additional degree of freedom representing the heat bath, and (ii) choosing a Hamiltonian so the equilibrium dynamics generates the proper canonical distribution in both momentum and configuration space. The Nosé extended system adds one coordinate $s$ (a time-scaling variable) and its conjugate momentum $p_s$. Real time intervals per step are unequal, which complicates Fourier-transform analysis. Hoover proposed an equivalent formulation with explicit real time and real velocities, allowing simpler dynamical-property calculation (FFTs).

### 3.1. Derivation (not implemented) — NEIMO–Nosé Equations of Motion

The equations of motion for the Nosé extended system in the constrained internal-coordinate framework are derived with the Lagrangian formalism. Virtual (Nosé) variables carry a prime; real variables have no prime (opposite the convention of Nosé). The system is placed in a bath at temperature $T_B$, with bath variables $s$ and $p_s$. In the Nosé formulation the real time $t$ relates to the virtual time $t_s$ by eq:3, and the virtual generalized variables are defined by eq:4a and eq:4b.

The extended-system Lagrangian in Nosé variables is eq:5a. Nosé showed that choosing $g = N + 1$ (eq:5b) yields a partition function that, integrated over the bath variables, gives the canonical partition function of the physical system, $Z_{\mathcal{N}} = \langle e^{-H_0/kT_B}\rangle$, so the extended-system trajectory yields a canonical distribution of the physical coordinates and momenta.

From eq:5a, the conjugate momentum is eq:6 (using $M_{jk} = M_{kj}$) and the conjugate-position derivative is eq:7. The Lagrangian equations of motion eq:8, with eq:6 and eq:7 substituted, give eq:9. Because the nonbond contributions to $\partial\Phi/\partial\theta_k$ are tedious in internal coordinates, the potential is split as $\Phi = \Phi_{\rm int} + \Phi_{\rm NB}$ (eq:10), where $\Phi_{\rm int}$ is the internal (e.g. torsional) contribution and $\Phi_{\rm NB}$ contains Coulomb, van der Waals, and external forces. The gradient is eq:11 (sum over the $3N$ Cartesian coordinates). Substituting eq:11 into eq:9 and rearranging gives the NEIMO–Nosé equation of motion eq:12, with Coriolis term eq:13, friction-like term eq:14, and generalized force eq:15. The friction-like $\mathbf{F}_N$ could be absorbed into the Coriolis term $\mathbf{C}_N$; the Coriolis term also carries the nonbond and external forces. NEIMO handles the nonbond/external gradient with respect to Cartesian coordinates (not internal coordinates), which is computationally advantageous.

The Lagrangian equation of motion for $s$ leads to eq:16. Using the physical kinetic energy eq:17 (defining an instantaneous temperature $T$), eq:16 becomes eq:18. Writing $Q$ in terms of the relaxation time $\tau_s$ (eq:19) and substituting gives eq:20; substituting $g = N+1$ gives eq:21. Equations eq:12 and eq:21 are the fundamental NEIMO–Nosé equations.

### 3.2. NEIMO–Hoover Dynamics Equations of Motion

The Nosé equations involve virtual time (unequal real-time steps). Hoover used eq:3 to transform them to real variables via eq:22. Substituting for the Nosé variables in eq:12 gives the NEIMO–Hoover equation of motion eq:23, with Coriolis term eq:24, friction-like term eq:25, and generalized force eq:26. The additional canonical friction term $\zeta \sum_j M_{kj}\dot{\theta}_j$ is included in the Coriolis terms.

Transforming eq:20 to real variables gives eq:27, where the friction coefficient is eq:28. Hoover showed (via the Liouville equation) that the probability density is conserved only if $g = N$ in eq:27, leading to eq:29. Equations eq:23 and eq:29 are the fundamental NEIMO–Hoover equations, solved for each cluster in the tree topology by the recursive algorithm of section 4.

## 4. The NEIMO–Hoover Equations of Motion: Recursive Solution

Since eq:23 has the same form as eq:1, all spatial-operator equations and factorizations hold. The innovations operator factorization provides a closed-form block $\mathbf{LDL}^T$ decomposition of $\mathcal{M}^{-1}$, eq:30. The factor $[\mathbf{I} - \mathbf{H}\psi\mathbf{K}]$ is square, block lower triangular, and nonsingular. $\mathbf{H}$ is the hinge matrix (relative-motion characteristics of the $m$ dof of each hinge; $m=1$ for torsion only); $\psi$ is a lower-diagonal spatial transformation matrix; $\mathbf{K}$ is a spatial operator defined in terms of body inertia and hinge characteristics; $\mathcal{D}$ is block diagonal with $\mathcal{N}$ subblocks of size $m \times m$ ($1\times1$ for torsion only).

For microcanonical NEIMO the accelerations solve eq:31. For NEIMO–Hoover the accelerations become eq:32, equivalently eq:33, where $-\zeta\dot{\boldsymbol{\theta}}$ is the friction-like correction to the microcanonical accelerations. Substituting eq:30 into eq:33 gives the operator expression eq:34, where $\mathbf{P}$ is the articulated-body inertia and $\mathbf{a}$, $\mathbf{b}$ are the Coriolis and gyroscopic forces.

Equation 34 is built from the recursive $O(\mathcal{N})$ sequence eq:34seq. The only difference from microcanonical NEIMO is the additional friction-like term $\zeta\dot{\boldsymbol{\theta}}$.

### Modified leapfrog-Verlet integrator

The leapfrog Verlet algorithm gives velocities and coordinates from accelerations. Because NEIMO needs the velocity at the current step (to compute Coriolis forces, eq:24) to compute the accelerations, the Verlet algorithm is modified:

i. Estimate the current-step velocity as eq:verlet-est.

ii. Use $\dot{\boldsymbol{\theta}}_n$ in the recursive NEIMO–Hoover solution to obtain the accelerations $\ddot{\boldsymbol{\theta}}_n$.

iii. Using $\ddot{\boldsymbol{\theta}}_n$ obtain $\dot{\boldsymbol{\theta}}_{n+1/2}$, then re-estimate $\dot{\boldsymbol{\theta}}_n$ via eq:verlet-recorrect.

iv. Repeat steps ii and iii until convergence. Convergence criterion: velocity difference less than 0.001 MD units (length = Å, time unit = 0.0488 ps). Convergence is generally reached after one or two iterations.

The Hoover variable $\zeta$ is integrated with Verlet, eq:35, with the half-step $D_{n+1/2}$ given by eq:36.

## 5. Optimization of the Nosé Mass Parameter

The mass parameter $Q$ does not affect the canonical distribution but affects the rate of equilibration. Prior Cartesian LJ studies show a broad range of $Q$ (~2 orders of magnitude) gives similar properties, but too small/large $Q$ destabilizes. $Q$ is expressed via the time constant $\tau_s$ (eq:19), the relaxation time of the Nosé variable $s$.

The lower limit on $\tau_s$ is set by the integration time step $\delta$: at least 10 steps per $\tau_s$ period, giving eq:37 hence eq:38 ($\tau_s \ge 1.6\delta$).

Two relaxation time scales are involved: (i) relaxation of interatomic forces (set by fundamental frequencies and mode coupling); (ii) relaxation of the bath variable $s$ (depends on $\tau_s$). Applying $\langle T\rangle = T_B$ to eq:21 with $\langle s\rangle \approx 1$ gives eq:39, a harmonic long-time fluctuation in $\langle s\rangle$ with characteristic time eq:40 ($\tau_{\rm slong} = \sqrt{\mathcal{N}}\,\tau_s$). For good averaging the total simulation time should exceed 20 periods, eq:41.

### 5.1. Tests on pe50

$pe50 = CH_3-(CH_2)_{48}-CH_3$ (eq:42). Each $CH_3$/$CH_2$ group is a cluster, giving 99 hinges and $\mathcal{N} = 99 + 6 = 105$ total dof (including the 6 base-body dof). Because linear and angular momenta are conserved, properties correspond to 99 independent dof. Velocity-autocorrelation analysis places torsion frequencies at 200–600 cm⁻¹ (periods 166–55 fs), suggesting time steps $\delta \approx$ 16–6 fs.

Simulations ran 400 ps at 300 K. For $\delta = 0.01$ ps, eq:38 gives $\tau_s > 0.016$ ps, and simulations blew up at $\tau_s = 0.007$ ps. For $\delta = 0.005$ ps, eq:38 gives $\tau_s > 0.008$ ps, and simulations blew up at $\tau_s = 0.007$ ps. From eq:41, $t_{\rm total} = 400$ ps requires $\tau_s < 0.32$ ps; the $\tau_s = 10.0$ ps run blew up (bath too slow / $Q$ too large). RMS total-energy fluctuations: for $\delta = 1$ fs energy is well conserved at all $\tau_s$; for $\delta = 5$ fs conservation is good except at $\tau_s = 0.01$; for $\delta = 10$ fs deviation is large only for $\tau_s < 0.05$. So $\tau_s = 0.05$ ps gives good energy conservation for all $\delta \le 10$ fs.

**5.1.2 Temperature distributions.** Statistical mechanics gives $\langle KE\rangle = \tfrac{\mathcal{N}}{2}kT_B$ (eq:KE-avg) and fluctuation $\langle(\delta KE)^2\rangle = \tfrac{\mathcal{N}}{2}(kT_B)^2$ (eq:KE-fluc). Hence $\langle T_{\rm calc}\rangle = T_B$ (eq:43) and the mean-square temperature deviation $\langle\delta T_{\rm calc}^2\rangle = \tfrac{2}{\mathcal{N}}T_B^2$ (eq:dT), so eq:44. Best properties for pe50 at $\delta = 5$ fs occur for $0.05\,{\rm ps} \le \tau_s \le 0.10\,{\rm ps}$ (eq:45); at $\delta = 10$ fs for $0.05\,{\rm ps} \le \tau_s \le 0.07\,{\rm ps}$ (eq:46).

### 5.2. PVDF-50

$PVDF50 = CF_3-(CH_2-CF_2)_{49}-CH_3$, MSXX force field, each $CF_3$/$CH_3$/$CH_2$/$CF_2$ a cluster connected by torsion-only hinges: 99 torsional dof. 400 ps, $\delta = 10$ fs. eq:38 demands $\tau_s > 0.016$; $\tau_s = 0.01$ blew up. $\langle T_{\rm calc}\rangle$ is close to $T_B$ except at $\tau_s = 10.0$ ps (eq:41 requires $t_{\rm total} > 12600$ ps). Best at $\tau_s = 0.1$ ps; for $\tau_s > 0.3$ ps the fluctuation deviates increasingly. Best range $0.05 < \tau_s < 0.10$ (eq:47).

### 5.3. Conclusion

Combining eq:45, eq:46, eq:47, the recommended range is $0.05 < \tau_s < 0.07$ ps (eq:48).

## 6. Applications of NEIMO–Hoover Dynamics

Applications use $\tau_s = 0.05$ ps (from eq:45); for $\delta = 25$ or 30 fs, $\tau_s = 0.07$ ps is used (consistent with eq:38 for $\delta < 44$ fs).

**6.1. Isolated polymer chains.** pe20, pe30, pe40, pe50 at 300 K for 400 ps, torsional dof only. Total dof (including base body): 45, 65, 85, 105 respectively. NEIMO–Hoover time steps of 25–30 fs give stable dynamics, while Cartesian–Hoover is limited to $\delta \le 2$–3 fs — about a 10× larger time step for the same total-energy conservation.

**6.2. Amorphous polymers (periodic boundary conditions).**
- **PVDF66** $= CF_3-(CH_2-CF_2)_{65}-CH_3$, one chain per unit cell, 132 carbons, cubic cell $a=b=c=18.0$ Å, 131 torsional dof, $\mathcal{N} = 137$. Each backbone carbon with ligands is a cluster; Karasawa force field; fast Ewald summation for nonbonds; energy minimized to rms force < 0.1 (kcal/mol)/Å. 300 K, $\tau_s \ge 0.05$ ps (eq:38: $\tau_s \ge 0.048$ ps for $\delta = 30$ fs). NEIMO–Hoover $\delta$ can be 10× that of Cartesian–Hoover for equal total-energy fluctuation.
- **(PVC20)₄** $= CH_3-(CCl-CH_2)_{19}-CH_2Cl$, four chains per unit cell, cell $a=22.8$ Å, $b=23.1$ Å, $c=12.1$ Å, $\alpha=94.82°$, $\beta=89.17°$, $\gamma=84.92°$; 488 atoms, 180 torsional dof (including 4 base bodies). Ewald summation; energy minimized to rms force < 0.1 (kcal/mol)/Å. $\tau_s = 0.05$ ps (0.07 ps for $\delta = 30$ fs). $\delta = 20$ fs gives stable dynamics; NEIMO time steps ~10× those of Cartesian dynamics.

## 7. Summary

NEIMO allows an order-of-magnitude speedup for polymers and proteins via larger integration time steps. The method is extended to sample the (T, V, N) canonical ensemble. Lagrangian equations of motion are derived for both NEIMO–Nosé and NEIMO–Hoover, adding a friction-like force to the Coriolis term, and the recursive $O(\mathcal{N})$ algorithm is retained. The heat-bath relaxation constant $\tau_s$ (equivalently mass $Q$) for torsional constrained MD is characterized; time steps of 20–30 fs give stable NEIMO–Hoover dynamics of amorphous polymers.
