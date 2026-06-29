# Robosample

Robosample is constraint-based, multiscale **HMC** molecular sampling software (C++17 + CUDA, shipped as a pybind11 `.so` plus a large Python molecular library). The goal of all work here is to generate configurations whose stationary distribution is `pi(x) ∝ exp(-beta·U(x))`, maximizing **basin hopping** and **intra-basin exploration**. Correctness is measured against the theory, not against intuition.

Each molecule is represented as a **robot** (kinematic tree of rigid bodies connected by joints, or *mobilizers*). For multiple molecules, this becomes a forest rooted at a shared single ground frame. An **O(n) articulated-body recursion** (Featherstone's articulated-body algorithm, in spatial-operator form) propagates positions, velocities, accelerations, and forces up and down the tree without forming or inverting the mass matrix M. **Per-atom** Cartesian forces (`-grad(U)`) computed by OpenMM are reduced to per-body spatial forces (net force + torque about the body origin) for the articulated-body solver.

The recursion is realized in a custom, data-oriented (structure-of-arrays) multibody engine; the global `M`, `M^-1`, `sqrt(M)`, and `det M` are never assembled, only applied as `O(n)` operators.

## Scientific background

### What problem are we solving?

Suppose a molecular system has coordinates `x`. The equilibrium probability of observing a configuration is `pi(x) proportional to e^-(beta * U(x))` where `U` is the potential energy and `beta = 1/(R*T)` is the Boltzmann constant. The entire purpose of the software is to generate a sequence of molecular configurations whose long-term distribution is `pi(x)`. Concretely, we are trying to maximize **basin hopping** and **intra-basin exploration**.

### Molecular system size

Up to 1M atoms clustered in up to 100k rigid bodies. Examples:

- Alanine dipeptide in vacuum
- Deca alanine in explicit solvent
- FFAR1 (GPCR) in implicit solvent, but with explicit membrane nanodisc
- Spliceosome in explicit solvent

## Definitions for default behavior

### **HAMILTONIAN GEOMETRY**

The space of all possible states `(p,q)` (position and momentum) is called the phase space. The Hamiltonian `Hamiltonian(q,p)` is the total energy written as a function of position and momentum and it dictates the complete motion through 2 equations:

- `dq/dt =  d(H)/dp`: Position is steered by how energy varies with momentum.
- `dp/dt = -d(H)/dq`: Momentum by how energy varies with position, but with a flipped sign.

The system does not move towards lower energy, but stays along iso-Hamiltonian levels. At every phase point, the symplectic form maps the energy gradient to a flow that conserves `H`.

### **MARKOV CHAIN**

The sampler produces `s0 -> s1 -> s2 -> s3 -> s...`. The next state depends on the current one: `prob(s(n+1)|s(n),s(n-1),...) = prob(s(n+1)|s(n))`. This is the Markov property. The Markov chain is the central object of the entire theory.

### **STATIONARY DISTRIBUTION**

A distribution `pi` is stationary if one step leaves it unchanged: `pi*TransitionKernel = pi`. If states are already distributed according to `pi`, another sampling step still gives `pi`.

### **DETAILED BALANCE**

Let `a` and `b` be two states and `TransitionKernel` the transition kernel between them. Then, `pi(a) * TransitionKernel(a->b) = pi(b) * TransitionKernel(b->a)`. Interpretation: The probability current flowing from `a` to `b` equals the current flowing backward. Detailed balance implies stationarity. Stationarity does not require detailed balance.
Example: let there be two conformations trans and gauche. At equilibrium:

- 100 transitions/day trans -> gauche
- 100 transitions/day gauche -> trans

The populations remain constant.

### **ERGODICITY**

The chain is ergodic if, starting from any `s0`, trajectory averages converge to expectations under `pi`: `(1/N) * sum f(s_n) -> integral f(x) * pi(x) * dx` where `pi` is defined over all degrees of freedom (not a subspace). Stationarity says the kernel preserves `pi`. Ergodicity says the chain actually reaches `pi`. Both are needed: preservation without reachability is useless. One-line counterexample: the identity kernel `TransitionKernel(s,s)=1` preserves every distribution but never moves. Stationary, not ergodic. In practice, the failure is reducibility: barriers the proposal cannot cross in the run budget. The chain samples one basin correctly and ignores the others. Example: deca-alanine has an alpha-helical basin and an extended beta/PPII basin separated by several kcal/mol. With short-trajectory Cartesian HMC at 300 K:

- 1000 intra-basin transitions/day in alpha
- 1000 intra-basin transitions/day in beta
- 0 inter-basin transitions/day
Detailed balance holds locally. Starting in alpha, the beta basin contributes zero weight instead of its true Boltzmann weight. The chain is non-ergodic with respect to the full pi.

### **HAMILTONIAN**

Let the Hamiltonian energy function be `Hamiltonian(q,p) = U(q) + KE(q,p) + U_Fixman(q) + U_Jacobian(q)`.

For generalized coordinates:

- `U(q)` is the potential energy computed on full-atom model.
- `KE(q,p) = 1/2* p^T *M(q)^-1* p` is the kinetic energy.
- `U_Fixman(q)` is a correction term (see below) and is null in Cartesian space.
- `U_Jacobian(q)` is `ln( sin^2(gamma2) )` summed over all `Free`-rooted molecules. Note that `gamma2` is the pitch of each root's absolute orientation. This term is null in Cartesian space.

For cartesian coordinates, `KE = 1/2 * p^T * M^-1 * p` where `p^T` is the transpose momenta vector.
The Hamiltonian determines the equilibrium distribution: `pi(q,p) proportional to e^(-beta*H(q,p))`

### **ACCEPTANCE**

A proposal is a candidate move: `s -> s'`. Example: draw momenta, integrate dynamics for 100 steps. The proposal is not yet accepted.

The Metropolis rule: `A=min(1,e^(-beta*deltaH))`. If `deltaH == 0`, then accept with probability `1`. If energy increased, accept probabilistically. Acceptance restores detailed balance.

Simple HMC example: draw momenta, integrate Hamilton's equations, compute `deltaH`, accept or reject. This produces a Markov kernel `T(s->s')`

### **REVERSIBILITY**

A proposal is reversible if the reverse path exists:

- Forward: `(q, p)` -> `(q', p')`.
- Reverse: `(q',-p)` -> `(q, -p)`.

Without reversibility, the backward transition probability cannot be computed.

### **VOLUME PRESERVATION**

The proposal must preserve phase-space volume `(dq,dp)`. No compression. No expansion. Otherwise the proposal density changes. The Jacobian would appear in the acceptance probability.

### **SYMPLECTICITY**

Symplectic integrators preserve Hamiltonian geometry. Consequences: phase-space volume preservation, near-energy conservation, shadow Hamiltonian. Energy errors remain bounded: `H(t) - H(0) = O(dt^2)`. This gives high acceptance rates. Symplecticity improves efficiency. It is not strictly required for correctness.

### **GIBBS SAMPLING**

A Gibbs step updates only part of the variables. Examples:

- `q=(q1, q2, q3)`. Update only `q2`. The other coordinates remain fixed.
- torsional dynamics: only torsions move; bond lengths and bond angles are fixed. This constitutes a torsional Gibbs block.

A Gibbs block defines: mobile coordinates, frozen coordinates, mass matrix, Hamiltonian, timestep. Each block is itself an HMC sampler.

### **COMPOSITION**

If `TransitionKernel1` and `TransitionKernel2` are transition kernels which preserve `pi`, then `TransitionKernel2 * TransitionKernel1` also preserves `pi`. This is why Cartesian blocks and torsional blocks can be mixed without introducing biases.

### **GENERALIZED COORDINATES**

Instead of Cartesian coordinates `(x,y,z)`, we use `q=(bond_length, bond_angle, torsion/dihedral)`. The mass matrix then becomes `M(q)`. Based on the accessible degrees of freedom, the space can be:

- Flat (linear Taylor):
  - `Pin` torsions live on the circle `S^1`, which has **zero intrinsic curvature** (a circle is a line made periodic, locally indistinguishable from `R`) and is integrated as an unwrapped real.
  - `Translation` and `Free`  coordinates lie in the flat, zero-curvature space `R^3`.
- Curved path (exponential map): The **orientation** of a `Free` or `Ball` body, stored as a unit quaternion living on **S^3, the unit 3-sphere in R^4** which is a compact manifold of **constant positive curvature** that **double-covers the rotation group SO(3)** (`q` and `-q` are the same physical rotation).

From a physics point of view, bonds and angles vary rapidly and contribute little to the overall RMSD. On the contrary, torsions are softer modes and heavily influence conformational transitions, but these transitions are dependent on all degrees of freedom. For example, torsional dynamics cannot sample alone cis to trans isomerization of alanine dipeptide which is gated by a 1-4 clash that requires angle relaxation. The transition occurs only when this torsions are sampled together with bond lengths.

### **MASS METRIC**

The generalized mass matrix acts as a metric `ds^2 = dq^T * M(q) * dq`. Nearby coordinate changes can have different physical meanings depending on configuration. This is the origin of: Fixman correction to correct for marginalizing over momenta, generalized equipartition and velocity-dependent forces.

### **EQUIPARTITION**

The correct momentum distribution is `p ~ N(0, R*T*M(q))`, not `p ~ N(0, I)`. Otherwise the kinetic temperature is wrong.

### **FIXMAN POTENTIAL**

The coordinate transformation introduces a metric bias. The correction is `U_Fixman(q) = 1/2 * R * T * ln(det(M(q)))`. It removes the sampling distortion caused by generalized coordinates.

Generalized-coordinate (e.g. torsional) sampling is designed so that the marginal distribution of configurations matches the marginal of the **fully-flexible Cartesian** Boltzmann distribution restricted to the sampled subspace; the Fixman term removes the mass-metric artifact that would otherwise distort this marginal, and the mixed Gibbs scan supplies the relaxation of the frozen coordinates.

### **CONSTRAINTS**

Constraints define a manifold `sigma(q)=0`. Examples: ring closure. SHAKE projects positions. RATTLE projects velocities.

## Authoritative references (read on demand, never wholesale)

- `.claude/CODING_RULES.md`: general coding rules and targeted hardware.

- `.claude/HPC_OPTIMIZE.md`: performance doctrine. Used only by the `optimizer` agent.

- `.claude/Z_MATRIX.md`: internal coordinates description and construction.

- `.claude/ROBOTICS.md`: rationale, physics and implementation details of the robotics library.

- `.claude/GCHMC.md`: implementation details and physics of Gibbs sampling coupled with Hamiltoninan Monte Carlo.

- `.claude/ENHANCED_SAMPLING.md`: comparison of enhanced sampling methods.

- `references/index.yaml`: summary of all papers Robosample is based on with link to their Markdown versions stored locally.

## The test gate (non-negotiable)

Every change must either pass the suite of tests and physical invariats stated here or change **deliberately and with justification**. A test encodes WHY behavior matters; "tests pass" is false if any were skipped.

- Fast smoke (<1 min, use while iterating):
  `python3 python/robosample/run.py 2ala tip3p/2ala.prmtop tip3p/2ala.rst7 6000 0 100 1 true`
- Full authoritative gate (before declaring done): `nox -s tests` - will build `cuda-release`.
- A `mamba`/`conda` env must be active (`CONDA_PREFIX` is required by both the presets and `noxfile`).

## Build (CUDA)

Can build `${CONFIG}` as production (`cuda-release` - no tests) or testing (`cuda-tests` - no `.so`/Python built).

- Configure: `cmake --preset ${CONFIG}$`
- Build & install into `python/robosample/`: `cmake --build --preset ${CONFIG}`
- Builds run **locally** on this machine's GPU (`CMAKE_CUDA_ARCHITECTURES=native`). This is why Remote Control (phone) works and Claude Code on the web does not - the cloud has no CUDA toolchain.

## How work is delegated

Research and design are read-only and produce a spec.

Implementation is surgical (Rules 2, 3, 7).

Validation runs the gate.

Review is hostile (Rule 12) and read-only.

See `.claude/agents/`. Optimization (`optimizer`) is a **separate, opt-in** pipeline - never invoke it inside a feature loop.

If in doubt, iteratively ask questions until the implementation plan is complete. After that, enter auto mode and finish coding independently. Do not ask further questions. Everything is permitted (read `.claude/settings.json`).
