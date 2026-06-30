# Robosample

Robosample is a molecular simulation software implementing blocked Gibbs Sampling coupled with Hamiltonian Monte Carlo [load paper]. Each molecule is represented as a **robot** (kinematic tree of rigid bodies connected by joints, or *mobilizers*) and each Gibbs block consists of the set of joints of a forest a forest of kinematic trees rooted at a shared single ground frame. Each robot is based on the previous robot i.e. Cartesian coordinates from the previous robot are assigned to the next robot. **Per-atom** forces computed by OpenMM are reduced to per-body spatial forces (net force + torque about the body origin) for the articulated-body solver.

Robosample is a `conda` package that ships:

- A molecular simulation library (C++17 + CUDA, shipped as a pybind11 `.so` plus a large Python molecular library).
- A Python library ().

The goal is to accelerate sampling of molecular conformations, maximizing sampling of **rare events** (basin hopping). Correctness is measured against the theory, not against intuition.

The systems it is supposed to simulate are to **1M atoms** clustered in up to **100k rigid bodies**. Examples:

- Alanine dipeptide in vacuum
- Deca alanine in explicit solvent
- FFAR1 (GPCR) in implicit solvent and with explicit membrane nanodisc
- Spliceosome in explicit solvent

## Agents and references

- `.claude/agents/researcher.md`

- `.claude/agents/coder.md`

- `.claude/agents/reviewer.md`

- `.claude/agents/optimizer.md`

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
