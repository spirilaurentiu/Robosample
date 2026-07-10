# Robosample

Robosample is a molecular simulation software implementing blocked Gibbs Sampling coupled with Hamiltonian Monte Carlo. Its goal is to enhance sampling of rare conformational transitions.

Atom positions can be represented in Cartesian (X/Y/Z) or generalized-coordinates. A common generalized-coordinate representation uses internal Bond/Angle/Torsion (BAT) coordinates. Each molecule is represented as an articulated robot - a kinematic tree of rigid bodies connected by joints. Each joint introduces one or more generalized-coordinates corresponding to bond, angle, or torsional motion, while omitted coordinates remain constrained. The molecular system is thus represented as a forest of such kinematic trees sharing a common ground frame.

Each Gibbs block is defined by a robot factorization, which exposes a particular subset of generalized-coordinates. Multiple robot factorizations are applied sequentially over the course of a Gibbs sweep. All robot factorizations represent the same molecular configuration; only the choice of active generalized-coordinates changes between Gibbs blocks. Each Gibbs block samples one subset of generalized-coordinates while conditioning on the remaining coordinates. By holding selected BAT coordinates fixed within a Gibbs block, the dimensionality of the sampled state space is reduced. Within each Gibbs block, Hamiltonian Monte Carlo evolves the active generalized-coordinates together with their conjugate momenta defined by the corresponding robot factorization.

The Cartesian coordinates produced at the end of one Gibbs block are used to initialize the robot representation for the next Gibbs block. **Per-atom** forces computed by OpenMM are accumulated into per-body spatial wrenches (net force and torque about the body origin) for the articulated-body solver. The articulated-body solver uses these body-level wrenches to compute the joint-space dynamics required for HMC integration.

Robosample is a `conda` package that ships:

- A molecular simulation engine written in C++17/CUDA and compiled with PyBind11 as an `.so`.
- A Python library built on top of it.

The build toolchain (`cmake`, `ninja`, `gcc`, `cuda` etc) lives in its own `conda` environment which the user must create and activate. Before running build commands, check that one of the environments in `envs/robo_cuda*.yaml` or `envs/robo_cpu.yaml` is active.

The target systems are up to **1M atoms** clustered in up to **100k rigid bodies** across **10k robots**. Examples:

- Alanine dipeptide in vacuum
- FFAR1 (GPCR) in implicit solvent and with explicit membrane nanodisc
- Spliceosome in explicit solvent
- Entire bacteriophages in explicit solvent

Correctness is measured not against intuition, but against theory (see `references/index.yaml`).

## Agents and references

Each agent's authoritative behavior lives in its own frontmatter under `.claude/agents/`:

- `.researcher.md` - **read-only.** Theory / sampling research -> a precise spec under `docs/specs/`. Runs first, before any physics-touching change. Never writes code.
- `.coder.md` - implements a reviewed spec or a self-contained task. **Auto-mode, surgical**. Correctness over performance - never optimizes speculatively.
- `.reviewer.md` - **independent, hostile, read-only.** Runs after the coder, before merge. Reviews science, conventions and implementation together. Hands confirmed findings back to `coder`, never patches.
- `optimizer.md` - **opt-in, user-invoked only; never inside a feature loop.** Static-analysis-guided source optimization measured with `perf`; hands every change to `reviewer`.
- `references/index.yaml`: summary of all papers Robosample is based on with link to their Markdown versions stored locally.

If in doubt, iteratively ask questions until the implementation plan is complete. After that, enter auto mode and finish coding independently. Do not ask further questions.

## Style

**Write and think** in clear, direct English. Prefer short, declarative sentences where they improve clarity. Introduce technical terms only when they improve precision. Avoid rhetorical flourish. Prefer precise, literal language over metaphor. Explain mechanisms explicitly rather than replacing them with slogans.

Authoritative examples for good and bad writing/thinking are under `styles/`:

- `architecture.md`: explanation of internals / how a system is built.
- `decision.md`: design proposals / rationale documents.
- `documentation.md`: tutorials, how-tos, and conceptual explanation - teachin.
- `issues.md`: bug reports, feature requests, and triage-ready problem writeups.
- `readme.md`: project overview / first-contact document.
- `reference.md`: lookup material - you arrive knowing what you want, you leave with the exact signature/behavior/return codes/errors/edge cases.
- `spec.md`: specification writing.
- `review.md`: code review comments and exchanges - giving and receiving.

Match their style unless the user explicitly requests otherwise.

Assume the reader is an experienced software engineer familiar with molecular simulation. Scale response length to the task. Lead with technical substance. Avoid performative tics and conversational filler:

- Unnecessary validation: *fair point*.
- Narrating the next move: *let me name them plainly*.
- Flagging significance: *this is the real issue*
- Advertising honesty: *to be honest*.
- Corporate jargon and metaphorical engineering slang: *load-bearing*, *blast radius*, *footgun*, *yak shaving*, *belt-and-suspenders*, *fan out*, *clique*, *bespoke*, *circuit-breaker*, *heavy-lifting*, *money shot*, *this lands*, *sidecar* etc.

Avoid excessive normative statements: *binding*, *mandatory*, *fail loud*, *must*, *resolved*, *almost* since everything becomes equally important. Instead, classify requirements:

- **Normative:** *SHALL*, *MUST*, *MUST NOT*.
- **Recommended:** *SHOULD*, *SHOULD NOT*.
- **Optional**: *MAY*.
- **Implementation notes:** *NOTE*.

## Building and testing

A `conda` env must be active since `CONDA_PREFIX` is required by both `cmake` presets and `nox`.

CMake presets `${CONFIG}` are `<platform>-<type>` where:

- `platform` is `reference`, `cpu`, `cuda`, `opencl`.
- `type` is `debug`, `release`, `relwithdebinfo`, `tests`, `pgo-train`, `pgo-use`.

Unless otherwise prompted, you will assume the platform is `cuda`. `tests` does not build the `.so` library, so you cannot run `python/robosample/run.py` or any other driver.

Configure: `cmake --preset ${CONFIG}`

Build & install: `cmake --build --preset ${CONFIG}`.

Every change must either pass the suite of tests and physical invariants stated here or change **deliberately and with justification**. A test encodes WHY behavior matters; "tests pass" is false if any were skipped. Since tests can run for a very long time, you will prompt the user whether to run or not. Regardless of the response, you will always run::

```bash
python3 python/robosample/run.py --name ala-dipeptide --prmtop examples/ala-dipeptide.prmtop --inprcrd examples/ala-dipeptide.rst7 --seed 6000 --equil_steps 0 --prod_steps 10 --write_freq 1 --validate true
```

If prmompted to run tests, execute:

```bash
nox -s tests`
```
