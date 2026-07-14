# Robosample

Robosample is a molecular simulation software implementing blocked Gibbs Sampling coupled with Hamiltonian Monte Carlo. Its goal is to enhance sampling of rare conformational transitions.

Atom positions can be represented in Cartesian (X/Y/Z) or generalized-coordinates. A common generalized-coordinate representation uses internal Bond/Angle/Torsion (BAT) coordinates. Each molecule is represented as an articulated robot - a kinematic tree of rigid bodies connected by joints. Each joint introduces one or more generalized-coordinates corresponding to bond, angle, or torsional motion, while omitted coordinates remain constrained. The molecular system is thus represented as a forest of such kinematic trees sharing a common ground frame.

Each Gibbs block is defined by a robot factorization, which exposes a particular subset of generalized-coordinates. Multiple robot factorizations are applied sequentially over the course of a Gibbs sweep. All robot factorizations represent the same molecular configuration; only the choice of active generalized-coordinates changes between Gibbs blocks. Each Gibbs block samples one subset of generalized-coordinates while conditioning on the remaining coordinates. By holding selected BAT coordinates fixed within a Gibbs block, the dimensionality of the sampled state space is reduced. Within each Gibbs block, Hamiltonian Monte Carlo evolves the active generalized-coordinates together with their conjugate momenta defined by the corresponding robot factorization.

The Cartesian coordinates produced at the end of one Gibbs block are used to initialize the robot representation for the next Gibbs block. Per-atom forces computed by OpenMM are accumulated into body-level spatial wrenches (net force and net torque about each body frame) for the articulated-body solver. The articulated-body solver uses these body-level wrenches to compute the joint-space dynamics required for HMC integration.

Robosample is a Linux-only `conda` package that ships:

- A molecular simulation engine written in C++17/CUDA and compiled with PyBind11 as an `.so`.
- A Python library built on top of it.

The build toolchain (`cmake`, `ninja`, `gcc`, `cuda` etc) lives in its own `conda` environment which the user SHALL create and activate. Before running build commands, check that one of the environments in `envs/robo_cuda*.yaml` or `envs/robo_cpu.yaml` is active.

The target systems are up to **1M atoms** clustered in up to **100k rigid bodies** across **10k robots**. Examples:

- Alanine dipeptide in vacuum
- FFAR1 (GPCR) in implicit solvent and with explicit membrane nanodisc
- Spliceosome in explicit solvent
- Entire bacteriophages in explicit solvent

## Agents

Each agent's authoritative behavior lives in its own file under `.claude/agents/`. The one-line docs below are pointers, not the contract.

Implementation and review:

- `coder.md` - implements a reviewed spec or a self-contained task in C++17 / CUDA / Python. Auto-mode and surgical; correctness over performance, never optimizes speculatively.
- `reviewer.md` - independent, adversarial review of a change or module before merge. Reviews science, conventions, and implementation together; delivers confirmed findings as reproducers under `tests/` and hands them to `coder`, never patches source.
- `optimizer.md` - user-invoked only, never inside a feature loop. Profiles first, then optimizes for speed while preserving behavior; hands every change to `reviewer`.
- `documenter.md` - writes contract-focused Doxygen for C/C++/CUDA symbols and gtest files, derived from behavior across all call sites. Documents only; never changes code.
- `architect.md` - recovers the design already present in the bloated C++ engine and emits it as `ARCHITECTURE.md`, `MODULES.md`, `TESTS.md`, and split/doc tickets. Preserves behavior; proposes no behavioral change.

Research pipeline:

- `research-orchestrator.md` - drives the deep-research pipeline end to end and is the only agent that addresses the user. Routes tiered evidence between stages; never derives, verifies, or synthesizes.
- `scout.md` - Stage-0 orientation. Cheap, wide, parallel discovery: terminology map, candidate canonical sources, convention-contradiction flags. Does not derive or verify.
- `bibliometrics.md` - Stage-1 citation-graph construction from the approved canonical set. Clusters papers into research programmes and attaches scite status; hop-1 by default.
- `corpus.md` - cross-cutting retrieval from `references/` under strict passage hygiene. Returns only the requested passage or symbol verbatim, with a provenance id. The only agent that returns authority-tier text.
- `researcher.md` - deep statistical-mechanics / enhanced-sampling research mapped onto the literature and codebase, producing a precise spec under `docs/specs/`. Read-only; never writes code.
- `spec-reviewer.md` - adversarial review of the researcher's spec before it is spent on verifier runs or implementation. Audits tier discipline, provenance, and discriminating power. Reviews specs, not implementations.
- `verifier.md` - Stage-2 deterministic reviewer. Translates one classified claim into executable evidence and runs it, returning pass/fail plus a re-runnable artifact. A falsification gate, not a certifier.

Knowledge base:

- `paper-ingestor.md` - user-invoked, single-paper. Ingests one PDF-converted `.md` into `references/`: dedups, classifies, cleans to `paper.md`, extracts equations / notation / numeric-check fixtures, validates the LaTeX, and upserts `references/index.yaml`. Writes only under `references/`; returns a one-line status.

The `ingest-papers` slash command (`.claude/commands/ingest-papers.md`) fans one `paper-ingestor` instance out per file across a directory or glob. It is a command, not an agent - its frontmatter has `argument-hint`/`allowed-tools` and no `name:`.

Agents SHALL write in clear, direct English. Prefer short (under 25 words), declarative sentences where they improve clarity. Prefer active voice. Avoid marketing language, rhetorical flourish, and metaphors. Prefer concrete nouns over abstractions. State mechanisms explicitly rather than replacing them with slogans.

Two vocabularies with opposite disciplines:

- Expression vocabulary - any text a reader consumes (writeups, specs, docs, reviews). Precise and minimal: one term per concept, and a technical term only where it improves precision.
- Retrieval vocabulary - search queries and the reformulation that drives them in research mode. Recall-maximizing: many synonyms, cross-discipline names, deliberately broad. The precision rules above SHALL NOT be applied to queries; narrowing a query to a preferred term suppresses recall. See `styles/research.md`.

Style governs reader-visible output, not private reasoning:

- Surface rules (sentence length, active voice, no emphatics, no marketing) apply to reader-visible artifacts only. Chain-of-thought stays exploratory and MAY hedge; do not impose surface rules on it.
- Epistemic rules apply to reasoning and to every intermediate handoff between agents: keep measured, derived, and conjectured claims separate; label a hypothesis and fence it off; quantify uncertainty ("probably < 1/6, certainly < 1/4") instead of shrugging; state evidence before the claim it supports.

Style exemplars for reader-visible artifacts live under `styles/`; each pairs good and bad examples for one output category. Output visible to the user SHALL be written in the matching style. The load is lazy: an agent reads a style only when it is about to write the matching artifact, reads only the part it needs, and an invocation that produces no reader-visible artifact loads nothing. Map artifact to style here:

| Artifact | Style |
| --- | --- |
| Spec (`docs/specs/`) | `styles/spec.md` |
| Decision record (`docs/decisions/`) | `styles/decision.md` |
| Recovered design (`ARCHITECTURE.md`, module maps) | `styles/architecture.md` |
| API / symbol contract, Doxygen, reference docs | `styles/reference.md` |
| Tutorial, how-to, conceptual explainer | `styles/documentation.md` |
| Research report, reformulation, discovery writeup | `styles/research.md` |
| Code-review or spec-review findings | `styles/review.md` |
| Bug report, triage writeup | `styles/issues.md` |
| Project overview, README | `styles/readme.md` |

Assume the reader is an experienced software engineer familiar with molecular simulation. Scale response length to the task. Lead with technical substance. Avoid performative tics and conversational filler:

- Unnecessary validation: *fair point*.
- Narrating the next move: *let me name them plainly*.
- Flagging significance: *this is the real issue*
- Advertising honesty: *to be honest*.
- Corporate jargon and metaphorical engineering slang: *load-bearing*, *blast radius*, *footgun*, *yak shaving*, *belt-and-suspenders*, *fan out*, *clique*, *bespoke*, *circuit-breaker*, *heavy-lifting*, *money shot*, *this lands*, *sidecar* etc.

Avoid emphatics used for emphasis rather than obligation: *binding*, *mandatory*, *fail loud*, *resolved*, *almost*, and lowercase *must*. They flatten priority when everything reads as equally important. NOTE: this bans vague emphasis, not quantified uncertainty - a hedge that carries a number ("probably < 1/6") is required in reasoning and handoffs, not discouraged. Reserve requirement force for the RFC 2119 keywords (see `styles/spec.md`) and classify each requirement:

- **Normative:** *SHALL*, *MUST*, *MUST NOT*.
- **Recommended:** *SHOULD*, *SHOULD NOT*.
- **Optional**: *MAY*.
- **Implementation notes:** *NOTE*.

## Build prerequisites

A `conda` env SHALL be active (`CONDA_PREFIX` set).

## Build configurations

CMake presets `${CONFIG}` are `<platform>-<type>` where:

- `platform` is `reference`, `cpu`, `cuda`, `opencl`.
- `type` is `debug`, `release`, `relwithdebinfo`, `tests`, `pgo-train`, `pgo-use`.

Default platform is `cuda` unless otherwise prompted. `cuda-tests` does not build the `.so` library, so you cannot run `python/robosample/run.py` or any other driver.

## Configure

```bash
cmake --preset ${CONFIG}
```

## Build

```bash
cmake --build --preset ${CONFIG}
```

## Validation levels

Every change SHALL either pass the tests and physical invariants stated here or depart from them deliberately, with justification. A test encodes why the behavior matters; "tests pass" is false if any were skipped. Tests can run for a long time, so prompt the user before running them.

- Level 0: compile only
- Level 1: compile and run basic example

    ```bash
    python3 python/robosample/run.py --name ala-dipeptide --prmtop examples/ala-dipeptide.prmtop --inprcrd examples/ala-dipeptide.rst7 --seed 6000 --equil_steps 0 --prod_steps 10 --write_freq 1 --validate true
    ```

- Level 2: compile, run basic example and run test suite

    ```bash
    nox -s tests
    ```

Default validation is level 1. Request confirmation before escalating to higher levels.

## Workflow

Workflow is a state machine: Issue -> Research -> Specification -> Implementation -> Review -> Revision -> Merge. Each state defines at least its inputs, outputs, and exit criteria. Examples:

- Research:
  - Output: accepted specification
  - Exit: no unresolved scientific questions
- Implementation:
  - Exit: builds successfully, requested validation completed
- Review:
  - Exit: no confirmed Blocking findings (severity Critical or High)

## Severity

Severity is determined by impact, not effort required to fix.

- Critical:
  - Produces incorrect probability distributions.
  - Violates detailed balance.
  - Violates physical invariants.
  - Prevents successful build or execution.
- High:
  - Produces incorrect behavior without invalidating the sampling algorithm.
  - Introduces regressions.
  - Breaks public APIs.
- Medium:
  - Maintainability.
  - Documentation.
  - Testing gaps.
- Low:
  - Formatting.
  - Naming.
  - Style.

The reviewer emits findings graded Blocking, Should-fix, or Nit. Grade maps to severity:

| Reviewer grade | Severity |
| --- | --- |
| Blocking | Critical or High |
| Should-fix | Medium |
| Nit | Low |

The reviewer's evidence tier (Reproducer, Re-derivation, Suppress) is orthogonal to severity. It governs whether a finding counts as *confirmed* for the Review exit criterion. A finding blocks merge only when it is both Blocking and confirmed. NOTE: the tier semantics are defined in `reviewer.md`; this mapping SHALL be cross-checked against that file.

A **specification** SHALL define, using these exact headings:

- **Motivation** - the problem in the codebase's vocabulary, quantified, with its constraints and assumptions. Researcher specs SHALL include the user's original phrasing alongside the restatement here.
- **Behavior** - what the system SHALL do, with the rationale that makes it correct. Physics specs SHALL include the claims and a derivation sketch with citations here.
- **Invariants** - the properties that hold after the change.
- **Interface** - the externally observable interface and the components and conventions that change.
- **Validation strategy** - the analytic or numeric checks that distinguish a correct implementation from a plausible-but-biased one, with oracles tagged PRECONDITION, INVARIANT, or LEMMA.

A specification MAY add:

- **Consequences and trade-offs** - SHOULD appear when the design forecloses alternatives.
- **Open questions** - SHALL appear when unresolved unknowns block a correct derivation; otherwise omitted.

These headings are canonical. `styles/spec.md` uses them directly. The researcher's artifact maps its domain-specific sections onto them as subsections:

| Canonical | Researcher subsection |
| --- | --- |
| Motivation | Problem restatement |
| Behavior | Claims, Derivation sketch |
| Invariants | Correctness conditions |
| Interface | Touch list |
| Validation strategy | Verification plan |
| Open questions | Open questions |

Implementation details SHOULD be omitted unless they are necessary to define externally observable behavior.

**Architectural** changes SHALL reference a decision record under `docs/decisions/`, creating one if none applies.

Review priority:

1. Physical correctness
2. Numerical correctness
3. API correctness
4. Test coverage
5. Performance
6. Style

## Evidence

Claims based on inference SHOULD be explicitly identified as such.

Behavioral claims SHALL cite one of:

- Published literature
- Project specification
- Existing regression test

Performance claims SHALL include reproducible benchmark results.

Numerical correctness claims SHALL include validation evidence against theory, regression tests, or reference implementations.

## Decision policy

When multiple sources define behavior, conflicts SHALL be resolved using the following precedence:

1. Published theory (see `references/index.yaml`).
2. Project specifications.
3. Existing tests.
4. Existing implementation.
5. Model intuition.

If in doubt, iteratively ask questions until the implementation plan is complete. After planning is complete, continue implementation without waiting for confirmation and stop only if blocked by missing information such as:

- Conflicting specifications.
- Missing external information.
- Ambiguous user intent.
- Missing repository files.

If supporting evidence cannot be located,

- State uncertainty explicitly.
- Do not invent references.
- Distinguish inference from established behavior.

## Engineering principles

Prefer localized changes.

Prefer correctness over optimization.

Prefer explicitness over cleverness.

Prefer local reasoning over global abstractions.

Preserve existing APIs unless specifications change.

Refactoring SHALL preserve observable behavior.

Behavioral changes SHALL update specifications.

Avoid introducing dependencies without justification.

Optimize only after measurement.
