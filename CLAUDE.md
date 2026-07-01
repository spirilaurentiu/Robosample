# Robosample

Robosample is a molecular simulation software implementing blocked Gibbs Sampling coupled with Hamiltonian Monte Carlo [load paper]. Each molecule is represented as a **robot** (kinematic tree of rigid bodies connected by joints, or *mobilizers*) and each Gibbs block consists of the set of joints of a forest a forest of kinematic trees rooted at a shared single ground frame. Each robot is based on the previous robot i.e. Cartesian coordinates from the previous robot are assigned to the next robot. **Per-atom** forces computed by OpenMM are reduced to per-body spatial forces (net force + torque about the body origin) for the articulated-body solver.

Robosample is a `conda` package that ships:

- A molecular simulation engine written in C++17/CUDA and compiled compiled with PyBind11 as an `.so`.
- A Python library built on top of it.

The build toolchain (`cmake`, `ninja`, `gcc`, `cuda` etc) lives in its own `conda` environment which the user must create and activate. Before running build commands, check if running any environment inside `envs/robo_cuda*.yaml` or `envs/robo_cpu.yaml`.

The goal is to accelerate sampling of molecular conformations, maximizing sampling of **rare events** (basin hopping). Correctness is measured against the theory, not against intuition.

The systems it is supposed to simulate are to **1M atoms** clustered in up to **100k rigid bodies**. Examples:

- Alanine dipeptide in vacuum
- Deca alanine in explicit solvent
- FFAR1 (GPCR) in implicit solvent and with explicit membrane nanodisc
- Spliceosome in explicit solvent

## Agents and references

Each agent's authoritative behavior lives in its own frontmatter; the line below is the routing gloss -
what it does, where it sits in the pipeline, and its binding constraint.

- `.claude/agents/researcher.md` - **read-only.** Theory / sampling research -> a precise spec under `docs/specs/`. Runs first, before any physics-touching change. Never writes code.
- `.claude/agents/coder.md` - implements a reviewed spec or a self-contained task. **Auto-mode, surgical** (Rules 2, 3, 7). Correctness over performance - never optimizes speculatively.
- `.claude/agents/reviewer.md` - **independent, hostile, read-only.** Runs after the coder, before merge. Reviews science + conventions + implementation together; hands confirmed findings back to `coder`, never patches.
- `.claude/agents/optimizer.md` - **opt-in, user-invoked only; never inside a feature loop.** Static-analysis-guided source optimization measured with `perf`; hands every change to `reviewer`.

- `references/index.yaml`: summary of all papers Robosample is based on with link to their Markdown versions stored locally.

## How to run shell commands (mandatory)

Permissions are **generous-allow + hard-deny**: almost every command runs without a prompt; a short
deny-list in `.claude/settings.json` (sudo, `rm -rf`, history/remote-mutating git, network tools,
secret reads) is the real boundary and cannot be overridden. So a prompt or block now means one of two
things - you tried a genuinely denied operation, or you wrote the command in a form the matcher can't
analyze. Never route around a block by asking to widen permissions; rewrite the command.

Claude Code splits compound commands on `&&`, `||`, `;`, `|`, and newlines and matches each piece
independently, and it auto-runs bare read-only commands (`ls`, `cat`, `grep`, `find`, `head`, `tail`,
`git log/status/diff/show`, ...). Write commands so they stay in that clean, analyzable form:

- **No `cd` prefixes.** Pass explicit paths: `ctest --test-dir build/cuda-tests -R <name>`, `cmake --preset <p>`, `grep -rn 'pat' tests/ src/`, `sed -n '1,40p' path/to/file`. A `cd X; ...` prefix trips the path-resolution guard and buys nothing.
- **Run git from the repo root, without `-C`.** Plain `git log -- references/`, `git status --short`, `git diff ...` - these hit the built-in read-only-git safelist and never prompt. `git -C /path log ...` does not match that safelist and will prompt; the `-C` is redundant since the repo is already the `cwd`.
- **No command substitution.** Never `$(...)` or backticks. For parallelism use `-j0` (`cmake` / `ninja` / `ctest` read 0 as *all cores* - maximum parallelism, the project default). Not `-j$(nproc)`.
- **No heredocs and no shell loops.** Never `python3 - <<'PY' ... PY`, never `for ...; do ... done`. Write a script to a file and run `python3 file.py`.
  For a file-existence sweep use one command, not a loop: `find references/papers -maxdepth 2 -name meta.yaml` or `ls references/papers/*/meta.yaml`.
  For a multi-pattern search use one `grep -rn 'a\|b\|c' <dir>` (or `rg`).
- **Avoid chaining unrelated commands with `;`/`&&`/`|`.** Each piece must match on its own, and a long chain is where analysis breaks. Prefer one command per step; when you must combine, keep every piece a plain allowed form.
- **Environment variables via `env`, not a bare prefix.** `env ROBOSAMPLE_SLOW_TESTS=1 ctest ...`, not `ROBOSAMPLE_SLOW_TESTS=1 ctest ...` - a bare `VAR=val` prefix stops the command matching its rule.
- **Run tests through `ctest`, not the binaries.** `ctest --test-dir build/cuda-tests -R <name> --output-on-failure`, not `./TestFoo`.

**Do not approve slips into `settings.local.json`.** Approving a blocked command writes a frozen exact string that never generalizes - the next variant prompts again. A block is a signal to rewrite the command (or that it is genuinely denied), never to allowlist it. Any existing `.claude/settings.local.json` accumulated this way should be deleted (it is git-ignored by default; if you created it yourself, add it to `.gitignore`).

## The test gate (non-negotiable)

Every change must either pass the suite of tests and physical invariats stated here or change **deliberately and with justification**. A test encodes WHY behavior matters; "tests pass" is false if any were skipped.

- Fast smoke (<1 min, use while iterating): `python3 python/robosample/roborun.py ala-dipeptide examples/ala-dipeptide.prmtop examples/ala-dipeptide.rst7 6000 0 100 1 true`
- Full authoritative gate (before declaring done): `nox -s tests` - will build `cuda-tests`.
- A `mamba`/`conda` env must be active (`CONDA_PREFIX` is required by both the presets and `noxfile`).
- **A red `tests` session is a failed gate even though coverage/badge were still produced.** The session tolerates ctest exit 8 (tests-failed) and pytest exit 5 (no tests collected) so coverage still runs, then `session.error`s at the end if anything failed. Judge the gate by the session's final status and `LastTestsFailed.log`, never by a raw exit code mid-run - a suite that failed still emits artifacts.

## Build (CUDA)

Can build `${CONFIG}` as production (`cuda-release` - no tests) or testing (`cuda-tests` - no `.so`/Python built).

- Configure: `cmake --preset ${CONFIG}$`
- Build & install into `python/robosample/`: `cmake --build --preset ${CONFIG}`
- Builds run **locally** on this machine's GPU (`CMAKE_CUDA_ARCHITECTURES=native`). This is why Remote Control (phone) works and Claude Code on the web does not - the cloud has no CUDA toolchain.

## Profiling and optimization

- Day-to-day optimization is **static-analysis-guided source optimization**, measured with `perf record` (`-e cycles:u -j any,u`), driven by the `optimizer` agent. This is the normal path.
- `nox -s build_optimized` is a **separate, opt-in PGO+BOLT facility** - not part of any feature or optimization loop. It:
  - **requires a human** to first run `sudo sysctl -w kernel.perf_event_paranoid=-1`; agents cannot (sudo is denied in `.claude/settings.json`), so an agent that needs it must stop and ask, not escalate;
  - **destructively overwrites** the installed `robo_bindings*.so` in `python/robosample/`;
  - currently trains PGO on the smallest system (`ala-dipeptide`) only - a **known stopgap, not a validated policy**; do not treat it as a convention to enforce or preserve.

## How work is delegated

Research and design are read-only and produce a spec.

Implementation is surgical (Rules 2, 3, 7).

Validation runs the gate.

Review is hostile (Rule 12) and read-only.

See `.claude/agents/`. Optimization (`optimizer`) is a **separate, opt-in** pipeline - never invoke it inside a feature loop.

If in doubt, iteratively ask questions until the implementation plan is complete. After that, enter auto mode and finish coding independently. Do not ask further questions. Everything is permitted (read `.claude/settings.json`).
