# Tool contract (shared)

Agents express **intent**. Each intent below has a fixed name and post-condition; the
implementation (which nox session, which CMake preset, which clang binary) MAY change without
changing the agent that calls it. No agent writes raw shell for anything this file names.

All static tools consume `build/latest/compile_commands.json` (the symlink CMake maintains to the
active binary dir when `CMAKE_EXPORT_COMPILE_COMMANDS` is on). A tool that needs the compile
database SHALL run `configure(cfg)` first if `build/latest` is stale or absent; it MUST NOT require
a full build.

## Build / configure

- `configure(cfg)` -> emits `compile_commands.json` for `cfg`. Cheap. Precedes every static tool.
- `build(cfg)` -> compiles `cfg`. `cfg ∈ {sanitize, tsan, coverage, release}`.

The four configs map to four purposes and MUST NOT be conflated:

| cfg | build type | preset | instrumentation | builds `.so`? | registers C++ tests? | purpose |
|-----|-----------|--------|-----------------|---------------|----------------------|---------|
| `sanitize` | Tests | `<plat>-tests` | ASan + UBSan, `-O1 -g` | no | yes | memory / UB on fast C++ tests |
| `tsan` | Tsan (to add) | `<plat>-tsan` (to add) | TSan + UBSan | no | yes | data races (cannot share a binary with ASan) |
| `coverage` | Coverage (to add) | `<plat>-coverage` (to add) | gcov (`--coverage`), `-O0 -g` | yes | yes | line/branch coverage over the whole suite |
| `release` | Release | `<plat>-release` | none, `-O3` | yes | no | authoritative slow-test numbers + benchmarks |

NOTE: `sanitize` cannot run anything that imports the Python module — the Tests build skips
`pybind11_add_module` because a sanitized `.so` cannot be `dlopen`'d by a normal interpreter. Slow
tests that go through the driver run under `coverage` (for gcov) and `release` (for correctness),
never under `sanitize`.

NOTE: `coverage` is a build type that does not yet exist. RelWithDebInfo carries `--coverage` but
gates the C++ tests off (`BUILD_TESTING` is ON only for the Tests build type, CMakeLists §C++ Tests),
so it registers no gtest cases and cannot be the coverage carrier as-is. The `Coverage` build type
adds `--coverage` without sanitizers, registers the tests (extend the `BUILD_TESTING` guard to accept
`Coverage`), and builds the `.so` (a gcov-instrumented module loads under a normal interpreter). One
`coverage` build then holds the gtest suite, the driver, and pytest, all gcov-instrumented; gcovr
merges the C++ `.gcda` with pytest-cov's XML. `tsan` likewise needs a new build type and presets.

## Runtime

- `test(scope, cfg)` -> runs tests. `scope ∈ {fast, slow, openmm}`. `slow` includes the `run.py`
  driver on the test systems. `openmm` is opt-in (third-party suite, off by default).
- `run_single(name)` -> one case by regex.
- `run_driver(system)` -> `run.py` on one system; used under `coverage` so driver lines are counted,
  and under `release` for correctness.

`asan()` / `ubsan()` are `test(fast, sanitize)` — they are the config, not separate steps.
`tsan()` is `test(fast, tsan)`.

## Static analysis (all diff-scoped by default; consume the compile db)

- `format()` -> `clang-format -i`. Runs first; canonicalizes before any tool reads structure.
- `lint()` -> clang-tidy over the diff. No full build needed.
- `iwyu()` -> include-what-you-use over the diff.
- `analyze()` -> clang static analyzer. Branches off `configure`, parallel to the sanitizer run.
- `coverage()` -> `build(coverage)` + full suite + driver + gcovr merge (C++ gcov + pytest-cov).
- `benchmark(target)` -> google-benchmark target or the profile session, under `release`/PGO.

## Code intelligence (replace grep entirely)

Two layers, distinct failure modes, both required.

clangd index (needs a warm index):

- `find_definition(sym)`, `find_symbol(sym)`, `find_callers(sym)`, `find_callees(sym)`,
  `find_implementations(sym)`.

clang-query AST matchers (needs a correct matcher):

- `find_overrides(sym)`, `find_virtual_overrides(sym)`, `find_template_specializations(sym)`.

Graph (Doxygen + dot):

- `call_graph(sym|module)`, `dependency_graph(module)`.

## Handoff

- `documentation(symbols)` -> Doxygen contract for changed symbols. Authoring is delegated to the
  `documenter` agent; the caller passes scope, not prose.
- `commit(message)` -> git.

## Gate order

`format -> configure -> lint -> iwyu -> build(sanitize) -> test(fast, sanitize) -> build(tsan) ->
test(fast, tsan) -> [analyze, parallel] -> build(release) -> test(slow, release) -> coverage() ->
benchmark -> model review`.

Compilation is the first gate for anything that runs. `format`/`lint`/`iwyu`/`analyze` gate on
`configure` only. Model review runs only when every deterministic gate is green.
