import csv
import os
import pathlib
import re
import shutil
import time
from pathlib import Path

import nox
from nox.command import CommandFailed

BUILD_DIR = Path("build")

TEST_PRESET = "cuda-tests"
BUILD_RELWITHDEBINFO = BUILD_DIR / TEST_PRESET

PROFILE_DIR = Path("profile-data")
SO_DIR = Path("python/robosample")
PYBIND_SO_PATTERN = "robo_bindings*.so"
BOLT_SO = "robo_bindings.bolt.so"

# ---- `nox -s profile` configuration -----------------------------------------
# Optimized-with-symbols build: -O3 -g -fno-omit-frame-pointer (see CMakeLists).
# perf/nsys need the line info this preset carries; cuda-release strips it.
PROFILE_PRESET = "cuda-relwithdebinfo"
PROFILE_DRIVER = SO_DIR / "prof_ffar1.py"
# The workload the user asked to profile: add_robotic_world 5 fs psi/psi torsion
# + external free DOFs on FFAR1, 0 equilibration / 100 production rounds.
PROFILE_SYSTEM = ("ffar1", "examples/ffar1.prmtop", "examples/ffar1.rst7")
PROFILE_EQUIL = 0
PROFILE_PROD = 100
# One DCD/CSV frame every 20 rounds -> 5 frames over a 100-round run. Enough to
# exercise the writer/position-gather path (~5% duty) without file I/O
# dominating a short profile. Set to PROFILE_PROD for a near-zero-I/O compute
# profile; lower it to mirror a production write cadence.
PROFILE_WRITE_FREQ = 20
PROFILE_NCU_LAUNCHES = 30  # ncu replays kernels serially; cap or it never ends.

SEED = 6000
EQUIL_STEPS = 0
PROD_STEPS = 1000
WRITE_FREQ = 1

TEST_SYSTEMS = [
    (
        "ala-dipeptide",
        "examples/ala-dipeptide.prmtop",
        "examples/ala-dipeptide.rst7",
        20,
    ),
    (
        "1APQ",
        "examples/1APQ.prmtop",
        "examples/1APQ.rst7",
        20,
    ),
    (
        "ffar1",
        "examples/ffar1.prmtop",
        "examples/ffar1.rst7",
        20,
    ),
    (
        "GfcDstrippedMin",
        "examples/GfcDstrippedMin.prmtop",
        "examples/GfcDstrippedMin.rst7",
        20,
    ),
]


@nox.session(python=False)
def tests(session):
    """Run CMake builds, pytest with xdist, and generate unified coverage."""
    if "CONDA_PREFIX" not in os.environ:
        session.error(
            "CONDA_PREFIX not found. Please activate your mamba environment first."
        )
    else:
        session.log(f"Using CONDA_PREFIX: {os.environ['CONDA_PREFIX']}")

    session.env["ROBOSAMPLE_SLOW_TESTS"] = "1"
    # The authoritative gate must FAIL loudly on missing oracles (OpenMM/ParmEd/
    # the compiled extension/2ala inputs) rather than silently skip them -- a
    # bare dev run may still skip (test_openmm_potential_energy.py checks this).
    session.env["ROBOSAMPLE_REQUIRE_OPENMM"] = "1"

    session.log("Cleaning old coverage data...")
    if os.path.exists("coverage"):
        shutil.rmtree("coverage")
    os.makedirs("coverage", exist_ok=True)
    shutil.rmtree(BUILD_RELWITHDEBINFO)

    session.log("Configuring and Building with CMake...")
    session.run("cmake", "--preset", TEST_PRESET)
    session.run("cmake", "--build", "--preset", TEST_PRESET)

    tests_failed = False

    # --- C++ tests: exit 8 = tests failed (tolerated). Any OTHER nonzero is a
    # real ctest error and should still abort. Detect "tests failed" from the
    # log ctest writes; the build dir was just wiped, so it can't be stale.
    session.log("Running C++ tests with CTest...")
    # -LE openmm: never run OpenMM's own third-party suite here. It is off by
    # default (built only under -DBUILD_OPENMM_TESTS=ON) and has its own opt-in
    # session, `nox -s openmm_tests`. The exclude is belt-and-suspenders in case
    # a build dir was configured with the option on.
    session.run(
        "ctest", "--test-dir", BUILD_RELWITHDEBINFO, "-j", "-LE", "openmm",
        success_codes=[0, 8],
    )
    if (
        BUILD_RELWITHDEBINFO / "Testing" / "Temporary" / "LastTestsFailed.log"
    ).exists():
        tests_failed = True
        session.warn("Some C++ tests failed (continuing to coverage).")

    # --- Python tests: catch the failure so coverage still runs. A pytest
    # error (usage/internal) also lands here and rightly keeps the session red.
    session.log("Running Python tests in parallel...")
    session.env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
    # Under the authoritative gate (ROBOSAMPLE_REQUIRE_OPENMM=1), exit 5 ("no
    # tests collected") must NOT be silently tolerated -- that is exactly the
    # "ran far fewer oracles and reported green" failure mode this gate exists
    # to catch. A bare dev run (var unset) keeps tolerating it.
    pytest_success_codes = [0] if session.env.get("ROBOSAMPLE_REQUIRE_OPENMM") else [0, 5]
    try:
        session.run(
            "pytest",
            "-p",
            "xdist",
            "-p",
            "pytest_cov",
            "-n",
            "auto",
            "--cov=python/robosample/",
            "--cov-report=xml:coverage/python_coverage.xml",
            *session.posargs,
            success_codes=pytest_success_codes,
        )  # 5 = no tests collected, treated as ok ONLY on a bare (non-gate) run
    except CommandFailed:
        tests_failed = True
        session.warn("Some Python tests failed (continuing to coverage).")

    # --- Coverage + badge: now always reached. ---
    session.log("Processing C++ coverage...")
    session.run(
        "gcovr",
        "-j",
        "0",
        "--exclude",
        r"tests/",
        "--xml",
        "coverage/cpp_coverage.xml",
        "--gcov-ignore-parse-errors",
        "negative_hits.warn",
        "--gcov-ignore-parse-errors",
        "suspicious_hits.warn",
        BUILD_RELWITHDEBINFO,
    )

    session.log("Merging coverage and generating HTML...")
    session.run(
        "gcovr",
        "--cobertura-add-tracefile",
        "coverage/python_coverage.xml",
        "--cobertura-add-tracefile",
        "coverage/cpp_coverage.xml",
        "--html-details",
        "coverage/index.html",
        "--json-summary-pretty",
        "-o",
        "coverage/summary.json",
    )

    session.log("Generating coverage badge...")
    badge_cmd = (
        "anybadge --value=$(jq '.line_percent' coverage/summary.json) "
        "--file=coverage.svg --label=Coverage --suffix='%' "
        "--overwrite 50=red 75=orange 90=yellow 102=green"
    )
    session.run("bash", "-c", badge_cmd, external=True)

    if tests_failed:
        session.error(
            "Test failures occurred — coverage and badge were still generated above."
        )


@nox.session(python=False)
def openmm_tests(session):
    """Run OpenMM's own test suite (opt-in; OFF by default, like slow tests).

    OpenMM's tests are third-party. We run them only when we have modified
    vendored OpenMM, to catch regressions. They are excluded from `nox -s tests`
    and are not even compiled unless BUILD_OPENMM_TESTS is set -- this session
    configures the Tests preset with -DBUILD_OPENMM_TESTS=ON, builds, and runs
    the "openmm" ctest label.

    The tests are SMALL systems (<= ~10k particles), so many share one GPU: each
    GPU test runs in single, mixed, AND double precision as independent, parallel
    ctest cases. -j is bounded (default 8) so VRAM is not oversubscribed; override
    with `--`, e.g.  nox -s openmm_tests -- -j 16 -R Cuda_Nonbonded.
    """
    if "CONDA_PREFIX" not in os.environ:
        session.error(
            "CONDA_PREFIX not found. Please activate your mamba environment first."
        )

    session.log("Configuring and building the Tests preset with OpenMM tests ON...")
    session.run("cmake", "--preset", TEST_PRESET, "-DBUILD_OPENMM_TESTS=ON")
    session.run("cmake", "--build", "--preset", TEST_PRESET)

    # Default to a bounded -j so the small GPU tests share one GPU without
    # oversubscribing VRAM; the user can override the whole ctest arg list via
    # posargs (e.g. a higher -j, a -R name filter, or --output-on-failure).
    ctest_args = list(session.posargs) or ["-j", "8", "--output-on-failure"]
    session.run(
        "ctest", "--test-dir", BUILD_RELWITHDEBINFO, "-L", "openmm", *ctest_args,
        success_codes=[0, 8],
    )


def is_perf_unrestricted():
    """Checks if kernel.perf_event_paranoid is set to -1."""
    path = pathlib.Path("/proc/sys/kernel/perf_event_paranoid")
    try:
        value = path.read_text().strip()
        return value == "-1"
    except FileNotFoundError:
        return False


@nox.session(reuse_venv=True)
def build_optimized(session):
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix is None:
        session.error(
            "CONDA_PREFIX is not set. Please activate your conda environment "
            "before running this nox session."
        )

    if not is_perf_unrestricted():
        session.error(
            "\n`perf_event_paranoid` is not -1.\n"
            "Please run the following command and try again:\n"
            "    sudo sysctl -w kernel.perf_event_paranoid=-1"
        )

    # Clean previous builds and profiles
    for path in BUILD_DIR.glob("cuda-pgo-*"):
        shutil.rmtree(path)
    shutil.rmtree(PROFILE_DIR, ignore_errors=True)

    # Profile data
    PROFILE_DIR.mkdir(exist_ok=True)
    so_path = next(SO_DIR.glob(PYBIND_SO_PATTERN))
    perf_data = PROFILE_DIR / "perf.data"
    # perf_lbr_data = PROFILE_DIR / "perf.lbr.data"
    fdata = PROFILE_DIR / "perf.fdata"

    # Build with PGO instrumentation
    session.run(
        "cmake",
        "--preset",
        "cuda-pgo-train",
        env={"CONDA_PREFIX": conda_prefix},
        external=True,
    )
    session.run(
        "cmake",
        "--build",
        "--preset",
        "cuda-pgo-train",
        env={"CONDA_PREFIX": conda_prefix},
        external=True,
    )

    # Local installation of our package
    session.run(
        "pip", "install", "-e", ".", env={"CONDA_PREFIX": conda_prefix}, external=True
    )

    # Run PGO on the smallest target to generate profile data
    session.run(
        "python3",
        "python/robosample/roborun.py",
        "ala-dipeptide",
        "examples/ala-dipeptide.prmtop",
        "examples/ala-dipeptide.rst7",
        str(SEED),
        str(EQUIL_STEPS),
        str(PROD_STEPS),
        str(WRITE_FREQ),
        env={"CONDA_PREFIX": conda_prefix},
        external=True,
    )

    # Build optimized
    session.run(
        "cmake",
        "--preset",
        "cuda-pgo-use",
        env={"CONDA_PREFIX": conda_prefix},
        external=True,
    )
    session.run(
        "cmake",
        "--build",
        "--preset",
        "cuda-pgo-use",
        env={"CONDA_PREFIX": conda_prefix},
        external=True,
    )

    # Profile with LBR enabled
    session.run(
        "perf",
        "record",
        "-o",
        str(perf_data),
        "-e",
        "cycles:u",
        "-j",
        "any,u",
        # "--call-graph",
        # "fp",
        "python3",
        "python/robosample/roborun.py",
        "ala-dipeptide",
        "examples/ala-dipeptide.prmtop",
        "examples/ala-dipeptide.rst7",
        str(SEED),
        str(EQUIL_STEPS),
        str(PROD_STEPS),
        str(WRITE_FREQ),
        env={"CONDA_PREFIX": conda_prefix},
        external=True,
    )

    # # Convert to LBR-augmented profile
    # session.run(
    #     "perf",
    #     "inject",
    #     "-j",
    #     "-i",
    #     str(perf_data),
    #     "-o",
    #     str(perf_lbr_data),
    #     env={"CONDA_PREFIX": conda_prefix},
    #     external=True,
    # )

    # Generate BOLT profile
    session.run(
        "perf2bolt",
        str(so_path),
        "-p",
        str(perf_data),  # perf_lbr_data
        "-o",
        str(fdata),
        "-nl",
        env={"CONDA_PREFIX": conda_prefix},
        external=True,
    )

    # Optimize with BOLT
    session.run(
        "llvm-bolt",
        str(so_path),
        "-o",
        BOLT_SO,
        f"-data={fdata}",
        "-reorder-blocks=ext-tsp",
        "-reorder-functions=hfsort+",  # hfsort+ is generally better than hfsort for large binaries
        "-split-functions",  # split hot/cold regions within functions
        "-split-eh",  # split exception handling code to cold
        "-dyno-stats",
        "-eliminate-unreachable",  # remove dead code BOLT can identify from the profile
        env={"CONDA_PREFIX": conda_prefix},
        external=True,
    )

    # Replace original SO
    session.run("mv", "-f", BOLT_SO, str(so_path), external=True)


# ---------------------------------------------------------------------------
# Profiling: three-layer capture (nsys -> perf -> ncu) into an agent-ingestible
# bundle under profile-data/<name>/. Layers, and why in this order:
#   1. Nsight Systems (nsys) -- CPU<->GPU timeline. Answers the FIRST question:
#      CPU-bound, GPU-kernel-bound, or transfer/sync-bound? perf is blind to the
#      GPU and to the getState(Forces) D2H copy that optimizer.md flags as the
#      top target, so nsys runs first and decides whether the other layers matter.
#   2. perf record -- host-CPU sampling attributed to C++ symbols (the ABA spine,
#      the force->body reduction, the pybind11 boundary). LBR form the noxfile
#      already standardizes on (-e cycles:u -j any,u), plus fp call graphs.
#   3. Nsight Compute (ncu) -- per-kernel counters. OPT-IN only (`--ncu`): it
#      replays every kernel serially (huge slowdown) and needs GPU-counter
#      permissions. Rarely your lever here -- all CUDA is OpenMM's, not ours.
# ---------------------------------------------------------------------------


def _csv_total_ns(path):
    """Sum the 'Total Time (ns)' column of an nsys stats CSV; ([], 0) on trouble.

    Returns (rows, total_ns) where rows is a list of (name, ns) sorted desc.
    Guarded: an nsys column-name change degrades to empty, never crashes the run.
    """
    if not path.exists():
        return [], 0
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            time_key = next(
                (k for k in (reader.fieldnames or []) if "Total Time" in k), None
            )
            name_key = next(
                (
                    k
                    for k in (reader.fieldnames or [])
                    if k in ("Name", "Operation", "Category")
                ),
                None,
            )
            if time_key is None:
                return [], 0
            rows = []
            for r in reader:
                try:
                    ns = float(r[time_key].replace(",", ""))
                except (ValueError, AttributeError, KeyError):
                    continue
                rows.append((r.get(name_key, "?") if name_key else "?", ns))
    except OSError:
        return [], 0
    rows.sort(key=lambda t: t[1], reverse=True)
    return rows, sum(ns for _, ns in rows)


def _write_manifest(session, out_dir, name, invocation, gpu_name, wall_s, ran_ncu):
    """Emit MANIFEST.md: the top-line an agent reads before the raw CSVs."""
    kern_rows, kern_ns = _csv_total_ns(out_dir / "nsys_cuda_gpu_kern_sum.csv")
    memcpy_rows, memcpy_ns = _csv_total_ns(out_dir / "nsys_cuda_gpu_mem_time_sum.csv")

    def us(ns):
        return f"{ns / 1e3:,.1f} us"

    def top(rows, n=5):
        return (
            "\n".join(f"  - {us(ns):>16}  {nm}" for nm, ns in rows[:n]) or "  - (none)"
        )

    gpu_ns = kern_ns + memcpy_ns
    lines = [
        f"# Profile manifest -- {name}",
        "",
        f"- generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        f"- gpu: {gpu_name or 'unknown'}",
        f"- build preset: {PROFILE_PRESET} (-O3 -g -fno-omit-frame-pointer)",
        f"- driver: {PROFILE_DRIVER}",
        f"- invocation: `{invocation}`",
        f"- wall clock (nsys run, driver-reported): "
        + (f"{wall_s:.2f} s" if wall_s else "unparsed"),
        "",
        "## Top-line split (GPU, from nsys)",
        "",
        f"- GPU kernel time (sum):  {us(kern_ns)}",
        f"- GPU memcpy time (sum):  {us(memcpy_ns)}"
        + (
            f"   ({100 * memcpy_ns / gpu_ns:.0f}% of on-GPU time -- the D2H/H2D bridge)"
            if gpu_ns
            else ""
        ),
        "",
        "  NOTE: host-CPU time is wall_clock - (overlapped GPU time); read "
        "perf.symbols.txt for its attribution. If memcpy dominates GPU time, the "
        "lever is fusing the force reduction on-device (optimizer.md, CUDA section).",
        "",
        "## Top CUDA kernels",
        "",
        top(kern_rows),
        "",
        "## Top CUDA memory ops",
        "",
        top(memcpy_rows),
        "",
        "## Files in this bundle",
        "",
        "- nsys.nsys-rep                  full Nsight Systems timeline (open in GUI)",
        "- nsys_cuda_gpu_kern_sum.csv     per-kernel GPU time",
        "- nsys_cuda_gpu_mem_time_sum.csv memcpy time by direction (the transfer story)",
        "- nsys_cuda_gpu_mem_size_sum.csv bytes moved by direction",
        "- nsys_cuda_api_sum.csv          host-side CUDA API time",
        "- nsys_osrt_sum.csv              OS-runtime (syscall/IO) time",
        "- perf.data                      raw perf sampling (open: hotspot / perf report)",
        "- perf.symbols.txt               flat host symbol table, >=0.5% samples",
        "- perf.callgraph.txt             host call-graph, >=0.5%",
    ]
    if ran_ncu:
        lines.append("- ncu.csv                        per-kernel deep counters")
    lines += [
        "",
        "## How to consume",
        "",
        "Hand this directory to the `optimizer` agent. It confirms the Amdahl "
        "ceiling from these numbers before touching source (optimizer.md, "
        "'Profile before touching anything').",
        "",
    ]
    (out_dir / "MANIFEST.md").write_text("\n".join(lines))
    session.log(f"Wrote {out_dir / 'MANIFEST.md'}")


@nox.session(python=False)
def profile(session):
    """Three-layer profile (nsys/perf[/ncu]) -> profile-data/<name>/ for an agent.

    Usage:
      nox -s profile                 # nsys + perf on the FFAR1 workload
      nox -s profile -- --ncu        # also run Nsight Compute (slow, needs perms)
      nox -s profile -- --build      # configure+build cuda-relwithdebinfo first
    """
    if "CONDA_PREFIX" not in os.environ:
        session.error("CONDA_PREFIX is not set. Activate the build/runtime conda env.")

    do_ncu = "--ncu" in session.posargs
    do_build = "--build" in session.posargs

    if do_build:
        session.log(f"Configuring and building {PROFILE_PRESET}...")
        session.run("cmake", "--preset", PROFILE_PRESET, external=True)
        session.run("cmake", "--build", "--preset", PROFILE_PRESET, external=True)

    so = next(SO_DIR.glob(PYBIND_SO_PATTERN), None)
    if so is None:
        session.error(
            f"No {PYBIND_SO_PATTERN} in {SO_DIR}. Build it first:\n"
            f"    cmake --preset {PROFILE_PRESET}\n"
            f"    cmake --build --preset {PROFILE_PRESET}\n"
            f"or re-run: nox -s profile -- --build"
        )

    name, prmtop, inpcrd = PROFILE_SYSTEM
    out_dir = PROFILE_DIR / name
    shutil.rmtree(out_dir, ignore_errors=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    driver_args = [
        str(PROFILE_DRIVER),
        name,
        prmtop,
        inpcrd,
        str(SEED),
        str(PROFILE_EQUIL),
        str(PROFILE_PROD),
        str(PROFILE_WRITE_FREQ),
    ]
    invocation = "python3 " + " ".join(driver_args)
    session.log(f"Workload: {invocation}")

    gpu_name = None
    try:
        gpu_name = (
            session.run(
                "nvidia-smi",
                "--query-gpu=name",
                "--format=csv,noheader",
                silent=True,
                external=True,
            )
            or ""
        ).strip()
    except CommandFailed:
        pass

    # --- Layer 1: Nsight Systems (also our lowest-overhead wall-clock source) --
    session.log("[1/3] nsys: CPU<->GPU timeline...")
    nsys_base = out_dir / "nsys"
    nsys_out = ""
    try:
        nsys_out = (
            session.run(
                "nsys",
                "profile",
                "-o",
                str(nsys_base),
                "--force-overwrite=true",
                "-t",
                "cuda,nvtx,osrt",
                "--sample=cpu",
                "--cuda-memory-usage=true",
                "python3",
                *driver_args,
                silent=True,
                external=True,
            )
            or ""
        )
        (out_dir / "run.log").write_text(nsys_out)
        for report in (
            "cuda_gpu_kern_sum",
            "cuda_gpu_mem_time_sum",
            "cuda_gpu_mem_size_sum",
            "cuda_api_sum",
            "osrt_sum",
        ):
            session.run(
                "nsys",
                "stats",
                "--report",
                report,
                "--format",
                "csv",
                "--output",
                str(nsys_base),
                f"{nsys_base}.nsys-rep",
                external=True,
                success_codes=[0, 1],
            )
    except CommandFailed:
        session.warn("nsys failed -- is Nsight Systems on PATH? Continuing to perf.")

    m = re.search(r"run_rex\(\) took ([\d.]+) seconds", nsys_out)
    wall_s = float(m.group(1)) if m else None

    # --- Layer 2: perf record (host-CPU attribution) --------------------------
    session.log("[2/3] perf: host-CPU sampling...")
    perf_data = out_dir / "perf.data"
    try:
        session.run(
            "perf",
            "record",
            "-F",
            "999",
            "-e",
            "cycles:u",
            "-j",
            "any,u",
            "-g",
            "--call-graph",
            "fp",
            "-o",
            str(perf_data),
            "--",
            "python3",
            *driver_args,
            external=True,
            env={"PYTHONPERFSUPPORT": "1"},
        )
        symbols = session.run(
            "perf", "report", "-i", str(perf_data), "--stdio", "-n",
            "--sort=overhead,symbol", "--percent-limit", "0.5",
            silent=True, external=True, success_codes=[0, 1],
        )
        (out_dir / "perf.symbols.txt").write_text(symbols or "")
        callgraph = session.run(
            "perf", "report", "-i", str(perf_data), "--stdio", "-n", "-g",
            "graph,0.5,caller",
            silent=True, external=True, success_codes=[0, 1],
        )
        (out_dir / "perf.callgraph.txt").write_text(callgraph or "")
    except CommandFailed:
        session.warn(
            "perf failed. LBR (-j any,u) may need a lower perf_event_paranoid:\n"
            "    sudo sysctl -w kernel.perf_event_paranoid=1"
        )

    # --- Layer 3: Nsight Compute (opt-in; slow; needs GPU-counter perms) ------
    if do_ncu:
        session.log(f"[3/3] ncu: per-kernel counters (<= {PROFILE_NCU_LAUNCHES} launches)...")
        try:
            session.run(
                "ncu",
                "--set",
                "full",
                "--launch-count",
                str(PROFILE_NCU_LAUNCHES),
                "--target-processes",
                "all",
                "--csv",
                "--log-file",
                str(out_dir / "ncu.csv"),
                "python3",
                *driver_args,
                external=True,
            )
        except CommandFailed:
            session.warn(
                "ncu failed. If ERR_NVGPUCTRPERM, GPU perf counters are admin-gated;\n"
                "a human must set the nvidia module option "
                "NVreg_RestrictProfilingToAdminUsers=0 (sudo)."
            )
    else:
        session.log("[3/3] ncu: skipped (pass -- --ncu to enable).")

    _write_manifest(session, out_dir, name, invocation, gpu_name, wall_s, do_ncu)
    session.log(f"Profile bundle ready: {out_dir}/  (start with MANIFEST.md)")
