import os
import pathlib
import shutil
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
    session.run("ctest", "--test-dir", BUILD_RELWITHDEBINFO, "-j", success_codes=[0, 8])
    if (
        BUILD_RELWITHDEBINFO / "Testing" / "Temporary" / "LastTestsFailed.log"
    ).exists():
        tests_failed = True
        session.warn("Some C++ tests failed (continuing to coverage).")

    # --- Python tests: catch the failure so coverage still runs. A pytest
    # error (usage/internal) also lands here and rightly keeps the session red.
    session.log("Running Python tests in parallel...")
    session.env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
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
            success_codes=[0, 5],
        )  # 5 = no tests collected, treat as ok
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
