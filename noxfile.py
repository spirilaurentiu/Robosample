import os
from pathlib import Path

import nox

BUILD_DIR = Path("build")

TEST_PRESET = "cuda-relwithdebinfo"
BUILD_RELWITHDEBINFO = BUILD_DIR / TEST_PRESET

PROFILE_DIR = Path("profile-data")
SO_DIR = Path("python/robosample")
PYBIND_SO_PATTERN = "robo_bindings*.so"
BOLT_SO = "robo_bindings.bolt.so"

TARGET = "1APQ"
PRMTOP = f"examples/{TARGET}.prmtop"
RST7 = f"examples/{TARGET}.rst7"

SEED = 6000
EQUIL_STEPS = 0
PROD_STEPS = 1000
WRITE_FREQ = 1


@nox.session(python=False)
def tests(session):
    """Run CMake builds, pytest with xdist, and generate unified coverage."""

    # 1. Verification: Ensure we are actually in the expected environment
    if "CONDA_PREFIX" not in os.environ:
        session.error(
            "CONDA_PREFIX not found. Please activate your mamba environment first."
        )
    else:
        session.log(f"Using CONDA_PREFIX: {os.environ['CONDA_PREFIX']}")

    conda_prefix = os.environ["CONDA_PREFIX"]
    gcov_exe = os.path.join(conda_prefix, "bin", "x86_64-conda-linux-gnu-gcov")

    # # 2. Build Process
    # session.log("Configuring and Building with CMake...")
    # session.run("cmake", "--preset", TEST_PRESET)
    # session.run("cmake", "--build", "--preset", TEST_PRESET)

    # # 3. Cleanup old coverage artifacts
    # session.log("Cleaning old coverage data...")
    # # Use python's shutil/os for faster/cleaner cleanup than spawning bash for rm
    # if os.path.exists("coverage"):
    #     shutil.rmtree("coverage")
    # os.makedirs("coverage", exist_ok=True)

    # session.run(
    #     "find", BUILD_RELWITHDEBINFO, "-name", "*.gcda", "-delete", external=True
    # )

    # # 4. Python Tests with Parallel Execution
    # session.log("Running Python tests in parallel...")
    # session.env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
    # session.run(
    #     "pytest",
    #     # "-p",
    #     # "xdist",
    #     "-p",
    #     "pytest_cov",
    #     # "-n",
    #     # "auto",
    #     "--cov=python/robosample/",
    #     "--cov-report=xml:coverage/python_coverage.xml",
    #     *session.posargs,  # Allows you to pass extra args to pytest via nox
    # )

    # 5. C++ / GCOV Coverage
    session.log("Processing C++ coverage...")
    session.run(
        "gcovr",
        "-j",
        "0",  # Use all cores
        "--gcov-executable",
        gcov_exe,
        "--xml",
        "coverage/cpp_coverage.xml",
        "--gcov-ignore-parse-errors",
        "negative_hits.warn",
        "--gcov-ignore-parse-errors",
        "suspicious_hits.warn",
        BUILD_RELWITHDEBINFO,
    )

    # 6. Merge & Report
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

    # 7. Badge Generation
    # We use a single bash string here to allow the $(jq ...) subshell to work
    session.log("Generating coverage badge...")
    badge_cmd = (
        "anybadge --value=$(jq '.line_percent' coverage/summary.json) "
        "--file=coverage.svg --label=Coverage --suffix='%' "
        "--overwrite 50=red 75=orange 90=yellow 102=green"
    )
    session.run("bash", "-c", badge_cmd, external=True)


@nox.session(reuse_venv=True)
def build_optimized(session):
    # Set up environment
    conda_prefix = os.environ.get("CONDA_PREFIX")

    # Profile data
    PROFILE_DIR.mkdir(exist_ok=True)
    so_path = next(SO_DIR.glob(PYBIND_SO_PATTERN))
    perf_data = PROFILE_DIR / "perf.data"
    fdata = PROFILE_DIR / "perf.fdata"

    # Clean previous builds and profiles
    for path in BUILD_DIR.glob("cuda-pgo-*"):
        session.run("rm", "-rf", str(path), external=True)
    session.run("rm", "-rf", str(PROFILE_DIR), external=True)

    # Build with PGO instrumentation
    session.run(
        "cmake", "--preset", "cuda-pgo-train", env={"CONDA_PREFIX": conda_prefix}
    )
    session.run("cmake", "--build", "--preset", "cuda-pgo-train")

    # Local installation of our package
    session.run("pip", "install", "-e", ".", external=True)

    # Run PGO
    session.run(
        "python3",
        "python/robosample/roborun.py",
        TARGET,
        PRMTOP,
        RST7,
        str(SEED),
        str(EQUIL_STEPS),
        str(PROD_STEPS),
        str(WRITE_FREQ),
    )

    # Build optimized
    session.run("cmake", "--preset", "cuda-pgo-use", env={"CONDA_PREFIX": conda_prefix})
    session.run("cmake", "--build", "--preset", "cuda-pgo-use")

    # Profile
    session.run(
        "perf",
        "record",
        "-o",
        str(perf_data),
        "-e",
        "cycles:u",
        "-j",
        "any,u",
        "python3",
        "python/robosample/roborun.py",
        TARGET,
        PRMTOP,
        RST7,
        str(SEED),
        str(EQUIL_STEPS),
        str(PROD_STEPS),
        str(WRITE_FREQ),
    )

    # Convert perf data to BOLT-friendly format
    session.run(
        "perf2bolt",
        str(so_path),
        "-p",
        str(perf_data),
        "-o",
        str(fdata),
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
    )

    # Replace original SO
    session.run("mv", "-f", BOLT_SO, str(so_path), external=True)
