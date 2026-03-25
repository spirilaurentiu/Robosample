import nox
import subprocess
import os
from pathlib import Path

BUILD_DIR = Path("build")
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

@nox.session
def tests(session):
    session.install("pytest", "pytest-cov", "gcovr", "anybadge", "jq")

    session.run("bash", "-c", "find . -name '*.gcd*' -delete")
    session.run("bash", "-c", "rm -rf coverage")

    session.run(
        "pytest",
        "--cov=python/robosample/",
        "--cov-report=xml:coverage/python_coverage.xml",
    )

    session.run(
        "gcovr",
        "-r", ".",
        "--xml", "coverage/cpp_coverage.xml",
        "--gcov-ignore-parse-errors", "negative_hits.warn",
    )

    session.run(
        "gcovr",
        "--cobertura-add-tracefile", "coverage/python_coverage.xml",
        "--cobertura-add-tracefile", "coverage/cpp_coverage.xml",
        "--html-details", "coverage/index.html",
        "--json-summary-pretty",
        "-o", "coverage/summary.json",
    )

    session.run(
        "bash", "-c",
        "anybadge --value=$(jq '.line_percent' coverage/summary.json) "
        "--file=coverage.svg --label=Coverage --suffix='%' "
        "--overwrite 50=red 75=orange 90=yellow 102=green"
    )

@nox.session
def build_optimized(session):
    # Set up environment
    conda_prefix = os.environ.get("CONDA_PREFIX")

    # Profile data
    PROFILE_DIR.mkdir(exist_ok=True)
    so_path = next(SO_DIR.glob(PYBIND_SO_PATTERN))
    perf_data = PROFILE_DIR / "perf.data"
    fdata = PROFILE_DIR / "perf.fdata"

    # Local installation of our package
    session.run("pip", "install", "-e", ".", external=True)
    
    # Clean previous builds and profiles
    for path in BUILD_DIR.glob("cuda-pgo-*"):
        session.run("rm", "-rf", str(path))
    session.run("rm", "-rf", str(PROFILE_DIR))

    # Build with PGO instrumentation
    session.run("cmake", "--preset", "cuda-pgo-train", env={"CONDA_PREFIX": conda_prefix})
    session.run("cmake", "--build", "--preset", "cuda-pgo-train")

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
        str(WRITE_FREQ)
    )

    # Build optimized
    session.run("cmake", "--preset", "cuda-pgo-use", env={"CONDA_PREFIX": conda_prefix})
    session.run("cmake", "--build", "--preset", "cuda-pgo-use")

    # Profile
    session.run(
        "perf",
        "record",
        "-o", str(perf_data),
        "-e", "cycles:u",
        "-j", "any,u",
        "python3",
        "python/robosample/roborun.py",
        TARGET,
        PRMTOP,
        RST7,
        str(SEED),
        str(EQUIL_STEPS),
        str(PROD_STEPS),
        str(WRITE_FREQ)
    )

    # Convert perf data to BOLT-friendly format
    session.run(
        "perf2bolt",
        str(so_path),
        "-p", str(perf_data),
        "-o", str(fdata),
    )

    # Optimize with BOLT
    session.run(
        "llvm-bolt",
        str(so_path),
        "-o", BOLT_SO,
        f"-data={fdata}",
        "-reorder-blocks=ext-tsp",
        "-reorder-functions=hfsort",
        "-dyno-stats",
    )

    # Replace original SO
    session.run("mv", "-f", BOLT_SO, str(so_path), external=True)
    