import nox
import subprocess
from pathlib import Path

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

BUILD_DIR = Path("build")
PROFILE_DIR = Path("profile-data")
SO_DIR = Path("python/robosample")
PYBIND_SO_PATTERN = "robo_bindings*.so"
BOLT_SO = "robo_bindings.bolt.so"
TARGET = "2ala"
PRMTOP = f"examples/{TARGET}.prmtop"
RST7 = f"examples/{TARGET}.rst7"
STEPS = 6000
SEED = 0
REPS = 100
THREADS = 1

@nox.session
def clean(session):
    """Clean old build and profile data."""
    for path in BUILD_DIR.glob("cuda-pgo-*"):
        session.run("rm", "-rf", str(path))
    session.run("rm", "-rf", str(PROFILE_DIR))

@nox.session
def build_pgo_train(session):
    """CMake build with PGO training preset."""
    session.run("cmake", "--preset", "cuda-pgo-train")
    session.run("cmake", "--build", "--preset", "cuda-pgo-train")

@nox.session
def train_pgo(session):
    """Run Python workload to collect PGO data."""
    session.run(
        "python3",
        "python/robosample/roborun.py",
        TARGET,
        PRMTOP,
        RST7,
        str(STEPS),
        str(SEED),
        str(REPS),
        str(THREADS),
    )

@nox.session
def build_pgo_use(session):
    """CMake build using collected PGO data."""
    session.run("cmake", "--preset", "cuda-pgo-use")
    session.run("cmake", "--build", "--preset", "cuda-pgo-use")

@nox.session
def profile_perf(session):
    """Profile Python workload and generate perf data."""
    PROFILE_DIR.mkdir(exist_ok=True)
    so_path = next(SO_DIR.glob(PYBIND_SO_PATTERN))
    perf_data = PROFILE_DIR / "perf.data"
    fdata = PROFILE_DIR / "perf.fdata"

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
        str(STEPS),
        str(SEED),
        str(REPS),
        str(THREADS),
    )

    # Convert perf data to BOLT-friendly format
    session.run(
        "perf2bolt",
        str(so_path),
        "-p", str(perf_data),
        "-o", str(fdata),
    )

@nox.session
def run_bolt(session):
    """Run LLVM BOLT to optimize the pybind11 module."""
    so_path = next(SO_DIR.glob(PYBIND_SO_PATTERN))
    fdata = PROFILE_DIR / "perf.fdata"

    session.run(
        "llvm-bolt",
        str(so_path),
        "-o", BOLT_SO,
        f"-data={fdata}",
        "-reorder-blocks=ext-tsp",
        "-reorder-functions=hfsort",
        "-dyno-stats",
    )

    # Replace original SO atomically
    session.run("mv", "-f", BOLT_SO, str(so_path))

@nox.session
def full_pipeline(session):
    """Run the full workflow end-to-end."""
    session.notify("clean")
    session.notify("build_pgo_train")
    session.notify("train_pgo")
    session.notify("build_pgo_use")
    session.notify("profile_perf")
    session.notify("run_bolt")
