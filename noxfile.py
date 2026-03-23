import nox

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