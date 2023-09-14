# noxfile.py

from pathlib import Path

import nox
from nox.sessions import Session

nox.options.sessions = "tests", "coverage"

locations = "src", "tests", "noxfile.py"


python_versions = ["3.10", "3.11"]


@nox.session(python=python_versions)
def tests(session: Session) -> None:
    session.install(".")
    session.install("coverage[toml]", "pytest", "pygments")
    try:
        session.run("coverage", "run", "--parallel", "-m", "pytest", *session.posargs)
    finally:
        if session.interactive:
            session.notify("coverage")


@nox.session(python=python_versions)
def coverage(session: Session) -> None:
    args = session.posargs or ["report"]

    session.install("coverage[toml]")

    if not session.posargs and any(Path().glob(".coverage.*")):
        session.run("coverage", "combine")

    session.run("coverage", *args)
