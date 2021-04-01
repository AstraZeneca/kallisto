# noxfile.py

import os
import tempfile

import nox
from nox.sessions import Session

nox.options.sessions = "lint", "mypy", "pytype", "tests"

locations = "src", "tests", "noxfile.py"


python_versions = ["3.9", "3.8", "3.7"]


@nox.session(python=python_versions)
def tests(session: Session) -> None:
    args = session.posargs or ["--cov", "-m", "not e2e"]
    session.run("poetry", "install", "--no-dev", external=True)
    install_with_constraints(
        session, "coverage[toml]", "pytest", "pytest-cov",
    )
    session.run("pytest", *args)


@nox.session(python=python_versions)
def lint(session: Session) -> None:
    args = session.posargs or locations
    install_with_constraints(
        session,
        "flake8",
        "flake8-bandit",
        "flake8-black",
        "flake8-bugbear",
        "flake8-import-order",
    )
    session.run("flake8", *args)


@nox.session(python="3.8")
def black(session: Session) -> None:
    args = session.posargs or locations
    install_with_constraints(session, "black")
    session.run("black", *args)


@nox.session(python="3.8")
def safety(session: Session) -> None:
    with tempfile.TemporaryDirectory() as tmpdirname:
        session.run(
            "poetry",
            "export",
            "--dev",
            "--format=requirements.txt",
            "--without-hashes",
            f"--output={os.path.join(tmpdirname, 'poetry-export.out')}",
            external=True,
        )
        install_with_constraints(session, "safety")
        session.run(
            "safety",
            "check",
            f"--file={os.path.join(tmpdirname, 'safety-check.out')}",
            "--full-report",
        )


@nox.session(python=python_versions)
def mypy(session: Session) -> None:
    args = session.posargs or locations
    install_with_constraints(session, "mypy")
    session.run("mypy", *args)


@nox.session(python="3.7")
def pytype(session: Session) -> None:
    """Run the static type checker."""
    args = session.posargs or ["--disable=import-error", *locations]
    install_with_constraints(session, "pytype")
    session.run("pytype", *args)


@nox.session(python="3.8")
def coverage(session: Session) -> None:
    """Upload coverage data."""
    install_with_constraints(session, "coverage[toml]", "codecov")
    session.run("coverage", "xml", "--fail-under=0")
    session.run("codecov", *session.posargs)


def install_with_constraints(session: Session, *args, **kwargs) -> None:
    with tempfile.TemporaryDirectory() as tmpdirname:
        session.run(
            "poetry",
            "export",
            "--dev",
            "--format=requirements.txt",
            "--without-hashes",
            f"--output={os.path.join(tmpdirname, 'poetry-export.out')}",
            external=True,
        )
        session.install(
            f"--constraint={os.path.join(tmpdirname, 'poetry-export.out')}",
            *args,
            **kwargs,
        )
