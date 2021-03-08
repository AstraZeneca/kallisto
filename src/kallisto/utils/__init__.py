# src/kallisto/utils/__init__.py

from math import gcd

import click

__all__ = ["basestring", "gcd"]

basestring = str


def silentPrinter(silent: bool, message: str, out: click.File):
    """Helper function to check for silent mode."""
    if not silent:
        click.echo(message, file=out)  # type: ignore


def errorbye(message: str, out: click.File):
    """Exit application due to error."""
    click.echo("", file=out)  # type: ignore
    click.echo(message)  # type: ignore
    click.echo(
        """-> kallisto was terminated due to an error.
        Please check the output.""",
        file=out,  # type: ignore
    )
    exit()


def goodbye(out: click.File):
    """Helper funtion to say goodbye."""
    click.echo("", file=out)  # type: ignore
    click.echo("All done.", file=out)  # type: ignore
