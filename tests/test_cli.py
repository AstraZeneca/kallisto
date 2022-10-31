# tests/test_cli.py

import os

import click.testing
import pytest

from kallisto.console import cli

# define global lineseperator
s = os.linesep


@pytest.fixture
def runner():
    return click.testing.CliRunner()


def test_main_succeeds(runner):
    result = runner.invoke(cli)
    got = result.exit_code
    want = 0
    assert got == want


# check verbose mode
def test_cli_print_help(runner):
    result = runner.invoke(cli)
    assert result.exit_code == 0
    assert "kallisto" in result.output


# test cli part for coordination number
def test_cli_cns_exp_silent(runner, pyridine_xyz):
    result = runner.invoke(cli, ["--silent", "cns", "--cntype", "exp", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_cns_exp(runner, pyridine_xyz):
    result = runner.invoke(cli, ["cns", "--cntype", "exp", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_cns_cov(runner, pyridine_xyz):
    result = runner.invoke(cli, ["cns", "--cntype", "cov", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_cns_erf(runner, pyridine_xyz):
    result = runner.invoke(cli, ["cns", "--cntype", "erf", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_cns_invalid(runner, pyridine_xyz):
    result = runner.invoke(cli, ["cns", "--cntype", "invalid", pyridine_xyz])
    assert result.exit_code == 1


# test cli part for bonds
def test_cli_bonds_silent(runner, pyridine_xyz):
    result = runner.invoke(cli, ["--silent", "bonds", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_bonds(runner, pyridine_xyz):
    result = runner.invoke(cli, ["bonds", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_bonds_partner(runner, pyridine_xyz):
    result = runner.invoke(cli, ["bonds", "--partner", "0", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_bonds_constrain(runner, pyridine_xyz):
    constrain = "constrain.inp"
    gotFile = os.path.isfile(constrain)
    assert gotFile is False
    # CLI runs without problems
    result = runner.invoke(cli, ["bonds", "--constrain", pyridine_xyz])
    assert result.exit_code == 0
    # check content of constrain.inp
    with open(constrain) as tmpfile:
        lines = tmpfile.readlines()
    assert "$constrain" in lines[0]
    assert "distance: 1, 2, auto" in lines[1]
    assert "distance: 1, 6, auto" in lines[2]
    assert "distance: 1, 7, auto" in lines[3]
    assert "distance: 2, 1, auto" in lines[4]
    assert "distance: 2, 3, auto" in lines[5]
    assert "distance: 2, 8, auto" in lines[6]
    # remove constrain.inp
    gotFile = os.path.isfile(constrain)
    assert gotFile is True
    if gotFile:
        os.remove(constrain)


# test cli part for sort
def test_cli_sort(runner, pyridine_xyz):
    result = runner.invoke(cli, ["sort", pyridine_xyz])
    assert result.exit_code == 0
    assert "11" in result.output
    assert "kallisto" in result.output


def test_cli_sort_with_start(runner, pyridine_xyz):
    result = runner.invoke(cli, ["sort", "--start", "5", pyridine_xyz])
    assert result.exit_code == 0
    assert "11" in result.output
    assert "kallisto\nN      0.6816    1.1960    0.0000" in result.output


# test cli part for eeq
def test_cli_eeq_silent(runner, pyridine_xyz):
    result = runner.invoke(cli, ["--silent", "eeq", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_eeq(runner, pyridine_xyz):
    result = runner.invoke(cli, ["eeq", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_eeq_cation(runner, pyridine_xyz):
    result = runner.invoke(cli, ["eeq", "--chrg", "1", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_eeq_anion(runner, pyridine_xyz):
    result = runner.invoke(cli, ["eeq", "--chrg", "-1", pyridine_xyz])
    assert result.exit_code == 0


# test cli part for alp
def test_cli_alp_silent(runner, pyridine_xyz):
    result = runner.invoke(cli, ["--silent", "alp", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_alp(runner, pyridine_xyz):
    result = runner.invoke(cli, ["alp", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_alp_cation(runner, pyridine_xyz):
    result = runner.invoke(cli, ["alp", "--chrg", "1", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_alp_anion(runner, pyridine_xyz):
    result = runner.invoke(cli, ["alp", "--chrg", "-1", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_alp_molecular(runner, pyridine_xyz):
    result = runner.invoke(cli, ["alp", "--molecular", pyridine_xyz])
    assert result.exit_code == 0


# test cli part for vdw
def test_cli_vdw_silent(runner, pyridine_xyz):
    result = runner.invoke(cli, ["--silent", "vdw", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_vdw(runner, pyridine_xyz):
    result = runner.invoke(cli, ["vdw", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_vdw_angstrom(runner, pyridine_xyz):
    result = runner.invoke(cli, ["vdw", "--angstrom", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_vdw_angstrom_truhlar(runner, pyridine_xyz):
    result = runner.invoke(
        cli,
        ["vdw", "--angstrom", "--vdwtype", "truhlar", pyridine_xyz],
    )
    assert result.exit_code == 0


def test_cli_vdw_angstrom_rahm(runner, pyridine_xyz):
    result = runner.invoke(
        cli, ["vdw", "--angstrom", "--vdwtype", "rahm", pyridine_xyz]
    )
    assert result.exit_code == 0


def test_cli_vdw_angstrom_cation(runner, pyridine_xyz):
    result = runner.invoke(cli, ["vdw", "--angstrom", "--chrg", "1", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_vdw_angstrom_anion(runner, pyridine_xyz):
    result = runner.invoke(cli, ["vdw", "--angstrom", "--chrg", "-1", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_vdw_invalid(runner, pyridine_xyz):
    result = runner.invoke(cli, ["vdw", "--vdwtype", "invalid", pyridine_xyz])
    assert result.exit_code == 1


# test cli part for rms
def test_cli_rms_silent(runner, pyridine_xyz):
    result = runner.invoke(cli, ["--silent", "rms", pyridine_xyz, pyridine_xyz])
    assert result.exit_code == 0


def test_cli_rms_nats_not_equal(runner, pyridine_xyz, ch_radical_xyz):
    result = runner.invoke(cli, ["rms", pyridine_xyz, ch_radical_xyz])
    assert result.exit_code == 1


def test_cli_rms(runner, ch_radical_xyz):
    result = runner.invoke(cli, ["rms", ch_radical_xyz, ch_radical_xyz])
    assert result.exit_code == 0


# test cli part for lig
def test_cli_lig_silent(runner, ch_radical_xyz):
    result = runner.invoke(cli, ["--silent", "lig", "--center", "1", ch_radical_xyz])
    assert result.exit_code == 0


def test_cli_lig_0(runner, ch_radical_xyz):
    result = runner.invoke(cli, ["lig", "--center", "0", ch_radical_xyz])
    assert result.exit_code == 0


# test cli part for exs
def test_cli_exs_silent(runner, pyridine_xyz):
    result = runner.invoke(
        cli,
        [
            "--silent",
            "exs",
            "--center",
            "0",
            "--subnr",
            "2",
            pyridine_xyz,
            pyridine_xyz,
        ],
    )
    assert result.exit_code == 0
    newstructure = "newstructure.xyz"
    gotFile = os.path.isfile(newstructure)
    if gotFile:
        os.remove(newstructure)
    constrain = "constrain.inp"
    gotFile = os.path.isfile(constrain)
    if gotFile:
        os.remove(constrain)


def test_cli_exs(runner, pyridine_xyz):
    result = runner.invoke(
        cli,
        ["exs", "--center", "0", "--subnr", "2", pyridine_xyz, pyridine_xyz],
    )
    assert result.exit_code == 0
    newstructure = "newstructure.xyz"
    gotFile = os.path.isfile(newstructure)
    assert gotFile is True
    if gotFile:
        os.remove(newstructure)
    constrain = "constrain.inp"
    gotFile = os.path.isfile(constrain)
    assert gotFile is True
    if gotFile:
        os.remove(constrain)


def test_cli_exs_with_rotation(runner, pyridine_xyz):
    result = runner.invoke(
        cli,
        [
            "exs",
            "--center",
            "0",
            "--subnr",
            "2",
            "--rotate",
            "180",
            pyridine_xyz,
            pyridine_xyz,
        ],
    )
    assert result.exit_code == 0
    newstructure = "newstructure.xyz"
    gotFile = os.path.isfile(newstructure)
    assert gotFile is True
    if gotFile:
        os.remove(newstructure)
    constrain = "constrain.inp"
    gotFile = os.path.isfile(constrain)
    assert gotFile is True
    if gotFile:
        os.remove(constrain)


# test cli part for stm
def test_cli_stm_silent(runner, iridiumcat_xyz):
    result = runner.invoke(
        cli,
        ["--silent", "stm", "--origin", "18", "--partner", "23", iridiumcat_xyz],
    )
    assert result.exit_code == 0


def test_cli_stm(runner, iridiumcat_xyz):
    result = runner.invoke(
        cli,
        ["stm", "--origin", "18", "--partner", "23", iridiumcat_xyz],
    )
    assert result.exit_code == 0


# test cli part for prox
def test_cli_prox_silent(runner, pyridine_xyz):
    result = runner.invoke(cli, ["--silent", "prox", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_prox(runner, pyridine_xyz):
    result = runner.invoke(cli, ["prox", pyridine_xyz])
    assert result.exit_code == 0


def test_cli_prox_invalid_size(runner, pyridine_xyz):
    result = runner.invoke(cli, ["prox", "--size", "3", "2", pyridine_xyz])
    assert result.exit_code == 1
