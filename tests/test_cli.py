# tests/test_cli.py

import os
import tempfile

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
def test_cli_cns_exp_silent(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(
            cli, ["--silent", "cns", "--inp", f.name, "--cntype", "exp"]
        )
        assert result.exit_code == 0


def test_cli_cns_exp(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["cns", "--inp", f.name, "--cntype", "exp"])
        assert result.exit_code == 0
        assert "0.98678927" in result.output


def test_cli_cns_cov(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["cns", "--inp", f.name, "--cntype", "cov"])
        assert result.exit_code == 0
        assert "0.9189476" in result.output


def test_cli_cns_erf(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["cns", "--inp", f.name, "--cntype", "erf"])
        assert result.exit_code == 0
        assert "0.9878465" in result.output


# test cli part for bonds
def test_cli_bonds_silent(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["--silent", "bonds", "--inp", f.name])
        assert result.exit_code == 0


def test_cli_bonds(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["bonds", "--inp", f.name])
        assert result.exit_code == 0


def test_cli_bonds_partner(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["bonds", "--inp", f.name, "--partner", "0"])
        assert result.exit_code == 0


def test_cli_bonds_constrain(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        constrain = "constrain.inp"
        gotFile = os.path.isfile(constrain)
        assert gotFile is False
        result = runner.invoke(cli, ["bonds", "--inp", f.name, "--constrain"])
        assert result.exit_code == 0
        gotFile = os.path.isfile(constrain)
        assert gotFile is True
        if gotFile:
            os.remove(constrain)


# test cli part for sort
def test_cli_sort(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz") as f:
        f.write("11" + s)
        f.write("Pyridine" + s)
        f.write("C 1.3603 0.0256 0.0000" + s)
        f.write("C 0.6971 -1.2020 0.0000" + s)
        f.write("C -0.6944 -1.2184 0.0000" + s)
        f.write("C -1.3895 -0.0129 0.0000" + s)
        f.write("C -0.6712 1.1834 0.0000" + s)
        f.write("N 0.6816 1.1960 0.0000" + s)
        f.write("H 2.4530 0.1083 0.0000" + s)
        f.write("H 1.2665 -2.1365 0.0000" + s)
        f.write("H -1.2365 -2.1696 0.0000" + s)
        f.write("H -2.4837 0.0011 0.0000" + s)
        f.write("H -1.1569 2.1657 0.0000" + s)
        f.flush()
        result = runner.invoke(cli, ["sort", "--inp", f.name])
        assert result.exit_code == 0
        assert "11" in result.output
        assert "kallisto" in result.output


def test_cli_sort_with_start(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz") as f:
        f.write("11" + s)
        f.write("Pyridine" + s)
        f.write("C 1.3603 0.0256 0.0000" + s)
        f.write("C 0.6971 -1.2020 0.0000" + s)
        f.write("C -0.6944 -1.2184 0.0000" + s)
        f.write("C -1.3895 -0.0129 0.0000" + s)
        f.write("C -0.6712 1.1834 0.0000" + s)
        f.write("N 0.6816 1.1960 0.0000" + s)
        f.write("H 2.4530 0.1083 0.0000" + s)
        f.write("H 1.2665 -2.1365 0.0000" + s)
        f.write("H -1.2365 -2.1696 0.0000" + s)
        f.write("H -2.4837 0.0011 0.0000" + s)
        f.write("H -1.1569 2.1657 0.0000" + s)
        f.flush()
        result = runner.invoke(cli, ["sort", "--inp", f.name, "--start", "5"])
        assert result.exit_code == 0
        assert "11" in result.output
        assert "kallisto\nN      0.6816    1.1960    0.0000" in result.output


# test cli part for eeq
def test_cli_eeq_silent(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["--silent", "eeq", "--inp", f.name])
        assert result.exit_code == 0


def test_cli_eeq(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["eeq", "--inp", f.name])
        assert result.exit_code == 0
        assert "-0.17166856" and "0.17166856" in result.output


def test_cli_eeq_cation(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["eeq", "--inp", f.name, "--chrg", "1"])
        assert result.exit_code == 0
        assert "0.59769359" and "0.40230641" in result.output


def test_cli_eeq_anion(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["eeq", "--inp", f.name, "--chrg", "-1"])
        assert result.exit_code == 0
        assert "-0.94103071" and "-0.058969287" in result.output


# test cli part for alp
def test_cli_alp_silent(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["--silent", "alp", "--inp", f.name])
        assert result.exit_code == 0


def test_cli_alp(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["alp", "--inp", f.name])
        assert result.exit_code == 0
        assert "6.5655467" and "1.7519379" in result.output


def test_cli_alp_cation(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["alp", "--inp", f.name, "--chrg", "1"])
        assert result.exit_code == 0
        assert "4.8142762" and "1.0744963" in result.output


def test_cli_alp_anion(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["alp", "--inp", f.name, "--chrg", "-1"])
        assert result.exit_code == 0
        assert "9.4232394" and "3.2506322" in result.output


def test_cli_alp_molecular(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["alp", "--inp", f.name, "--molecular"])
        assert result.exit_code == 0
        assert "8.317" in result.output


# test cli part for vdw
def test_cli_vdw_silent(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["--silent", "vdw", "--inp", f.name])
        assert result.exit_code == 0


def test_cli_vdw(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["vdw", "--inp", f.name])
        assert result.exit_code == 0
        assert "3.29019696" and "2.5041682" in result.output


def test_cli_vdw_angstrom(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["vdw", "--inp", f.name, "--angstrom"])
        assert result.exit_code == 0
        assert "1.74109" and "1.32514" in result.output


def test_cli_vdw_angstrom_truhlar(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(
            cli, ["vdw", "--inp", f.name, "--angstrom", "--vdwtype", "truhlar"],
        )
        assert result.exit_code == 0
        assert "1.565228" and "0.946534" in result.output


def test_cli_vdw_angstrom_rahm(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(
            cli, ["vdw", "--inp", f.name, "--angstrom", "--vdwtype", "rahm"],
        )
        assert result.exit_code == 0
        assert "1.741096" and "1.325148" in result.output


def test_cli_vdw_angstrom_cation(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(
            cli, ["vdw", "--inp", f.name, "--angstrom", "--chrg", "1"]
        )
        assert result.exit_code == 0
        assert "1.665613" and "1.235759" in result.output


def test_cli_vdw_angstrom_anion(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(
            cli, ["vdw", "--inp", f.name, "--angstrom", "--chrg", "-1"]
        )
        assert result.exit_code == 0
        assert "1.833333" and "1.447486" in result.output


# test cli part for rms
def test_cli_rms_silent(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f1:
        f1.write("$coord" + s)
        f1.write("0 0 -1.06176 c" + s)
        f1.write("0 0  1.06176 h" + s)
        f1.write("$end")
        f1.flush()
        f2 = tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8")
        f2.write("$coord" + s)
        f2.write("1 1 -0.06176 c" + s)
        f2.write("1 1  2.06176 h" + s)
        f2.write("$end")
        f2.flush()
        result = runner.invoke(cli, ["--silent", "rms", "--compare", f1.name, f2.name])
        assert result.exit_code == 0
        f2.close()


def test_cli_rms_nats_not_equal(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f1:
        f1.write("$coord" + s)
        f1.write("0 0 -1.06176 c" + s)
        f1.write("0 0  1.06176 h" + s)
        f1.write("$end")
        f1.flush()
        f2 = tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8")
        f2.write("$coord" + s)
        f2.write("1 1 -0.06176 c" + s)
        f2.write("$end")
        f2.flush()
        result = runner.invoke(cli, ["rms", "--compare", f1.name, f2.name])
        assert result.exit_code == 0
        assert "Error" in result.output
        f2.close()


def test_cli_rms(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f1:
        f1.write("$coord" + s)
        f1.write("0 0 -1.06176 c" + s)
        f1.write("0 0  1.06176 h" + s)
        f1.write("$end")
        f1.flush()
        f2 = tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8")
        f2.write("$coord" + s)
        f2.write("1 1 -0.06176 c" + s)
        f2.write("1 1  2.06176 h" + s)
        f2.write("$end")
        f2.flush()
        result = runner.invoke(cli, ["rms", "--compare", f1.name, f2.name])
        assert result.exit_code == 0
        assert "[1. 0. 0.]" and "[0. 1. 0.]" and "[0. 0. 1.]" in result.output
        f2.close()


# test cli part for lig
def test_cli_lig_silent(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0 0 -1.06176 c" + s)
        f.write("0 0  1.06176 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(
            cli, ["--silent", "lig", "--inp", f.name, "--center", "1"]
        )
        assert result.exit_code == 0


def test_cli_lig_0(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0 0 -1.06176 c" + s)
        f.write("0 0  1.06176 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["lig", "--inp", f.name, "--center", "0"])
        assert result.exit_code == 0
        assert "Substructure" and "[1]" in result.output


def test_cli_lig_1(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0 0 -1.06176 c" + s)
        f.write("0 0  1.06176 h" + s)
        f.write("$end")
        f.flush()
        result = runner.invoke(cli, ["lig", "--inp", f.name, "--center", "1"])
        assert result.exit_code == 0
        assert "Substructure" and "[0]" in result.output


# test cli part for exs
def test_cli_exs_silent(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz") as f1:
        f1.write("96" + s)
        f1.write(s)
        f1.write("N    -1.3672999   -1.4398999    0.1359000" + s)
        f1.write("C    -2.4911998   -0.6808999    0.1396000" + s)
        f1.write("C    -3.6534996   -1.1211999   -0.5090000" + s)
        f1.write("C    -3.6468996   -2.3434998   -1.1725999" + s)
        f1.write("C    -2.4848998   -3.1187997   -1.1555999" + s)
        f1.write("C    -1.3670999   -2.6316997   -0.4883000" + s)
        f1.write("H    -0.4373000   -3.1872997   -0.4306000" + s)
        f1.write("H    -2.4432998   -4.0866996   -1.6463998" + s)
        f1.write("H    -4.5575996   -0.5223999   -0.4887000" + s)
        f1.write("C    -2.4206998    0.5908999    0.8954999" + s)
        f1.write("N    -1.2879999    0.7903999    1.6181998" + s)
        f1.write("C    -1.1378999    1.9348998    2.3084998" + s)
        f1.write("C    -2.1077998    2.9319997    2.3219998" + s)
        f1.write("C    -3.2770997    2.7402997    1.5819998" + s)
        f1.write("C    -3.4330997    1.5600998    0.8608999" + s)
        f1.write("H    -4.3267996    1.4057999    0.2659000" + s)
        f1.write("H    -1.9411998    3.8419996    2.8913997" + s)
        f1.write("H    -0.1872000    2.0459998    2.8181997" + s)
        f1.write("Ir    0.4009000   -0.6061999    1.1172999" + s)
        f1.write("C    -1.2690999   -3.8143996    3.7856996" + s)
        f1.write("C    -0.1664000   -4.5494996    4.2269996" + s)
        f1.write("C     1.1218999   -4.0950996    3.9273996" + s)
        f1.write("C     1.2993999   -2.9384997    3.1675997" + s)
        f1.write("C     0.2001000   -2.2075998    2.6786997" + s)
        f1.write("C    -1.0849999   -2.6466997    3.0382997" + s)
        f1.write("H    -1.9573998   -2.0870998    2.7090997" + s)
        f1.write("H     0.8509999   -0.7173999    2.6636997" + s)
        f1.write("H     2.3007998   -2.5989997    2.9226997" + s)
        f1.write("H    -0.3087000   -5.4547995    4.8119995" + s)
        f1.write("B     0.6392999    0.6220999   -0.5923999" + s)
        f1.write("O    -0.0586000    0.3754000   -1.7751998" + s)
        f1.write("C     0.0637000    1.5387999   -2.6275997" + s)
        f1.write("C     0.0955000    1.0794999   -4.0821996" + s)
        f1.write("H     0.8716999    0.3276000   -4.2397996" + s)
        f1.write("H     0.2802000    1.9248998   -4.7547995" + s)
        f1.write("H    -0.8681999    0.6330999   -4.3491996" + s)
        f1.write("C    -1.1760999    2.4077998   -2.3666998" + s)
        f1.write("H    -1.2042999    3.2867997   -3.0193997" + s)
        f1.write("H    -2.0717998    1.8058998   -2.5499998" + s)
        f1.write("H    -1.2019999    2.7410997   -1.3247999" + s)
        f1.write("C     1.3891999    2.1923998   -2.1029998" + s)
        f1.write("O     1.3915999    1.7859998   -0.7128999" + s)
        f1.write("C     2.6481997    1.5975998   -2.7492997" + s)
        f1.write("H     2.6573997    0.5124000   -2.6283997" + s)
        f1.write("H     2.7186997    1.8556998   -3.8108996" + s)
        f1.write("H     3.5309997    1.9918998   -2.2375998" + s)
        f1.write("C     1.4299999    3.7172996   -2.1670998" + s)
        f1.write("H     0.6241999    4.1645996   -1.5814998" + s)
        f1.write("H     2.3812998    4.0763996   -1.7610998" + s)
        f1.write("H     1.3476999    4.0651996   -3.2032997" + s)
        f1.write("B     2.0756998    0.4378000    1.7666998" + s)
        f1.write("O     3.3654997   -0.0810000    1.8671998" + s)
        f1.write("C     4.2709996    1.0315999    2.0579998" + s)
        f1.write("C     5.4819995    0.5533999    2.8527997" + s)
        f1.write("H     5.1838995    0.0640000    3.7825996" + s)
        f1.write("H     6.0490994   -0.1689000    2.2567998" + s)
        f1.write("H     6.1442994    1.3928999    3.0939997" + s)
        f1.write("C     4.6909995    1.5017999    0.6583999" + s)
        f1.write("H     5.4398995    2.2997998    0.7063999" + s)
        f1.write("H     5.1090995    0.6483999    0.1171000" + s)
        f1.write("H     3.8181996    1.8505998    0.1015000" + s)
        f1.write("C     3.3502997    2.0656998    2.7942997" + s)
        f1.write("O     2.0458998    1.7442998    2.2507998" + s)
        f1.write("C     3.6590996    3.5300997    2.4959998" + s)
        f1.write("H     3.5358997    3.7464996    1.4332999" + s)
        f1.write("H     2.9753997    4.1760996    3.0565997" + s)
        f1.write("H     4.6847995    3.7776996    2.7930997" + s)
        f1.write("C     3.2796997    1.8273998    4.3083996" + s)
        f1.write("H     3.0708997    0.7747999    4.5225996" + s)
        f1.write("H     4.2123996    2.1093998    4.8080995" + s)
        f1.write("H     2.4671998    2.4302998    4.7258995" + s)
        f1.write("B     1.7917998   -1.7489998    0.1412000" + s)
        f1.write("O     1.8500998   -3.1467997    0.2110000" + s)
        f1.write("C     3.0569997   -3.5802997   -0.4612000" + s)
        f1.write("C     4.1632996   -3.6178996    0.6029999" + s)
        f1.write("H     4.3272996   -2.6173997    1.0149999" + s)
        f1.write("H     3.8420996   -4.2739996    1.4174999" + s)
        f1.write("H     5.1055995   -3.9992996    0.1957000" + s)
        f1.write("C     2.8261997   -4.9693995   -1.0466999" + s)
        f1.write("H     2.6792997   -5.6906994   -0.2364000" + s)
        f1.write("H     3.6907996   -5.2904995   -1.6389998" + s)
        f1.write("H     1.9392998   -4.9911995   -1.6839998" + s)
        f1.write("C     3.2640997   -2.4330998   -1.5048999" + s)
        f1.write("O     2.7238997   -1.2952999   -0.7998999" + s)
        f1.write("C     4.7166995   -2.1413998   -1.8718998" + s)
        f1.write("H     5.1881995   -3.0188997   -2.3293998" + s)
        f1.write("H     4.7565995   -1.3170999   -2.5912997" + s)
        f1.write("H     5.2941995   -1.8488998   -0.9925999" + s)
        f1.write("C     2.4197998   -2.6279997   -2.7728997" + s)
        f1.write("H     1.3752999   -2.8206997   -2.5121998" + s)
        f1.write("H     2.4501998   -1.7101998   -3.3667997" + s)
        f1.write("H     2.7924997   -3.4536997   -3.3878997" + s)
        f1.write("H    -2.2764998   -4.1481996    4.0262996" + s)
        f1.write("H     1.9905998   -4.6454995    4.2840996" + s)
        f1.write("H    -4.5414996   -2.6926997   -1.6821998" + s)
        f1.write("H    -4.0522996    3.5020997    1.5576998" + s)
        f1.flush()
        f2 = tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz")
        f2.write("10" + s)
        f2.write(s)
        f2.write("  C      1.3603      0.0256      0.0000" + s)
        f2.write("  C      0.6971     -1.2020      0.0000" + s)
        f2.write("  C     -0.6944     -1.2184      0.0000" + s)
        f2.write("  C     -1.3895     -0.0129      0.0000" + s)
        f2.write("  C     -0.6712      1.1834      0.0000" + s)
        f2.write("  N      0.6816      1.1960      0.0000" + s)
        f2.write("  H      1.2665     -2.1365      0.0000" + s)
        f2.write("  H     -1.2365     -2.1696      0.0000" + s)
        f2.write("  H     -2.4837      0.0011      0.0000" + s)
        f2.write("  H     -1.1569      2.1657      0.0000" + s)
        f2.flush()
        result = runner.invoke(
            cli,
            [
                "--silent",
                "exs",
                "--inp",
                f1.name,
                f2.name,
                "--center",
                "2",
                "--subnr",
                "2",
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


def test_cli_exs(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz") as f1:
        f1.write("96" + s)
        f1.write(s)
        f1.write("N    -1.3672999   -1.4398999    0.1359000" + s)
        f1.write("C    -2.4911998   -0.6808999    0.1396000" + s)
        f1.write("C    -3.6534996   -1.1211999   -0.5090000" + s)
        f1.write("C    -3.6468996   -2.3434998   -1.1725999" + s)
        f1.write("C    -2.4848998   -3.1187997   -1.1555999" + s)
        f1.write("C    -1.3670999   -2.6316997   -0.4883000" + s)
        f1.write("H    -0.4373000   -3.1872997   -0.4306000" + s)
        f1.write("H    -2.4432998   -4.0866996   -1.6463998" + s)
        f1.write("H    -4.5575996   -0.5223999   -0.4887000" + s)
        f1.write("C    -2.4206998    0.5908999    0.8954999" + s)
        f1.write("N    -1.2879999    0.7903999    1.6181998" + s)
        f1.write("C    -1.1378999    1.9348998    2.3084998" + s)
        f1.write("C    -2.1077998    2.9319997    2.3219998" + s)
        f1.write("C    -3.2770997    2.7402997    1.5819998" + s)
        f1.write("C    -3.4330997    1.5600998    0.8608999" + s)
        f1.write("H    -4.3267996    1.4057999    0.2659000" + s)
        f1.write("H    -1.9411998    3.8419996    2.8913997" + s)
        f1.write("H    -0.1872000    2.0459998    2.8181997" + s)
        f1.write("Ir    0.4009000   -0.6061999    1.1172999" + s)
        f1.write("C    -1.2690999   -3.8143996    3.7856996" + s)
        f1.write("C    -0.1664000   -4.5494996    4.2269996" + s)
        f1.write("C     1.1218999   -4.0950996    3.9273996" + s)
        f1.write("C     1.2993999   -2.9384997    3.1675997" + s)
        f1.write("C     0.2001000   -2.2075998    2.6786997" + s)
        f1.write("C    -1.0849999   -2.6466997    3.0382997" + s)
        f1.write("H    -1.9573998   -2.0870998    2.7090997" + s)
        f1.write("H     0.8509999   -0.7173999    2.6636997" + s)
        f1.write("H     2.3007998   -2.5989997    2.9226997" + s)
        f1.write("H    -0.3087000   -5.4547995    4.8119995" + s)
        f1.write("B     0.6392999    0.6220999   -0.5923999" + s)
        f1.write("O    -0.0586000    0.3754000   -1.7751998" + s)
        f1.write("C     0.0637000    1.5387999   -2.6275997" + s)
        f1.write("C     0.0955000    1.0794999   -4.0821996" + s)
        f1.write("H     0.8716999    0.3276000   -4.2397996" + s)
        f1.write("H     0.2802000    1.9248998   -4.7547995" + s)
        f1.write("H    -0.8681999    0.6330999   -4.3491996" + s)
        f1.write("C    -1.1760999    2.4077998   -2.3666998" + s)
        f1.write("H    -1.2042999    3.2867997   -3.0193997" + s)
        f1.write("H    -2.0717998    1.8058998   -2.5499998" + s)
        f1.write("H    -1.2019999    2.7410997   -1.3247999" + s)
        f1.write("C     1.3891999    2.1923998   -2.1029998" + s)
        f1.write("O     1.3915999    1.7859998   -0.7128999" + s)
        f1.write("C     2.6481997    1.5975998   -2.7492997" + s)
        f1.write("H     2.6573997    0.5124000   -2.6283997" + s)
        f1.write("H     2.7186997    1.8556998   -3.8108996" + s)
        f1.write("H     3.5309997    1.9918998   -2.2375998" + s)
        f1.write("C     1.4299999    3.7172996   -2.1670998" + s)
        f1.write("H     0.6241999    4.1645996   -1.5814998" + s)
        f1.write("H     2.3812998    4.0763996   -1.7610998" + s)
        f1.write("H     1.3476999    4.0651996   -3.2032997" + s)
        f1.write("B     2.0756998    0.4378000    1.7666998" + s)
        f1.write("O     3.3654997   -0.0810000    1.8671998" + s)
        f1.write("C     4.2709996    1.0315999    2.0579998" + s)
        f1.write("C     5.4819995    0.5533999    2.8527997" + s)
        f1.write("H     5.1838995    0.0640000    3.7825996" + s)
        f1.write("H     6.0490994   -0.1689000    2.2567998" + s)
        f1.write("H     6.1442994    1.3928999    3.0939997" + s)
        f1.write("C     4.6909995    1.5017999    0.6583999" + s)
        f1.write("H     5.4398995    2.2997998    0.7063999" + s)
        f1.write("H     5.1090995    0.6483999    0.1171000" + s)
        f1.write("H     3.8181996    1.8505998    0.1015000" + s)
        f1.write("C     3.3502997    2.0656998    2.7942997" + s)
        f1.write("O     2.0458998    1.7442998    2.2507998" + s)
        f1.write("C     3.6590996    3.5300997    2.4959998" + s)
        f1.write("H     3.5358997    3.7464996    1.4332999" + s)
        f1.write("H     2.9753997    4.1760996    3.0565997" + s)
        f1.write("H     4.6847995    3.7776996    2.7930997" + s)
        f1.write("C     3.2796997    1.8273998    4.3083996" + s)
        f1.write("H     3.0708997    0.7747999    4.5225996" + s)
        f1.write("H     4.2123996    2.1093998    4.8080995" + s)
        f1.write("H     2.4671998    2.4302998    4.7258995" + s)
        f1.write("B     1.7917998   -1.7489998    0.1412000" + s)
        f1.write("O     1.8500998   -3.1467997    0.2110000" + s)
        f1.write("C     3.0569997   -3.5802997   -0.4612000" + s)
        f1.write("C     4.1632996   -3.6178996    0.6029999" + s)
        f1.write("H     4.3272996   -2.6173997    1.0149999" + s)
        f1.write("H     3.8420996   -4.2739996    1.4174999" + s)
        f1.write("H     5.1055995   -3.9992996    0.1957000" + s)
        f1.write("C     2.8261997   -4.9693995   -1.0466999" + s)
        f1.write("H     2.6792997   -5.6906994   -0.2364000" + s)
        f1.write("H     3.6907996   -5.2904995   -1.6389998" + s)
        f1.write("H     1.9392998   -4.9911995   -1.6839998" + s)
        f1.write("C     3.2640997   -2.4330998   -1.5048999" + s)
        f1.write("O     2.7238997   -1.2952999   -0.7998999" + s)
        f1.write("C     4.7166995   -2.1413998   -1.8718998" + s)
        f1.write("H     5.1881995   -3.0188997   -2.3293998" + s)
        f1.write("H     4.7565995   -1.3170999   -2.5912997" + s)
        f1.write("H     5.2941995   -1.8488998   -0.9925999" + s)
        f1.write("C     2.4197998   -2.6279997   -2.7728997" + s)
        f1.write("H     1.3752999   -2.8206997   -2.5121998" + s)
        f1.write("H     2.4501998   -1.7101998   -3.3667997" + s)
        f1.write("H     2.7924997   -3.4536997   -3.3878997" + s)
        f1.write("H    -2.2764998   -4.1481996    4.0262996" + s)
        f1.write("H     1.9905998   -4.6454995    4.2840996" + s)
        f1.write("H    -4.5414996   -2.6926997   -1.6821998" + s)
        f1.write("H    -4.0522996    3.5020997    1.5576998" + s)
        f1.flush()
        f2 = tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz")
        f2.write("10" + s)
        f2.write(s)
        f2.write("  C      1.3603      0.0256      0.0000" + s)
        f2.write("  C      0.6971     -1.2020      0.0000" + s)
        f2.write("  C     -0.6944     -1.2184      0.0000" + s)
        f2.write("  C     -1.3895     -0.0129      0.0000" + s)
        f2.write("  C     -0.6712      1.1834      0.0000" + s)
        f2.write("  N      0.6816      1.1960      0.0000" + s)
        f2.write("  H      1.2665     -2.1365      0.0000" + s)
        f2.write("  H     -1.2365     -2.1696      0.0000" + s)
        f2.write("  H     -2.4837      0.0011      0.0000" + s)
        f2.write("  H     -1.1569      2.1657      0.0000" + s)
        f2.flush()
        newstructure = "newstructure.xyz"
        gotFile = os.path.isfile(newstructure)
        assert gotFile is False
        result = runner.invoke(
            cli, ["exs", "--inp", f1.name, f2.name, "--center", "18", "--subnr", "2"],
        )
        assert result.exit_code == 0
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
def test_cli_stm_silent(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz") as f:
        f.write("96" + s)
        f.write(s)
        f.write("N    -1.3672999   -1.4398999    0.1359000" + s)
        f.write("C    -2.4911998   -0.6808999    0.1396000" + s)
        f.write("C    -3.6534996   -1.1211999   -0.5090000" + s)
        f.write("C    -3.6468996   -2.3434998   -1.1725999" + s)
        f.write("C    -2.4848998   -3.1187997   -1.1555999" + s)
        f.write("C    -1.3670999   -2.6316997   -0.4883000" + s)
        f.write("H    -0.4373000   -3.1872997   -0.4306000" + s)
        f.write("H    -2.4432998   -4.0866996   -1.6463998" + s)
        f.write("H    -4.5575996   -0.5223999   -0.4887000" + s)
        f.write("C    -2.4206998    0.5908999    0.8954999" + s)
        f.write("N    -1.2879999    0.7903999    1.6181998" + s)
        f.write("C    -1.1378999    1.9348998    2.3084998" + s)
        f.write("C    -2.1077998    2.9319997    2.3219998" + s)
        f.write("C    -3.2770997    2.7402997    1.5819998" + s)
        f.write("C    -3.4330997    1.5600998    0.8608999" + s)
        f.write("H    -4.3267996    1.4057999    0.2659000" + s)
        f.write("H    -1.9411998    3.8419996    2.8913997" + s)
        f.write("H    -0.1872000    2.0459998    2.8181997" + s)
        f.write("Ir    0.4009000   -0.6061999    1.1172999" + s)
        f.write("C    -1.2690999   -3.8143996    3.7856996" + s)
        f.write("C    -0.1664000   -4.5494996    4.2269996" + s)
        f.write("C     1.1218999   -4.0950996    3.9273996" + s)
        f.write("C     1.2993999   -2.9384997    3.1675997" + s)
        f.write("C     0.2001000   -2.2075998    2.6786997" + s)
        f.write("C    -1.0849999   -2.6466997    3.0382997" + s)
        f.write("H    -1.9573998   -2.0870998    2.7090997" + s)
        f.write("H     0.8509999   -0.7173999    2.6636997" + s)
        f.write("H     2.3007998   -2.5989997    2.9226997" + s)
        f.write("H    -0.3087000   -5.4547995    4.8119995" + s)
        f.write("B     0.6392999    0.6220999   -0.5923999" + s)
        f.write("O    -0.0586000    0.3754000   -1.7751998" + s)
        f.write("C     0.0637000    1.5387999   -2.6275997" + s)
        f.write("C     0.0955000    1.0794999   -4.0821996" + s)
        f.write("H     0.8716999    0.3276000   -4.2397996" + s)
        f.write("H     0.2802000    1.9248998   -4.7547995" + s)
        f.write("H    -0.8681999    0.6330999   -4.3491996" + s)
        f.write("C    -1.1760999    2.4077998   -2.3666998" + s)
        f.write("H    -1.2042999    3.2867997   -3.0193997" + s)
        f.write("H    -2.0717998    1.8058998   -2.5499998" + s)
        f.write("H    -1.2019999    2.7410997   -1.3247999" + s)
        f.write("C     1.3891999    2.1923998   -2.1029998" + s)
        f.write("O     1.3915999    1.7859998   -0.7128999" + s)
        f.write("C     2.6481997    1.5975998   -2.7492997" + s)
        f.write("H     2.6573997    0.5124000   -2.6283997" + s)
        f.write("H     2.7186997    1.8556998   -3.8108996" + s)
        f.write("H     3.5309997    1.9918998   -2.2375998" + s)
        f.write("C     1.4299999    3.7172996   -2.1670998" + s)
        f.write("H     0.6241999    4.1645996   -1.5814998" + s)
        f.write("H     2.3812998    4.0763996   -1.7610998" + s)
        f.write("H     1.3476999    4.0651996   -3.2032997" + s)
        f.write("B     2.0756998    0.4378000    1.7666998" + s)
        f.write("O     3.3654997   -0.0810000    1.8671998" + s)
        f.write("C     4.2709996    1.0315999    2.0579998" + s)
        f.write("C     5.4819995    0.5533999    2.8527997" + s)
        f.write("H     5.1838995    0.0640000    3.7825996" + s)
        f.write("H     6.0490994   -0.1689000    2.2567998" + s)
        f.write("H     6.1442994    1.3928999    3.0939997" + s)
        f.write("C     4.6909995    1.5017999    0.6583999" + s)
        f.write("H     5.4398995    2.2997998    0.7063999" + s)
        f.write("H     5.1090995    0.6483999    0.1171000" + s)
        f.write("H     3.8181996    1.8505998    0.1015000" + s)
        f.write("C     3.3502997    2.0656998    2.7942997" + s)
        f.write("O     2.0458998    1.7442998    2.2507998" + s)
        f.write("C     3.6590996    3.5300997    2.4959998" + s)
        f.write("H     3.5358997    3.7464996    1.4332999" + s)
        f.write("H     2.9753997    4.1760996    3.0565997" + s)
        f.write("H     4.6847995    3.7776996    2.7930997" + s)
        f.write("C     3.2796997    1.8273998    4.3083996" + s)
        f.write("H     3.0708997    0.7747999    4.5225996" + s)
        f.write("H     4.2123996    2.1093998    4.8080995" + s)
        f.write("H     2.4671998    2.4302998    4.7258995" + s)
        f.write("B     1.7917998   -1.7489998    0.1412000" + s)
        f.write("O     1.8500998   -3.1467997    0.2110000" + s)
        f.write("C     3.0569997   -3.5802997   -0.4612000" + s)
        f.write("C     4.1632996   -3.6178996    0.6029999" + s)
        f.write("H     4.3272996   -2.6173997    1.0149999" + s)
        f.write("H     3.8420996   -4.2739996    1.4174999" + s)
        f.write("H     5.1055995   -3.9992996    0.1957000" + s)
        f.write("C     2.8261997   -4.9693995   -1.0466999" + s)
        f.write("H     2.6792997   -5.6906994   -0.2364000" + s)
        f.write("H     3.6907996   -5.2904995   -1.6389998" + s)
        f.write("H     1.9392998   -4.9911995   -1.6839998" + s)
        f.write("C     3.2640997   -2.4330998   -1.5048999" + s)
        f.write("O     2.7238997   -1.2952999   -0.7998999" + s)
        f.write("C     4.7166995   -2.1413998   -1.8718998" + s)
        f.write("H     5.1881995   -3.0188997   -2.3293998" + s)
        f.write("H     4.7565995   -1.3170999   -2.5912997" + s)
        f.write("H     5.2941995   -1.8488998   -0.9925999" + s)
        f.write("C     2.4197998   -2.6279997   -2.7728997" + s)
        f.write("H     1.3752999   -2.8206997   -2.5121998" + s)
        f.write("H     2.4501998   -1.7101998   -3.3667997" + s)
        f.write("H     2.7924997   -3.4536997   -3.3878997" + s)
        f.write("H    -2.2764998   -4.1481996    4.0262996" + s)
        f.write("H     1.9905998   -4.6454995    4.2840996" + s)
        f.write("H    -4.5414996   -2.6926997   -1.6821998" + s)
        f.write("H    -4.0522996    3.5020997    1.5576998" + s)
        f.flush()
        result = runner.invoke(
            cli,
            ["--silent", "stm", "--inp", f.name, "--origin", "18", "--partner", "23"],
        )
        assert result.exit_code == 0


def test_cli_stm(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz") as f:
        f.write("96" + s)
        f.write(s)
        f.write("N    -1.3672999   -1.4398999    0.1359000" + s)
        f.write("C    -2.4911998   -0.6808999    0.1396000" + s)
        f.write("C    -3.6534996   -1.1211999   -0.5090000" + s)
        f.write("C    -3.6468996   -2.3434998   -1.1725999" + s)
        f.write("C    -2.4848998   -3.1187997   -1.1555999" + s)
        f.write("C    -1.3670999   -2.6316997   -0.4883000" + s)
        f.write("H    -0.4373000   -3.1872997   -0.4306000" + s)
        f.write("H    -2.4432998   -4.0866996   -1.6463998" + s)
        f.write("H    -4.5575996   -0.5223999   -0.4887000" + s)
        f.write("C    -2.4206998    0.5908999    0.8954999" + s)
        f.write("N    -1.2879999    0.7903999    1.6181998" + s)
        f.write("C    -1.1378999    1.9348998    2.3084998" + s)
        f.write("C    -2.1077998    2.9319997    2.3219998" + s)
        f.write("C    -3.2770997    2.7402997    1.5819998" + s)
        f.write("C    -3.4330997    1.5600998    0.8608999" + s)
        f.write("H    -4.3267996    1.4057999    0.2659000" + s)
        f.write("H    -1.9411998    3.8419996    2.8913997" + s)
        f.write("H    -0.1872000    2.0459998    2.8181997" + s)
        f.write("Ir    0.4009000   -0.6061999    1.1172999" + s)
        f.write("C    -1.2690999   -3.8143996    3.7856996" + s)
        f.write("C    -0.1664000   -4.5494996    4.2269996" + s)
        f.write("C     1.1218999   -4.0950996    3.9273996" + s)
        f.write("C     1.2993999   -2.9384997    3.1675997" + s)
        f.write("C     0.2001000   -2.2075998    2.6786997" + s)
        f.write("C    -1.0849999   -2.6466997    3.0382997" + s)
        f.write("H    -1.9573998   -2.0870998    2.7090997" + s)
        f.write("H     0.8509999   -0.7173999    2.6636997" + s)
        f.write("H     2.3007998   -2.5989997    2.9226997" + s)
        f.write("H    -0.3087000   -5.4547995    4.8119995" + s)
        f.write("B     0.6392999    0.6220999   -0.5923999" + s)
        f.write("O    -0.0586000    0.3754000   -1.7751998" + s)
        f.write("C     0.0637000    1.5387999   -2.6275997" + s)
        f.write("C     0.0955000    1.0794999   -4.0821996" + s)
        f.write("H     0.8716999    0.3276000   -4.2397996" + s)
        f.write("H     0.2802000    1.9248998   -4.7547995" + s)
        f.write("H    -0.8681999    0.6330999   -4.3491996" + s)
        f.write("C    -1.1760999    2.4077998   -2.3666998" + s)
        f.write("H    -1.2042999    3.2867997   -3.0193997" + s)
        f.write("H    -2.0717998    1.8058998   -2.5499998" + s)
        f.write("H    -1.2019999    2.7410997   -1.3247999" + s)
        f.write("C     1.3891999    2.1923998   -2.1029998" + s)
        f.write("O     1.3915999    1.7859998   -0.7128999" + s)
        f.write("C     2.6481997    1.5975998   -2.7492997" + s)
        f.write("H     2.6573997    0.5124000   -2.6283997" + s)
        f.write("H     2.7186997    1.8556998   -3.8108996" + s)
        f.write("H     3.5309997    1.9918998   -2.2375998" + s)
        f.write("C     1.4299999    3.7172996   -2.1670998" + s)
        f.write("H     0.6241999    4.1645996   -1.5814998" + s)
        f.write("H     2.3812998    4.0763996   -1.7610998" + s)
        f.write("H     1.3476999    4.0651996   -3.2032997" + s)
        f.write("B     2.0756998    0.4378000    1.7666998" + s)
        f.write("O     3.3654997   -0.0810000    1.8671998" + s)
        f.write("C     4.2709996    1.0315999    2.0579998" + s)
        f.write("C     5.4819995    0.5533999    2.8527997" + s)
        f.write("H     5.1838995    0.0640000    3.7825996" + s)
        f.write("H     6.0490994   -0.1689000    2.2567998" + s)
        f.write("H     6.1442994    1.3928999    3.0939997" + s)
        f.write("C     4.6909995    1.5017999    0.6583999" + s)
        f.write("H     5.4398995    2.2997998    0.7063999" + s)
        f.write("H     5.1090995    0.6483999    0.1171000" + s)
        f.write("H     3.8181996    1.8505998    0.1015000" + s)
        f.write("C     3.3502997    2.0656998    2.7942997" + s)
        f.write("O     2.0458998    1.7442998    2.2507998" + s)
        f.write("C     3.6590996    3.5300997    2.4959998" + s)
        f.write("H     3.5358997    3.7464996    1.4332999" + s)
        f.write("H     2.9753997    4.1760996    3.0565997" + s)
        f.write("H     4.6847995    3.7776996    2.7930997" + s)
        f.write("C     3.2796997    1.8273998    4.3083996" + s)
        f.write("H     3.0708997    0.7747999    4.5225996" + s)
        f.write("H     4.2123996    2.1093998    4.8080995" + s)
        f.write("H     2.4671998    2.4302998    4.7258995" + s)
        f.write("B     1.7917998   -1.7489998    0.1412000" + s)
        f.write("O     1.8500998   -3.1467997    0.2110000" + s)
        f.write("C     3.0569997   -3.5802997   -0.4612000" + s)
        f.write("C     4.1632996   -3.6178996    0.6029999" + s)
        f.write("H     4.3272996   -2.6173997    1.0149999" + s)
        f.write("H     3.8420996   -4.2739996    1.4174999" + s)
        f.write("H     5.1055995   -3.9992996    0.1957000" + s)
        f.write("C     2.8261997   -4.9693995   -1.0466999" + s)
        f.write("H     2.6792997   -5.6906994   -0.2364000" + s)
        f.write("H     3.6907996   -5.2904995   -1.6389998" + s)
        f.write("H     1.9392998   -4.9911995   -1.6839998" + s)
        f.write("C     3.2640997   -2.4330998   -1.5048999" + s)
        f.write("O     2.7238997   -1.2952999   -0.7998999" + s)
        f.write("C     4.7166995   -2.1413998   -1.8718998" + s)
        f.write("H     5.1881995   -3.0188997   -2.3293998" + s)
        f.write("H     4.7565995   -1.3170999   -2.5912997" + s)
        f.write("H     5.2941995   -1.8488998   -0.9925999" + s)
        f.write("C     2.4197998   -2.6279997   -2.7728997" + s)
        f.write("H     1.3752999   -2.8206997   -2.5121998" + s)
        f.write("H     2.4501998   -1.7101998   -3.3667997" + s)
        f.write("H     2.7924997   -3.4536997   -3.3878997" + s)
        f.write("H    -2.2764998   -4.1481996    4.0262996" + s)
        f.write("H     1.9905998   -4.6454995    4.2840996" + s)
        f.write("H    -4.5414996   -2.6926997   -1.6821998" + s)
        f.write("H    -4.0522996    3.5020997    1.5576998" + s)
        f.flush()
        result = runner.invoke(
            cli, ["stm", "--inp", f.name, "--origin", "18", "--partner", "23"],
        )
        assert result.exit_code == 0
        # L value in Bohr and Angstrom
        assert "14.08" and "7.45" in result.output
        # Bmin value in Bohr and Angstrom
        assert "11.18" and "5.92" in result.output
        # Bmax value in Bohr and Angstrom
        assert "14.62" and "7.74" in result.output


# test cli part for cnsp
def test_cli_cnsp_silent(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz") as f:
        f.write("15" + s)
        f.write("toluene" + s)
        f.write("c 1.2264 0.0427 0.0670" + s)
        f.write("c 1.0031 -1.3293 0.0600" + s)
        f.write("c -0.2945 -1.8256 -0.0060" + s)
        f.write("c -1.3704 -0.9461 -0.0646" + s)
        f.write("c -1.1511 0.4266 -0.0578" + s)
        f.write("c 0.1497 0.9292 0.0066" + s)
        f.write("c 0.3871 2.3956 -0.0022" + s)
        f.write("h 2.2495 0.4310 0.1211" + s)
        f.write("h 1.8510 -2.0202 0.1071" + s)
        f.write("h -0.4688 -2.9062 -0.0109" + s)
        f.write("h -2.3926 -1.3347 -0.1157" + s)
        f.write("h -2.0006 1.1172 -0.1021" + s)
        f.write("h 0.5024 2.7582 -1.0330" + s)
        f.write("h 1.2994 2.6647 0.5466" + s)
        f.write("h -0.4475 2.9470 0.4506" + s)
        f.flush()
        result = runner.invoke(cli, ["--silent", "cnsp", "--inp", f.name])
        assert result.exit_code == 0


def test_cli_cnsp(runner):
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz") as f:
        f.write("15" + s)
        f.write("toluene" + s)
        f.write("c 1.2264 0.0427 0.0670" + s)
        f.write("c 1.0031 -1.3293 0.0600" + s)
        f.write("c -0.2945 -1.8256 -0.0060" + s)
        f.write("c -1.3704 -0.9461 -0.0646" + s)
        f.write("c -1.1511 0.4266 -0.0578" + s)
        f.write("c 0.1497 0.9292 0.0066" + s)
        f.write("c 0.3871 2.3956 -0.0022" + s)
        f.write("h 2.2495 0.4310 0.1211" + s)
        f.write("h 1.8510 -2.0202 0.1071" + s)
        f.write("h -0.4688 -2.9062 -0.0109" + s)
        f.write("h -2.3926 -1.3347 -0.1157" + s)
        f.write("h -2.0006 1.1172 -0.1021" + s)
        f.write("h 0.5024 2.7582 -1.0330" + s)
        f.write("h 1.2994 2.6647 0.5466" + s)
        f.write("h -0.4475 2.9470 0.4506" + s)
        f.flush()
        result = runner.invoke(cli, ["cnsp", "--inp", f.name])
        assert result.exit_code == 0
        assert "4.38" and "3.36" in result.output
