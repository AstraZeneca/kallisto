# tests/test_eeq.py

from kallisto.units import Bohr
from tests.store import acetylene


# global epsilon
epsilon = 1e-06


def test_vdw_rahm():
    charge = 0
    vdwtype = "rahm"
    scale = 1
    mol = acetylene()
    vdw = mol.get_vdw(charge, vdwtype, scale)
    assert (vdw[0] - 3.29019696) < epsilon
    assert (vdw[1] - 2.5041682) < epsilon


def test_vdw_rahm_angstrom():
    charge = 0
    vdwtype = "rahm"
    scale = Bohr
    mol = acetylene()
    vdw = mol.get_vdw(charge, vdwtype, scale)
    assert (vdw[0] - 1.74109725) < epsilon
    assert (vdw[1] - 1.32514874) < epsilon


def test_vdw_rahm_cation():
    charge = 1
    vdwtype = "rahm"
    scale = 1
    mol = acetylene()
    vdw = mol.get_vdw(charge, vdwtype, scale)
    assert (vdw[0] - 3.14755509) < epsilon
    assert (vdw[1] - 2.33524772) < epsilon


def test_vdw_rahm_anion():
    charge = -1
    vdwtype = "rahm"
    scale = 1
    mol = acetylene()
    vdw = mol.get_vdw(charge, vdwtype, scale)
    assert (vdw[0] - 3.46449846) < epsilon
    assert (vdw[1] - 2.73535295) < epsilon


def test_vdw_truhlar():
    charge = 0
    vdwtype = "truhlar"
    scale = 1
    mol = acetylene()
    vdw = mol.get_vdw(charge, vdwtype, scale)
    assert (vdw[0] - 2.95785383) < epsilon
    assert (vdw[1] - 1.78869157) < epsilon


def test_vdw_truhlar_angstrom():
    charge = 0
    vdwtype = "truhlar"
    scale = Bohr
    mol = acetylene()
    vdw = mol.get_vdw(charge, vdwtype, scale)
    assert (vdw[0] - 1.56522884) < epsilon
    assert (vdw[1] - 0.94653482) < epsilon


def test_vdw_truhlar_cation():
    charge = 1
    vdwtype = "truhlar"
    scale = 1
    mol = acetylene()
    vdw = mol.get_vdw(charge, vdwtype, scale)
    assert (vdw[0] - 2.82962024) < epsilon
    assert (vdw[1] - 1.66803409) < epsilon


def test_vdw_truhlar_anion():
    charge = -1
    vdwtype = "truhlar"
    scale = 1
    mol = acetylene()
    vdw = mol.get_vdw(charge, vdwtype, scale)
    assert (vdw[0] - 3.11454912) < epsilon
    assert (vdw[1] - 1.95382353) < epsilon
