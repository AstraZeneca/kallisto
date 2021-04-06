# tests/test_alp.py

from tests.store import acetylene

# global epsilon
epsilon = 1e-06


def test_alp():
    charge = 0
    mol = acetylene()
    alp = mol.get_alp(charge)
    assert (alp[0] - 6.56554674) < epsilon
    assert (alp[1] - 1.75193793) < epsilon


def test_alp_cation():
    charge = 1
    mol = acetylene()
    alp = mol.get_alp(charge)
    assert (alp[0] - 4.81427623) < epsilon
    assert (alp[1] - 1.07449632) < epsilon


def test_alp_anion():
    charge = -1
    mol = acetylene()
    alp = mol.get_alp(charge)
    assert (alp[0] - 9.42323948) < epsilon
    assert (alp[1] - 3.25063226) < epsilon
