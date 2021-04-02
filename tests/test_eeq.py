# tests/test_eeq.py

from tests.store import acetylene

# global epsilon
epsilon = 1e-06


def test_eeq():
    charge = 0
    mol = acetylene()
    eeq = mol.get_eeq(charge)
    assert (eeq[0] - -0.17166856) < epsilon
    assert (eeq[1] - 0.17166856) < epsilon


def test_eeq_cation():
    charge = 1
    mol = acetylene()
    eeq = mol.get_eeq(charge)
    assert (eeq[0] - 0.59769359) < epsilon
    assert (eeq[1] - 0.40230641) < epsilon


def test_eeq_anion():
    charge = -1
    mol = acetylene()
    eeq = mol.get_eeq(charge)
    assert (eeq[0] - -0.94103071) < epsilon
    assert (eeq[1] - -0.05896929) < epsilon
