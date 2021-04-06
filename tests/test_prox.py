# tests/test_prox.py

from tests.store import toluene

# global epsilon
epsilon = 1e-06


def test_prox():
    mol = toluene()
    size = (2, 3)
    prox = mol.get_prox(size)
    assert (prox[0] - 4.381228496800993) < epsilon
    assert (prox[1] - 3.3677521592350548) < epsilon
    assert (prox[2] - 2.8495247858657446) < epsilon


def test_prox_different_size():
    mol = toluene()
    size = (1, 2)
    prox = mol.get_prox(size)
    assert (prox[0] - 5.841080507811956) < epsilon
    assert (prox[1] - 5.501735547766682) < epsilon
    assert (prox[2] - 5.037814140144763) < epsilon
