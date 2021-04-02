# tests/test_stm.py

from kallisto.sterics import getClassicalSterimol
from tests.store import toluene

# global epsilon
epsilon = 1e-06


def test_stm():
    mol = toluene()
    origin = 6
    partner = 5
    L, bmin, bmax = getClassicalSterimol(mol, origin, partner)
    assert (L - 12.714385) < epsilon
    assert (bmin - 3.539068) < epsilon
    assert (bmax - 6.640342) < epsilon
