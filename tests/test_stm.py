# tests/test_stm.py
import numpy as np
from tests.store import toluene

from kallisto.sterics import getClassicalSterimol


def test_stm():
    mol = toluene()
    origin = 6
    partner = 5
    L, bmin, bmax = getClassicalSterimol(mol, origin, partner)
    assert np.isclose(L, 12.714385)
    assert np.isclose(bmin, 3.539068)
    assert np.isclose(bmax, 6.640342)
