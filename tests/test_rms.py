# tests/test_rms.py
import numpy as np
from tests.store import propanolIntermediate
from tests.store import propanolLowest

from kallisto.rmsd import rmsd


def test_rms():
    mol1 = propanolLowest()
    nat1 = mol1.get_number_of_atoms()
    coord1 = mol1.get_positions()
    mol2 = propanolIntermediate()
    coord2 = mol2.get_positions()
    _, u = rmsd(nat1, coord1, coord2)
    assert np.isclose(u[0, 0], 0.98139458)
    assert np.isclose(u[0, 1], -0.04965545)
    assert np.isclose(u[0, 2], -0.18546973)
    assert np.isclose(u[1, 0], 0.06170977)
    assert np.isclose(u[1, 1], 0.9963015)
    assert np.isclose(u[1, 2], 0.05979323)
    assert np.isclose(u[2, 0], 0.18181471)
    assert np.isclose(u[2, 1], -0.07012604)
    assert np.isclose(u[2, 2], 0.98082911)
