# tests/test_rms.py

from kallisto.rmsd import rmsd
from tests.store import propanolIntermediate, propanolLowest

# global epsilon
epsilon = 1e-06


def test_rms():
    mol1 = propanolLowest()
    nat1 = mol1.get_number_of_atoms()
    coord1 = mol1.get_positions()
    mol2 = propanolIntermediate()
    coord2 = mol2.get_positions()
    error, u = rmsd(nat1, coord1, coord2)
    assert (error - 1.12070194) < epsilon
    assert (u[0, 0] - 0.98139458) < epsilon
    assert (u[0, 1] - -0.04965545) < epsilon
    assert (u[0, 2] - -0.18546973) < epsilon
    assert (u[1, 0] - 0.06170977) < epsilon
    assert (u[1, 1] - 0.9963015) < epsilon
    assert (u[1, 2] - 0.05979323) < epsilon
    assert (u[2, 0] - 0.18181471) < epsilon
    assert (u[2, 1] - -0.07012604) < epsilon
    assert (u[2, 2] - 0.98082911) < epsilon
