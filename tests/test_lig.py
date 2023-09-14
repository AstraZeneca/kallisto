# tests/test_lig.py
from tests.store import neopentane

from kallisto.rmsd import recursiveGetSubstructures


def test_lig_neopentane():
    center = 0
    mol = neopentane()
    nat = mol.get_number_of_atoms()
    bonds = mol.get_bonds()
    substructures = recursiveGetSubstructures(nat, bonds, center)  # type: ignore
    assert 1 in substructures[0]
    assert 6 in substructures[0]
    assert 8 in substructures[0]
    assert 12 in substructures[0]
    center = 1
    substructures = recursiveGetSubstructures(nat, bonds, center)  # type: ignore
    assert 6 in substructures[1]
    assert 8 in substructures[2]
    assert 12 in substructures[3]
