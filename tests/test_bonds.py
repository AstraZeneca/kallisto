# tests/test_bonds.py
from tests.store import ch_radical
from tests.store import pyridine
from tests.store import toluene


def test_bonds_ch_radical():
    mol = ch_radical()
    bonds = mol.get_bonds()
    assert bonds[0] == [1]
    assert bonds[1] == [0]


def test_bonds_pyridine():
    mol = pyridine()
    bonds = mol.get_bonds()
    assert bonds[0] == [1, 5, 6]
    assert bonds[1] == [0, 2, 7]
    assert bonds[4] == [3, 5, 10]


def test_bonds_toluene_partner():
    mol = toluene()
    partner = 1
    bonds = mol.get_bonds(partner)
    assert bonds == [0, 2, 8]
    partner = 5
    bonds = mol.get_bonds(partner)
    assert bonds == [0, 4, 6]
