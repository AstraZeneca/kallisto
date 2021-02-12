# test/test_atom.py

from kallisto.atom import Atom
from kallisto.molecule import Molecule


def test_user_can_create_an_atom():
    symbol = "C"
    atomPosition = [0, 0, 0]
    atom = Atom(symbol=symbol, position=atomPosition)
    assert atom.symbol == "C"


def test_user_can_set_an_atom_symbol():
    symbol = "C"
    atomPosition = [0, 0, 0]
    atom = Atom(symbol=symbol, position=atomPosition)
    atom.set("symbol", "N")
    assert atom.symbol == "N"


def test_user_can_get_an_atom_symbol():
    atom = Atom(symbol="C", position=[0, 0, 0])
    got = atom.get("symbol")
    want = "C"
    assert got == want


def test_user_can_create_atom_from_element_string():
    atom = Atom("H")
    assert atom.symbol == "H"


def test_user_can_modify_atom_via_element_number():
    atom = Atom("H")
    atom.number = 6
    assert atom.symbol == "C"


def test_molecule_is_not_none():
    atoms = []
    for i in range(2):
        atoms.append(Atom(symbol="H", position=(0, 0, i)))
    molecule = Molecule(atoms)
    atom = Atom(symbol="H", position=(0, 0, 4), molecule=molecule)
    assert atom.molecule is not None
