# tests/test_molecule.py

import os
import tempfile

import numpy as np
import pytest

from kallisto.atom import Atom
from kallisto.molecule import Molecule
import kallisto.reader.strucreader as ksr

# define global lineseperator
s = os.linesep

# Create global molecule
atoms = []
for i in range(3):
    atoms.append(Atom(symbol="H", position=(0, 0, i)))


# turbomole coord files can be read
def test_a_user_can_compute_coordination_numbers_for_molecule(fluoromethane_coord):
    fileObject = open(fluoromethane_coord, "r+")
    atoms = ksr.read(fileObject)
    molecule = Molecule(symbols=atoms)
    cns = molecule.get_cns(cntype="exp")
    c1 = np.around(cns, decimals=1)
    got = c1[0]
    want = 4.0
    assert got == want


def test_user_can_calculate_eeq_atomic_charges(lithium_hydride_coord):
    fileObject = open(lithium_hydride_coord, "r+")
    atoms = ksr.read(fileObject)
    molecule = Molecule(symbols=atoms)
    eeq = molecule.get_eeq(charge=0)
    want = [0.51925854, -0.35007273, -0.16918582]
    difference = sum([a - b for a, b in zip(want, eeq)])
    assert difference < 1e-6


def test_user_can_pass_information_from_mol_to_mol():
    reference = Molecule(atoms)
    ref = reference.get_number_of_atoms()
    molecule = Molecule(symbols=reference)
    nat = molecule.get_number_of_atoms()
    assert nat == ref
    positions = molecule.get_positions()
    assert positions[0][2] == 0
    assert positions[1][2] == 1
    assert positions[2][2] == 2
    numbers = molecule.get_array("numbers")
    assert numbers[0] == 1


def test_user_can_set_an_array():
    reference = Molecule(atoms)
    reference.set_array("numbers", [1, 2, 3])
    numbers = reference.get_array("numbers")
    assert numbers[0] == 1
    assert numbers[1] == 2
    assert numbers[2] == 3
    reference.set_array("numbers", None)
    reference.set_array("numbers", None)
    with pytest.raises(Exception) as error:
        reference.get_array("numbers")
    assert error.value.args[0] == "numbers"
    reference.set_array("numbers", [1, 2, 3], None)
    numbers = reference.get_array("numbers")
    assert numbers[0] == 1


def test_user_can_create_mol_without_positions():
    molecule = Molecule(symbols=1, positions=None)
    assert len(molecule.get_positions()) == 0


def test_user_cannot_create_existing_array():
    reference = Molecule(atoms)
    reference.new_array("symbol", ["H", "H", "C"])
    with pytest.raises(Exception) as error:
        reference.new_array("symbol", ["H", "H", "C"])
    assert error.value.args[0] == 'Array "symbol" already present'


def test_user_cannot_create_array_with_different_size():
    reference = Molecule(atoms)
    with pytest.raises(Exception) as error:
        reference.new_array("symbol", [])
    assert error.value.args[0] == 'Array "symbol" has wrong length: 0 != 3.'


def test_user_cannot_set_array_with_different_sizes():
    reference = Molecule(atoms)
    with pytest.raises(Exception) as error:
        reference.set_array("symbols", [1, 2, 3, 4])
    assert "4 != 3" in error.value.args[0]


def test_user_can_copy_an_array():
    reference = Molecule(atoms)
    numbers = reference.get_array("numbers", True)
    assert numbers[0] == 1


def test_user_can_get_atomic_numbers():
    reference = Molecule(atoms)
    numbers = reference.get_atomic_numbers()
    assert numbers[0] == 1


def test_user_can_copy_a_mol():
    reference = Molecule(atoms)
    new = reference.copy()
    numbers = new.get_atomic_numbers()
    assert numbers[0] == 1


def test_user_can_get_number_of_atoms():
    reference = Molecule(atoms)
    nat = reference.get_number_of_atoms()
    assert nat == 3
