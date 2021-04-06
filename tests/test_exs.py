# tests/test_exs.py

import os

from kallisto.rmsd import exchangeSubstructure
from kallisto.units import Bohr
from tests.store import iridiumCatalyst, pyridine_mH

# global epsilon
epsilon = 1e-4


def test_exs():
    # Iridium atom complex
    center = 18
    # Exchange benzene with pyridine
    subnr = 2
    # set output name
    name = "iridium_pyridine"
    # define rotation angle along covalent bond
    rotate = 0
    # do not constrain bonding atom of substrate
    exclude = False

    ref = iridiumCatalyst()
    refBonds = ref.get_bonds()
    refNat = ref.get_number_of_atoms()
    exchanger = pyridine_mH()
    exchangerBonds = exchanger.get_bonds()

    mol = exchangeSubstructure(
        refNat,
        center,
        subnr,
        refBonds,
        ref,
        exchanger,
        exchangerBonds,
        name,
        rotate,
        exclude,
    )  # type: ignore

    at = mol.get_atomic_numbers()
    assert at[85] == 6
    coords = mol.get_positions()
    # position of Carbon bonding to Iridium
    assert (coords[85, 0] - 0.2001 / Bohr) < epsilon
    assert (coords[85, 1] - -2.2076 / Bohr) < epsilon
    assert (coords[85, 2] - 2.6787 / Bohr) < epsilon
    # position of Nitrogen
    assert (coords[90, 0] - -1.0435 / Bohr) < epsilon
    assert (coords[90, 1] - -2.4569 / Bohr) < epsilon
    assert (coords[90, 2] - 3.1496 / Bohr) < epsilon
    constrain = "constrain.inp"
    gotFile = os.path.isfile(constrain)
    assert gotFile is True
    if gotFile:
        os.remove(constrain)
    name += ".xyz"
    gotFile = os.path.isfile(name)
    assert gotFile is True
    if gotFile:
        os.remove(name)


def test_exs_rotate_180_degrees():
    # Iridium atom complex
    center = 18
    # Exchange benzene with pyridine
    subnr = 2
    # set output name
    name = "iridium_pyridine_180_rotated"
    # define rotation angle along covalent bond
    rotate = 180
    # do not constrain bonding atom of substrate
    exclude = False

    ref = iridiumCatalyst()
    refBonds = ref.get_bonds()
    refNat = ref.get_number_of_atoms()
    exchanger = pyridine_mH()
    exchangerBonds = exchanger.get_bonds()

    mol = exchangeSubstructure(
        refNat,
        center,
        subnr,
        refBonds,  # type: ignore
        ref,
        exchanger,
        exchangerBonds,  # type: ignore
        name,
        rotate,
        exclude,
    )  # type: ignore

    at = mol.get_atomic_numbers()
    assert at[85] == 6
    coords = mol.get_positions()
    # position of Carbon bonding to Iridium
    assert (coords[85, 0] - 0.2001 / Bohr) < epsilon
    assert (coords[85, 1] - -2.2076 / Bohr) < epsilon
    assert (coords[85, 2] - 2.6787 / Bohr) < epsilon
    # position of Nitrogen
    assert (coords[90, 0] - 1.3335 / Bohr) < epsilon
    assert (coords[90, 1] - -2.8374 / Bohr) < epsilon
    assert (coords[90, 2] - 3.0650 / Bohr) < epsilon
    constrain = "constrain.inp"
    gotFile = os.path.isfile(constrain)
    assert gotFile is True
    if gotFile:
        os.remove(constrain)
    name += ".xyz"
    gotFile = os.path.isfile(name)
    assert gotFile is True
    if gotFile:
        os.remove(name)
