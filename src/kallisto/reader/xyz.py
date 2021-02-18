# src/kallisto/reader/xyz.py

from kallisto.atom import Atom
from kallisto.units import Bohr


def read(fileObject):
    """Method to read XYZ files.

    atomic coord in Angstroem, atom types in lowercase, format:
    <number of atoms>
    <comment line>
    <atomType x y z>
    <atomType x y z>
    """

    lines = fileObject.readlines()

    atoms = []

    # by convention the first thing in the xyz is the number of atoms
    nat = int(lines[0])

    # loop over the nat lines
    for line in lines[2 : nat + 2]:
        symbolRaw, x, y, z = line.split()[:4]
        symbolShort = symbolRaw.strip()
        atomSymbol = symbolShort[0].upper() + symbolShort[1:].lower()
        atomPosition = [float(x) / Bohr, float(y) / Bohr, float(z) / Bohr]
        atom = Atom(symbol=atomSymbol, position=atomPosition)
        atoms.append(atom)

    return atoms
