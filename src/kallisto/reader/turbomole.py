# src/kallisto/reader/turbomole.py

from kallisto.atom import Atom


def read(fileObject):
    """Method to read Turbomole coord files.

    coord in Bohr, atom types in lowercase, format:
    $coord
    x y z atomType
    x y z atomType
    $end
    """

    lines = fileObject.readlines()

    atoms = None

    if atoms is None:
        atoms = []

    start = 0

    # find "$coord" marker
    for i, l in enumerate(lines):
        if l.strip().startswith("$coord"):
            start = i
            break
    for line in lines[start + 1 :]:
        # check if new section begins
        if line.startswith("$"):
            break
        else:
            x, y, z, symbolRaw = line.split()[:4]
            symbolShort = symbolRaw.strip()
            atomSymbol = symbolShort[0].upper() + symbolShort[1:].lower()
            atomPosition = [float(x), float(y), float(z)]
            atom = Atom(symbol=atomSymbol, position=atomPosition)
            atoms.append(atom)

    return atoms
