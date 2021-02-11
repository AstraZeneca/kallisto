# src/kallisto/reader/strucreader.py

from kallisto.atom import Atom
from kallisto.molecule import Molecule
import kallisto.reader.turbomole as tm
import kallisto.reader.xyz as xyz


def constructMolecule(geometry: str) -> Molecule:
    """Helper function to construct a Molecule."""
    try:
        with open(geometry, "r+") as fileObject:
            # read atoms from input file
            atoms = read(fileObject)

            # create molecule from atoms
            molecule = Molecule(symbols=atoms)
    except FileNotFoundError:
        print(" Error: Input file '" + geometry + "' not found.")
        exit()

    return molecule


def construct_atom():
    return Atom()


def read(fileObject):
    """Method to first check the file type and then read

    the structure accordingly
    The returned atom coordinates will be in Bohr
    """

    fname = fileObject.name.lower()
    filetp = "unknown"

    lines = fileObject.readlines()

    if fname.endswith((".xyz")):
        filetp = "xyz"
    else:
        for line in lines:
            if line.strip().startswith("$coord"):
                filetp = "turbomole"
                break

    fileObject.close()

    atoms = None

    if atoms is None:
        atoms = []

    newfile = open(fname, "r+")

    if filetp == "turbomole":
        atoms = tm.read(newfile)
    elif filetp == "xyz":
        atoms = xyz.read(newfile)
    else:
        print("Error: file type of input '" + fname + "' not recognized")
        quit()

    newfile.close()

    return atoms
