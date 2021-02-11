# tests/test_molecule.py

import os
import tempfile

import numpy as np

from kallisto.molecule import Molecule
import kallisto.reader.strucreader as ksr

# define global lineseperator
s = os.linesep


# turbomole coord files can be read
def test_a_user_can_compute_coordination_numbers_for_molecule():
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("  1.87167924 -0.101043656  0.1596818582  c" + s)
        f.write("  4.43543289 -0.101043656  0.1596818582  f" + s)
        f.write("  1.20847986 -1.386321988  1.6312493924  h" + s)
        f.write("  1.20847986 -0.732816897 -1.6891694984  h" + s)
        f.write("  1.20846096  1.816007916  0.5369845779  h" + s)
        f.write("$end")
        f.flush()
        fileObject = open(f.name, "r+")
        atoms = ksr.read(fileObject)
        molecule = Molecule(symbols=atoms)
        cns = molecule.get_cns(cntype="exp")
        c1 = np.around(cns, decimals=1)
        got = c1[0]
        want = 4.0
        assert got == want


def test_user_can_calculate_eeq_atomic_charges():
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0 0 0 li" + s)
        f.write("0 0 2 h" + s)
        f.write("0 0 4 h" + s)
        f.write("$end")
        f.flush()
        fileObject = open(f.name, "r+")
        atoms = ksr.read(fileObject)
        molecule = Molecule(symbols=atoms)
        eeq = molecule.get_eeq(charge=0)
        want = [0.51925854, -0.35007273, -0.16918582]
        difference = sum([a - b for a, b in zip(want, eeq)])
        assert difference < 1e-6
