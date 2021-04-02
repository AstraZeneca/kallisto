# tests/test_strucreader.py

import os
import tempfile

from kallisto.atom import Atom
import kallisto.reader.strucreader as ksr
from kallisto.reader.turbomole import read as tmreader
from kallisto.reader.xyz import read as xyzreader

# define global lineseperator
s = os.linesep


# turbomole coord files can be read
def test_a_user_can_read_coord():
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
        got = type(atoms[0])
        want = Atom
        assert got is want


def test_user_can_read_xyz():
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz") as f:
        f.write("11" + s)
        f.write("Pyridine" + s)
        f.write("C 1.3603 0.0256 0.0000" + s)
        f.write("C 0.6971 -1.2020 0.0000" + s)
        f.write("C -0.6944 -1.2184 0.0000" + s)
        f.write("C -1.3895 -0.0129 0.0000" + s)
        f.write("C -0.6712 1.1834 0.0000" + s)
        f.write("N 0.6816 1.1960 0.0000" + s)
        f.write("H 2.4530 0.1083 0.0000" + s)
        f.write("H 1.2665 -2.1365 0.0000" + s)
        f.write("H -1.2365 -2.1696 0.0000" + s)
        f.write("H -2.4837 0.0011 0.0000" + s)
        f.write("H -1.1569 2.1657 0.0000" + s)
        f.flush()
        fname = open(f.name)
        atoms = xyzreader(fname)
        fname.close()
        assert len(atoms) == 11


def test_user_can_read_coord():
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        fname = open(f.name)
        atoms = tmreader(fname)
        fname.close()
        assert len(atoms) == 2


def test_ignore_coordinates_after_marker():
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz") as f:
        f.write("test" + s)
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.write("0.00000000000000 0.00000000000000  2.06176434496059 h" + s)
        f.flush()
        fname = open(f.name)
        atoms = tmreader(fname)
        fname.close()
        assert len(atoms) == 2
        position = atoms[0].get("position")
        assert position[0] == 0.0


def test_user_finds_coordnates_in_coord():
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz") as f:
        f.write("test" + s)
        f.write("test" + s)
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        fname = open(f.name)
        atoms = tmreader(fname)
        fname.close()
        assert len(atoms) == 2
