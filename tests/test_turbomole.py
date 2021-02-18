# tests/test_turbomole.py

import os
import tempfile

from kallisto.reader.turbomole import read

# define global lineseperator
s = os.linesep


def test_user_can_read_coord():
    with tempfile.NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".xyz") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()
        fname = open(f.name)
        atoms = read(fname)
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
        atoms = read(fname)
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
        atoms = read(fname)
        fname.close()
        assert len(atoms) == 2
