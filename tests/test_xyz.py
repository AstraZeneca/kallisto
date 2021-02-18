# tests/test_xyz.py

import os
import tempfile

from kallisto.reader.xyz import read

# define global lineseperator
s = os.linesep


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
        atoms = read(fname)
        fname.close()
        assert len(atoms) == 11
