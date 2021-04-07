# tests/conftest.py

import os

import pytest

from tests.store import ch_radical, iridiumCatalyst, pyridine

s = os.linesep


@pytest.fixture(scope="session")
def pyridine_xyz(tmpdir_factory):
    """Write a pyridine molecule to a data file in a temporary directory."""
    name = "pyridine.xyz"
    dn = tmpdir_factory.mktemp("test_pyridine", numbered=True)
    fn = dn.join(name)
    dn = str(dn)
    fn = str(fn)

    mol = pyridine()
    mol.writeMolecule(name=name, path=dn)

    return fn


@pytest.fixture(scope="session")
def ch_radical_xyz(tmpdir_factory):
    """Write a ch radical molecule to a data file in a temporary directory."""
    name = "ch_radical.xyz"
    dn = tmpdir_factory.mktemp("test_ch_radical", numbered=True)
    fn = dn.join(name)
    dn = str(dn)
    fn = str(fn)

    mol = ch_radical()
    mol.writeMolecule(name=name, path=dn)

    return fn


@pytest.fixture(scope="session")
def iridiumcat_xyz(tmpdir_factory):
    """Write an Iridium catalyst to a data file in a temporary directory."""
    name = "iridium.xyz"
    dn = tmpdir_factory.mktemp("test_iridium", numbered=True)
    fn = dn.join(name)
    dn = str(dn)
    fn = str(fn)

    mol = iridiumCatalyst()
    mol.writeMolecule(name=name, path=dn)

    return fn


@pytest.fixture(scope="session")
def fluoromethane_coord(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data").join("fluoromethane.coord")

    with open(fn, mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("  1.87167924 -0.101043656  0.1596818582  c" + s)
        f.write("  4.43543289 -0.101043656  0.1596818582  f" + s)
        f.write("  1.20847986 -1.386321988  1.6312493924  h" + s)
        f.write("  1.20847986 -0.732816897 -1.6891694984  h" + s)
        f.write("  1.20846096  1.816007916  0.5369845779  h" + s)
        f.write("$end")
        f.flush()

    return fn


@pytest.fixture(scope="session")
def lithium_hydride_coord(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data").join("lithium_hydride.coord")

    with open(fn, mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0 0 0 li" + s)
        f.write("0 0 2 h" + s)
        f.write("0 0 4 h" + s)
        f.write("$end")
        f.flush()

    return fn


@pytest.fixture(scope="session")
def ch_coord(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data").join("ch.coord")

    with open(fn, mode="w+", encoding="utf-8") as f:
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()

    return fn


@pytest.fixture(scope="session")
def ch_coord_2(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data").join("ch-2.coord")

    with open(fn, mode="w+", encoding="utf-8") as f:
        f.write("test" + s)
        f.write("test" + s)
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.flush()

    return fn


@pytest.fixture(scope="session")
def ch_coord_ignore(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data").join("ch_ignore.coord")

    with open(fn, mode="w+", encoding="utf-8") as f:
        f.write("test" + s)
        f.write("$coord" + s)
        f.write("0.00000000000000 0.00000000000000 -1.06176434496059 c" + s)
        f.write("0.00000000000000 0.00000000000000  1.06176434496059 h" + s)
        f.write("$end")
        f.write("0.00000000000000 0.00000000000000  2.06176434496059 h" + s)
        f.flush()

    return fn
