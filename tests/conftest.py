# tests/conftest.py

import pytest

from tests.store import acetylene, iridiumCatalyst, pyridine


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
def acetylene_xyz(tmpdir_factory):
    """Write an acetylene molecule to a data file in a temporary directory."""
    name = "acetylene.xyz"
    dn = tmpdir_factory.mktemp("test_acetylene", numbered=True)
    fn = dn.join(name)
    dn = str(dn)
    fn = str(fn)

    mol = acetylene()
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
