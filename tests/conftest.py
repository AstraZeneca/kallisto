# tests/conftest.py

import pytest

from tests.store import acetylene, iridiumCatalyst, pyridine, pyridine_minusH


@pytest.fixture(scope="session")
def pyridine_xyz(tmpdir_factory):
    """Write a pyridine molecule to a data file in a temporary directory."""
    dn = tmpdir_factory.mktemp("test_pyridine", numbered=True)
    fn = dn.join("pyridine.xyz")
    dn = str(dn)
    fn = str(fn)

    mol = pyridine()
    mol.writeMolecule(name="pyridine.xyz", path=dn)

    return fn


@pytest.fixture(scope="session")
def pyridine_minusH_xyz(tmpdir_factory):
    """Write a pyridine molecule to a data file in a temporary directory."""
    dn = tmpdir_factory.mktemp("test_pyridine_minusH", numbered=True)
    fn = dn.join("pyridine_minusH.xyz")
    dn = str(dn)
    fn = str(fn)

    mol = pyridine_minusH()
    mol.writeMolecule(name="pyridine_minusH.xyz", path=dn)

    return fn


@pytest.fixture(scope="session")
def acetylene_xyz(tmpdir_factory):
    """Write a pyridine molecule to a data file in a temporary directory."""
    dn = tmpdir_factory.mktemp("test_acetylene", numbered=True)
    fn = dn.join("acetylene.xyz")
    dn = str(dn)
    fn = str(fn)

    mol = acetylene()
    mol.writeMolecule(name="acetylene.xyz", path=dn)

    return fn


@pytest.fixture(scope="session")
def iridiumcat_xyz(tmpdir_factory):
    """Write a pyridine molecule to a data file in a temporary directory."""
    dn = tmpdir_factory.mktemp("test_iridium", numbered=True)
    fn = dn.join("iridium.xyz")
    dn = str(dn)
    fn = str(fn)

    mol = iridiumCatalyst()
    mol.writeMolecule(name="iridium.xyz", path=dn)

    return fn
