# tests/test_cns.py

from tests.store import acetylene

# global epsilon
epsilon = 1e-06


def test_cns_exp_acetylene():
    cntype = "exp"
    mol = acetylene()
    cns = mol.get_cns(cntype)
    assert (cns[0] - 0.9867892) < epsilon
    assert (cns[1] - 0.9867892) < epsilon


def test_cns_erf_acetylene():
    cntype = "erf"
    mol = acetylene()
    cns = mol.get_cns(cntype)
    assert (cns[0] - 0.9878465) < epsilon
    assert (cns[1] - 0.9878465) < epsilon


def test_cns_cov_acetylene():
    cntype = "cov"
    mol = acetylene()
    cns = mol.get_cns(cntype)
    assert (cns[0] - 0.9878465) < epsilon
    assert (cns[1] - 0.9878465) < epsilon
