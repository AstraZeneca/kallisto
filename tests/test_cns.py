# tests/test_cns.py
import numpy as np
from tests.store import ch_radical


def test_cns_exp_ch_radical():
    cntype = "exp"
    mol = ch_radical()
    cns = mol.get_cns(cntype)
    assert np.isclose(cns[0], 0.9867892708145132)
    assert np.isclose(cns[1], 0.9867892708145132)


def test_cns_erf_ch_radical():
    cntype = "erf"
    mol = ch_radical()
    cns = mol.get_cns(cntype)
    assert np.isclose(cns[0], 0.9878465867658608)
    assert np.isclose(cns[1], 0.9878465867658608)


def test_cns_cov_ch_radical():
    cntype = "cov"
    mol = ch_radical()
    cns = mol.get_cns(cntype)
    assert np.isclose(cns[0], 0.9189476178185281)
    assert np.isclose(cns[1], 0.9189476178185281)
