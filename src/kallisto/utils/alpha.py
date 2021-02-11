# src/kallisto/utils/alpha.py

import numpy as np


def zeta(a: float, c: float, qref: float, q: float) -> float:
    """Charge scaling function for polarizabilities."""

    if q <= 0:
        zeta = np.exp(a)
    else:
        zeta = np.exp(a * (1.0 - np.exp(c * (1.0 - qref / q))))

    return zeta


def cngw(wf: float, cn: float, cnref: float) -> float:
    """Gaussian weighting factor for reference systems."""

    return np.exp(-wf * (cn - cnref) ** 2)
