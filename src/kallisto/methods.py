# src/kallisto/methods.py

from typing import Tuple

import numpy as np


def getCoordinationNumbers(
    at: np.ndarray, coords: np.ndarray, cntype: str, threshold: float
):
    """A method to compute coordination numbers (cns).

    CN values are calculated for a given structure and are returned as an
    array. Choose functional type by "cn" defining standard (exp), covalent (cov),
    or error (err)."""

    from kallisto.data import covalent_radius as rcov
    from kallisto.data import pauling_en
    from scipy import special

    nat = len(at)
    cns = np.zeros(shape=(nat,), dtype=np.float64)

    if cntype == "exp":
        k1 = 16.0
        for i in range(nat):
            ia = at[i] - 1
            for j in range(nat):
                if i is j:
                    continue
                dx = coords[j][0] - coords[i][0]
                dy = coords[j][1] - coords[i][1]
                dz = coords[j][2] - coords[i][2]
                rSquared = dx * dx + dy * dy + dz * dz
                if rSquared > threshold:
                    continue
                ja = at[j] - 1
                r = np.sqrt(rSquared)
                rco = rcov[ia] + rcov[ja]
                rr = rco / r
                k = -k1 * (rr - 1.0)
                damp = 1.0 / (1.0 + np.exp(k))
                cns[i] += damp

    if cntype == "erf":
        kn = 7.50
        for i in range(nat):
            ia = at[i] - 1
            for j in range(nat):
                if i is j:
                    continue
                ja = at[j] - 1
                dx = coords[j][0] - coords[i][0]
                dy = coords[j][1] - coords[i][1]
                dz = coords[j][2] - coords[i][2]
                rSquared = dx * dx + dy * dy + dz * dz
                if rSquared > threshold:
                    continue
                r = np.sqrt(rSquared)
                rco = rcov[ia] + rcov[ja]
                damp = 0.5 * (1.0 + special.erf(-kn * (r - rco) / rco))
                cns[i] += damp

    if cntype == "cov":
        # Fitted to match Wiberg bond orders of diatomic molecules
        k4 = 4.10451
        k5 = 19.08857
        k6 = 2 * 11.28174**2

        kn = 7.50
        for i in range(nat):
            ia = at[i] - 1
            for j in range(nat):
                if i is j:
                    continue
                ja = at[j] - 1
                dx = coords[j][0] - coords[i][0]
                dy = coords[j][1] - coords[i][1]
                dz = coords[j][2] - coords[i][2]
                rSquared = dx * dx + dy * dy + dz * dz
                if rSquared > threshold:
                    continue
                r = np.sqrt(rSquared)
                rco = rcov[ia] + rcov[ja]
                eni = pauling_en[ia]
                enj = pauling_en[ja]
                den = k4 * np.exp(-((np.abs(eni - enj) + k5) ** 2) / k6)
                damp = den * 0.5 * (1 + special.erf(-kn * (r - rco) / rco))
                cns[i] += damp

    return cns


def getProximityShells(
    at: np.ndarray, coords: np.ndarray, size: Tuple[int, int], threshold: float
):
    """A method to compute atomic proximity shells (prox)."""

    from kallisto.data import covalent_radius as rcov
    from kallisto.data import pauling_en
    from scipy import special

    nat = len(at)
    prox1 = np.zeros(shape=(nat,), dtype=np.float64)
    prox2 = np.zeros(shape=(nat,), dtype=np.float64)

    # Fitted to match Wiberg bond orders of diatomic molecules
    k4 = 4.10451
    k5 = 19.08857
    k6 = 2 * 11.28174**2

    kn = 7.50
    threshold = 800.0

    # unpack tuple in sizes
    scale1, scale2 = size

    for i in range(nat):
        ia = at[i] - 1
        for j in range(nat):
            if i is j:
                continue
            ja = at[j] - 1
            dx = coords[j][0] - coords[i][0]
            dy = coords[j][1] - coords[i][1]
            dz = coords[j][2] - coords[i][2]
            rSquared = dx * dx + dy * dy + dz * dz
            if rSquared > threshold:
                continue
            r = np.sqrt(rSquared)
            eni = pauling_en[ia]
            enj = pauling_en[ja]
            den = k4 * np.exp(-((np.abs(eni - enj) + k5) ** 2) / k6)

            # smaller border
            rco = scale1 * (rcov[ia] + rcov[ja])
            damp = den * 0.5 * (1 + special.erf(-kn * (r - rco) / rco))
            prox1[i] += damp
            # larger border
            rco = scale2 * (rcov[ia] + rcov[ja])
            damp = den * 0.5 * (1 + special.erf(-kn * (r - rco) / rco))
            prox2[i] += damp

    return prox2 - prox1


def getAtomicPartialCharges(
    at: np.ndarray, coords: np.ndarray, cns: np.ndarray, charge: int
):
    """A method to compute atomic electronegativity equilibration partial
    charges (eeqs).

    EEQ values are calculated for a given structure and are returned as an
    array."""

    from kallisto.data import eeq_alp, eeq_cnfak, eeq_en, eeq_gamm
    from numpy import linalg as LA
    from scipy import special

    # parameter
    sqrt2pi = np.sqrt(2.0 / np.pi)
    nat = len(at)

    # Lagragian space is +1 in dimensionality
    m = nat + 1

    # convert lists to numpy arrays for vectorization
    eeq_alp = np.array(eeq_alp)
    eeq_cnfak = np.array(eeq_cnfak)
    eeq_en = np.array(eeq_en)
    eeq_gamm = np.array(eeq_gamm)

    # setup parameter arrays
    z = at - 1
    xi = eeq_en[z]
    gam = eeq_gamm[z]
    kappa = eeq_cnfak[z]
    alpha = np.power(eeq_alp[z], 2)

    """Set up A matrix and X vector

            αi -> alpha(i), ENi -> xi(i), κi -> kappa(i), Jii -> gam(i)
            γij = 1/√(αi+αj)
            Xi  = -ENi + κi·√CNi
            Aii = Jii + 2/√π·γii
            Aij = erf(γij·Rij)/Rij = 2/√π·F0(γ²ij·R²ij)"""

    # A matrix
    A = np.zeros(shape=(m, m), dtype=np.float64)
    for i in range(nat):
        xyzi = coords[i]
        A[i][i] = gam[i] + sqrt2pi / np.sqrt(alpha[i])
        for j in range(nat):
            if i == j:
                continue
            xyzj = coords[j]
            r = LA.norm(xyzj - xyzi)
            gamij = 1.0 / np.sqrt(alpha[i] + alpha[j])
            A[j][i] = special.erf(gamij * r) / r
            A[i][j] = A[j][i]

    # X vector
    X = np.zeros(shape=(m,), dtype=np.float64)
    X[:nat] = -xi + kappa * np.sqrt(cns)

    # setup Lagragian constraints
    A[:, nat] = 1.0
    A[nat, :] = 1.0
    A[nat, nat] = 0.0
    X[nat] = charge

    # get eeq charges
    qs = np.linalg.solve(A, X)

    return qs[:-1]


def getPolarizabilities(at: np.ndarray, covcn: np.ndarray, qs: np.ndarray, charge: int):
    """A method to compute atomic-charge dependent dynamic atomic polarizabilities (alps).

    ALP values are calculated for a given structure and are returned as an
    array. For the charge dependency EEQ atomic partial charges are used
    in an empirical scaling function as used in the dftd4 program."""

    from kallisto.data.alpha import refx, refh, hcount, ascale, refn
    from kallisto.data.alpha import refcn, refsys, alphaiw, zeff
    from kallisto.data.alpha import sscale, seciw, gam
    from kallisto.utils.alpha import zeta, cngw

    # parameter
    g_a = 3.0
    g_c = 2.0

    nat = len(at)

    # get dimensionality
    ndim = 0
    for i in range(nat):
        ndim += refn[at[i] - 1]

    # index table
    itbl = np.zeros(shape=(7, nat), dtype=np.int64)
    k = 0
    for i in range(nat):
        ia = at[i] - 1
        for ii in range(refn[ia]):
            itbl[ii][i] = k
            k += 1

    # setup ncount and charge scale polarizabilities
    ncount = np.zeros(shape=(7, 86), dtype=np.int64)
    alpha = np.zeros(shape=(23,), dtype=np.float64)
    alphar = np.zeros(shape=(23, 7, 86), dtype=np.float64)

    for i in range(nat):
        cncount = np.zeros(shape=(18,), dtype=int)
        cncount[0] = 1
        ia = at[i] - 1
        for j in range(refn[ia]):
            refis = refsys[j][ia] - 1
            refiz = zeff[refis]
            for jj in range(23):
                alpha[jj] = (
                    sscale[refis]
                    * seciw[jj][refis]
                    * zeta(g_a, g_c * gam[refis], refiz, refh[j][ia] + refiz)
                )
            icn = np.rint(refcn[j][ia])
            cncount[int(icn)] += 1
            for jj in range(23):
                alphar[jj][j][ia] = np.maximum(
                    ascale[j][ia] * (alphaiw[jj][j][ia] - hcount[j][ia] * alpha[jj]),
                    0,
                )
        for j in range(refn[ia]):
            icn = cncount[int(np.rint(refcn[j][ia]))]
            ncount[j][ia] = icn * (icn + 1) / 2

    # weigths
    gw = np.zeros(shape=(ndim,), dtype=np.float64)
    wf = 6.0
    for i in range(nat):
        ia = at[i] - 1
        norm = 0.0
        for ii in range(refn[ia]):
            for iii in range(ncount[ii][ia]):
                twf = (iii + 1) * wf
                norm = norm + cngw(twf, covcn[i], refcn[ii][ia])
        norm = 1.0 / norm
        for ii in range(refn[ia]):
            k = itbl[ii][i]
            for iii in range(ncount[ii][ia]):
                twf = (iii + 1) * wf
                gw[k] += cngw(twf, covcn[i], refcn[ii][ia]) * norm

    # polarizabilities
    zetvec = np.zeros(shape=(ndim), dtype=np.float64)
    aw = np.zeros(shape=(23, nat), dtype=np.float64)
    for i in range(nat):
        ia = at[i] - 1
        iz = zeff[ia]
        for ii in range(refn[ia]):
            k = itbl[ii][i]
            zetvec[k] = gw[k] * zeta(g_a, g_c * gam[ia], refx[ii][ia] + iz, qs[i] + iz)
            for iii in range(23):
                aw[iii][i] += zetvec[k] * alphar[iii][ii][ia]

    atomicAiw = np.zeros(shape=(nat,), dtype=np.float64)
    for i in range(nat):
        atomicAiw[i] = aw[0][i]

    return atomicAiw


def getCovalentBondingPartner(
    at: np.ndarray,
    coords: np.ndarray,
    partner: str,
    thresholdBond: float,
    thresholdCN: float,
):
    """A method to compute an index table for covalent bonding partner.

    thresholdBond defines the treshold for a covalent bond."""

    from kallisto.data import covalent_radius as rcov

    # parameter
    k1 = 16.0

    nat = len(at)

    btbl = np.zeros(shape=(nat, nat), dtype=np.int32)

    for i in range(nat):
        ia = at[i] - 1
        for j in range(nat):
            if i is j:
                continue
            dx = coords[j][0] - coords[i][0]
            dy = coords[j][1] - coords[i][1]
            dz = coords[j][2] - coords[i][2]
            rSquared = dx * dx + dy * dy + dz * dz
            if rSquared > thresholdCN:
                continue
            ja = at[j] - 1
            r = np.sqrt(rSquared)
            rco = rcov[ia] + rcov[ja]
            rr = rco / r
            alpha = -k1 * (rr - 1.0)
            damp = 1.0 / (1.0 + np.exp(alpha))
            if damp > thresholdBond:
                btbl[i][j] = 1

    if partner == "X":
        # Get covalent bonding partners for all atoms
        covalentList = []
        for i in range(nat):
            covalentPartner = []
            k = 0
            for elem in btbl[i][:]:
                if elem == 1:
                    covalentPartner.append(k)
                k += 1
            covalentList.append(covalentPartner)
        return covalentList
    else:
        # Get covalent bonding partners of atom #partner
        covalentPartner = []
        k = 0
        for elem in btbl[int(partner)][:]:
            if elem == 1:
                covalentPartner.append(k)
            k += 1
        return covalentPartner


def getVanDerWaalsRadii(
    nat: int, at: np.ndarray, aw: np.ndarray, vdwtype: str, scale: float
):
    """A method to compute environment and charge dependent van der Waals radii (vdws).

    VDW values are calculated from atomic polarizabilities for a given structure
    and are returned as an array."""

    from kallisto.data import chemical_symbols
    from kallisto.data.vdw import rahm, truhlar

    vdw = np.zeros(shape=(nat,), dtype=np.float64)

    osev = 1.0 / 7.0
    # Empirical scaling and theta_a value from DOI:
    # 10.1103/PhysRevLett.121.183401
    theta_a = 2.54

    if vdwtype == "truhlar":
        # Truhlar: theta_b fitted to match radii from DOI:
        # 10.1021/jp8111556
        for i in range(nat):
            ia = at[i]
            isym = chemical_symbols[ia]
            theta_b = truhlar[isym]
            vdw[i] = scale * theta_b * theta_a * np.power(aw[i], osev)
    elif vdwtype == "rahm":
        # Rahm: theta_b fitted to match radii from DOI:
        # 10.1002/chem.201700610
        for i in range(nat):
            ia = at[i]
            isym = chemical_symbols[ia]
            theta_b = rahm[isym]
            vdw[i] = scale * theta_b * theta_a * np.power(aw[i], osev)

    return vdw
