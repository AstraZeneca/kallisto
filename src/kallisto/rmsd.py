# src/kallisto/rmsd.py

from collections import Counter
from typing import Tuple

import numpy as np

from kallisto.atom import Atom
from kallisto.molecule import Molecule


def rmsd(n: int, coord1: np.ndarray, coord2: np.ndarray) -> Tuple[float, np.ndarray]:
    """Calculate the least square rmsd in Angstrom for two
    coordinate sets coord1(n,3) and coord2(n,3) using a method based on
    quaternions."""

    from kallisto.units import Bohr

    # copy original coordinates
    x = np.zeros(shape=(n, 3), dtype=np.float64)
    y = np.zeros(shape=(n, 3), dtype=np.float64)

    x[0:n, :] = coord1[0:n, :]
    y[0:n, :] = coord2[0:n, :]

    # calculate the barycenters, centroidal coordinates, and the norms
    x_norm = 0.0
    y_norm = 0.0
    x_center = np.zeros(shape=(3,), dtype=np.float64)
    y_center = np.zeros(shape=(3,), dtype=np.float64)
    xi = np.zeros(shape=(n,), dtype=np.float64)
    yi = np.zeros(shape=(n,), dtype=np.float64)

    for i in range(3):
        for j in range(n):
            xi[j] = x[j, i] * Bohr
            yi[j] = y[j, i] * Bohr
        x_center[i] = np.sum(xi) / float(n)
        y_center[i] = np.sum(yi) / float(n)
        xi[:] = xi[:] - x_center[i]
        yi[:] = yi[:] - y_center[i]
        x[:, i] = xi[:]
        y[:, i] = yi[:]
        x_norm += np.dot(xi, xi)
        y_norm += np.dot(yi, yi)

    # calculate the R matrix
    rmat = np.zeros(shape=(3, 3), dtype=np.float64)
    for i in range(3):
        for j in range(3):
            rmat[i, j] = np.dot(x[:, i], y[:, j])

    # calculate the S matrix (quaternions)
    fac = -1
    smat = np.zeros(shape=(4, 4), dtype=np.float64)
    smat[0, 0] = rmat[0, 0] + rmat[1, 1] + rmat[2, 2]
    smat[1, 0] = fac * (rmat[1, 2] - rmat[2, 1])
    smat[2, 0] = fac * (rmat[2, 0] - rmat[0, 2])
    smat[3, 0] = fac * (rmat[0, 1] - rmat[1, 0])

    smat[0, 1] = smat[1, 0]
    smat[1, 1] = rmat[0, 0] - rmat[1, 1] - rmat[2, 2]
    smat[2, 1] = rmat[0, 1] + rmat[1, 0]
    smat[3, 1] = rmat[0, 2] + rmat[2, 0]

    smat[0, 2] = smat[2, 0]
    smat[1, 2] = smat[2, 1]
    smat[2, 2] = -rmat[0, 0] + rmat[1, 1] - rmat[2, 2]
    smat[3, 2] = rmat[1, 2] + rmat[2, 1]

    smat[0, 3] = smat[3, 0]
    smat[1, 3] = smat[3, 1]
    smat[2, 3] = smat[3, 2]
    smat[3, 3] = -rmat[0, 0] - rmat[1, 1] + rmat[2, 2]

    # calculate largest eigenvalue and eigenvector
    eigenval = 0.0
    eigenvec = np.zeros(shape=(4,), dtype=np.float64)
    eigenval, eigenvec = dstmev(smat)

    # convert quaternion eigenvec to rotation matrix U
    u = np.zeros(shape=(3, 3), dtype=np.float64)
    u = rotationMatrix(eigenvec)

    # root mean squared deviation
    error = 0.0
    error = np.sqrt(np.maximum(0.0, ((x_norm + y_norm) - 2 * eigenval) / float(n),))

    return error, u


def rotationMatrix(q: np.ndarray) -> np.ndarray:
    """Constructs rotation matrix U from quaternion q."""
    q0 = q[0]
    q1 = q[1]
    q2 = q[2]
    q3 = q[3]

    b0 = 2.0 * q0
    b1 = 2.0 * q1
    b2 = 2.0 * q2
    b3 = 2.0 * q3

    q00 = b0 * q0 - 1.0
    q01 = b0 * q1
    q02 = b0 * q2
    q03 = b0 * q3

    q11 = b1 * q1
    q12 = b1 * q2
    q13 = b1 * q3

    q22 = b2 * q2
    q23 = b2 * q3

    q33 = b3 * q3

    u = np.zeros(shape=(3, 3), dtype=np.float64)

    u[0, 0] = q00 + q11
    u[0, 1] = q12 - q03
    u[0, 2] = q13 + q02

    u[1, 0] = q12 + q03
    u[1, 1] = q00 + q22
    u[1, 2] = q23 - q01

    u[2, 0] = q13 - q02
    u[2, 1] = q23 + q01
    u[2, 2] = q00 + q33

    return u


def dstmev(smat: np.ndarray):
    """Compute the leading eigenvalue and eigenvector of a symmetric,
    traceless 4x4 matrix A by an inverse power iteration:
    1) convert matrix to tridiagonal form by rotations
     V*A*V' = T
    2) use Gershgorin's theorem to estimate a lower bound
    for the leading negative eigenvalue:
     lambda > g=min(T00-t01,-t10+T11-t12,-t21+T22-t23,-t23+T33)
     where tij=abs(Tij)
    3) create positive definite matrix
     B = T-gI
    4) apply single value decomposition to compute eigenvalues and
    eigenvectors for SPD matrix B
    5) shift the spectrum back and keep leading singular vector and largest
    eigenvalue
    6) convert eigenvalue to original matrix A by multiplication with V'"""

    # convert triogonal form
    tmat, vmat = givens4(smat)

    # 2) estimate lower bond of smallest eigenvalue by Gershgorin's theorem
    eigenval = np.amin(
        [
            tmat[0, 0] - abs(tmat[0, 1]),
            -abs(tmat[1, 0]) + tmat[1, 1] - abs(tmat[1, 2]),
            -abs(tmat[2, 1]) + tmat[2, 2] - abs(tmat[2, 3]),
            -abs(tmat[3, 2]) + tmat[3, 3],
        ]
    )

    # 3) form positive definite matrix
    #    T = lambda * I - T
    for i in range(4):
        tmat[i, i] = tmat[i, i] - eigenval

    # 4) compute singular values/vectors of SPD matrix
    from scipy import linalg

    U, s, Vh = linalg.svd(tmat)

    # 5) get maximum eigenvalue and shift spectrum back
    max_loc = np.argmax(s)
    eigenval += s[max_loc]

    # 6) convert eigenvalue to original matrix A
    eigenvec = np.matmul(vmat, Vh[max_loc])

    return eigenval, eigenvec


def givens4(A: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Perform givens rotations to reduce the
    symmetric 4x4 matrix to tridiagonal form."""

    tmat = np.zeros(shape=(4, 4), dtype=np.float64)
    vmat = np.zeros(shape=(4, 4), dtype=np.float64)

    # initialize
    tmat = A

    # zero out entries [3,0] and [0,3]
    # compute cos and sin of rotation angle in the 3-4 plane
    r1 = pythag(tmat[2, 0], tmat[3, 0])
    if r1 != 0.0:
        c1 = tmat[2, 0] / r1
        s1 = tmat[3, 0] / r1
        vmat[2, 2] = c1
        vmat[2, 3] = s1
        vmat[3, 2] = -s1
        vmat[3, 3] = c1
        tmat[2, 0] = r1
        tmat[3, 0] = 0.0
        tmat[2:, 1:] = np.matmul(vmat[2:, 2:], tmat[2:, 1:])
        tmat[0:2, 2:] = np.transpose(tmat[2:, 0:2])
        tmat[2:, 2:] = np.matmul(tmat[2:, 2:], np.transpose(vmat[2:, 2:]))
    else:
        c1 = 1.0
        s1 = 0.0

    # zero out entries [2,0] and [0,2]
    # compute cos and sin of rotation angle in the 2-3 plane
    r2 = pythag(tmat[2, 0], tmat[1, 0])
    if r2 != 0.0:
        c2 = tmat[1, 0] / r2
        s2 = tmat[2, 0] / r2
        vmat[1, 1] = c2
        vmat[1, 2] = s2
        vmat[2, 1] = -s2
        vmat[2, 2] = c2
        tmat[1, 0] = r2
        tmat[2, 0] = 0.0
        tmat[1:3, 1:4] = np.matmul(vmat[1:3, 1:3], tmat[1:3, 1:4])
        tmat[0, 1:3] = tmat[1:3, 0]
        tmat[3, 1:3] = tmat[1:3, 3]
        tmat[1:3, 1:3] = np.matmul(tmat[1:3, 1:3], np.transpose(vmat[1:3, 1:3]))
    else:
        c2 = 1.0
        s2 = 0.0

    # zero out entries [3,1] and [1,3]
    # compute cos and sin of rotation angle in the 3-4 plane
    r3 = pythag(tmat[3, 1], tmat[2, 1])
    if r3 != 0.0:
        c3 = tmat[2, 1] / r3
        s3 = tmat[3, 1] / r3
        vmat[2, 2] = c3
        vmat[2, 3] = s3
        vmat[3, 2] = -s3
        vmat[3, 3] = c3
        tmat[2, 1] = r3
        tmat[3, 1] = 0.0
        tmat[2:, 2:] = np.matmul(vmat[2:, 2:], tmat[2:, 2:])
        tmat[0:3, 2:] = np.transpose(tmat[2:, 0:3])
        tmat[2:, 2:] = np.matmul(tmat[2:, 2:], np.transpose(vmat[2:, 2:]))
    else:
        c3 = 1.0
        s3 = 0.0

    # compute net rotation matrix (accumulate similarity for evec. computation)
    # to save transposing later, This is the transpose!
    vmat[0, 0] = 1.0
    vmat[0, 1:] = 0.0
    vmat[1:, 0] = 0.0
    vmat[1, 1] = c2
    vmat[2, 1] = c1 * s2
    vmat[3, 1] = s1 * s2

    c1c2 = c1 * c2
    s1c2 = s1 * c2

    vmat[1, 2] = -s2 * c3
    vmat[2, 2] = c1c2 * c3 - s1 * s3
    vmat[3, 2] = s1c2 * c3 + c1 * s3
    vmat[1, 3] = s2 * s3
    vmat[2, 3] = -c1c2 * s3 - s1 * c3
    vmat[3, 3] = -s1c2 * s3 + c1 * c3

    return tmat, vmat


def pythag(a: np.ndarray, b: np.ndarray) -> float:

    aAbs = abs(a)
    bAbs = abs(b)
    result = 0.0

    if aAbs > bAbs:
        frac = bAbs / aAbs
        result = aAbs * np.sqrt(1 + np.power(frac, 2))
    else:
        if bAbs == 0.0:
            result = 0.0
        else:
            frac = aAbs / bAbs
            result = bAbs * np.sqrt(1 + np.power(frac, 2))

    return result


def recursiveGetSubstructures(n: int, bonds: np.ndarray, center: int):
    """Recursively call getSubstructures to get all covalent substructres of center atom."""

    paths = []

    for i in range(len(bonds[center])):
        path = np.zeros(shape=(n,), dtype=np.int32)
        path.fill(-1)
        partner = bonds[center][i]
        count = Counter()  # type: ignore
        count["step"] = -1
        getSubstructures(n, bonds, center, partner, path, count=count)
        # get rid of -1 elements in path
        path = [x for x in path if x != -1]
        paths.append(path)

    return paths


def getSubstructures(
    n: int,
    bonds: np.ndarray,
    center: int,
    partner: int,
    path: np.ndarray,
    count: Counter,
):
    """Recursively determine all atoms that belong to the ligands of center atom."""

    # partner atom matches center
    if center == partner:
        return
    # actual atom is already in path
    if partner in path:
        return
    # count is larger than number of atoms
    if (count["step"] + 1) > (n - 1):
        return

    count["step"] += 1
    path[count["step"]] = partner

    if len(bonds[partner]) > 1:
        for i in range(len(bonds[partner])):
            p = bonds[partner][i]
            # enter recursive call
            getSubstructures(n, bonds, center, p, path, count=count)


def exchangeSubstructure(
    n: int,
    center: int,
    subnr: int,
    bonds: np.ndarray,
    ref: Molecule,
    newsub: Molecule,
    newSubBonds: np.ndarray,
    name: str,
    rotate: int,
    exclude: bool,
):
    """Exchange substructure (subnr) from ref with substrate."""

    # setup initial complex
    refxyz = ref.get_positions()
    refat = ref.get_atomic_numbers()
    refnat = ref.get_number_of_atoms()

    # match new and old substrate
    newat = newsub.get_atomic_numbers()
    newnat = newsub.get_number_of_atoms()

    # extract coordinates of central atom
    centralAtom = np.zeros(shape=(1, 3), dtype=np.float64)
    centralAtom = refxyz[center, :]

    for i in range(len(bonds[center])):
        if i == subnr:
            path = np.zeros(shape=(n,), dtype=np.int32)
            # set all elements to -1
            path.fill(-1)
            partner = bonds[center][i]
            count = Counter()  # type: ignore
            count["step"] = -1
            getSubstructures(n, bonds, center, partner, path, count=count)
            # get rid of -1 elements
            path = [x for x in path if x != -1]
            path = np.array(path)
            # create molecule from given path
            oldsub = getSubstructureFromPath(ref, path)
            # get all bonding partner
            oldSubBonds = oldsub.get_bonds(partner="X")
            outxyz = np.zeros(shape=(newnat, 3), dtype=np.float64)
            outxyz = matchSubstrates(
                bonds, newsub, newSubBonds, oldsub, oldSubBonds, centralAtom,
            )

            # atoms from complex excluding old substrate
            atoms = []
            for i in range(refnat):
                if i in path:
                    continue
                atom = Atom(symbol=refat[i], position=refxyz[i, :])
                atoms.append(atom)

            # atoms from new substrate

            # Should we rotate the substrate around covalent bond?
            if rotate:
                theta = rotate
                origin = refxyz[center, :]
                partner = outxyz[0, :]
                outxyz2 = np.zeros(shape=outxyz.shape, dtype=np.float64)
                # reference shift
                refShift = getRodriguezRotation(partner, origin, partner, theta)
                shift = outxyz[0, :] - refShift
                for i in range(newnat):
                    outxyz2[i, :] = getRodriguezRotation(
                        outxyz[i, :], origin, partner, theta,
                    )
                    outxyz2[i, :] += shift

                    atom = Atom(symbol=newat[i], position=outxyz2[i, :])
                    atoms.append(atom)
            else:
                for i in range(newnat):
                    atom = Atom(symbol=newat[i], position=outxyz[i, :])
                    atoms.append(atom)

            # create molecule from atoms
            mol = Molecule(symbols=atoms)

            # write new structure in xyz format
            mol.writeMolecule(name + ".xyz")

            # write constrain file
            refnat = ref.get_number_of_atoms()
            oldnat = oldsub.get_number_of_atoms()
            nat = refnat - oldnat
            writeTransitionMetalConstrains(nat, newnat, newSubBonds, exclude)


def getRodriguezRotation(
    point: np.ndarray, origin: np.ndarray, partner: np.ndarray, theta: float
):
    """Let 'unit' be a unit vector defining a rotation axis, and let 'point' be
    any vector to rotate about 'unit' by angle 'theta' (right hand rule).
    Using the dot and cross products, the vector 'point' can be decomposed into
    components parallel (vector projection) and perpendicular
    (vector rejection) to the axis 'unit'."""

    # Define vector of covalent bond and normalize
    unit = origin + partner
    unit /= np.linalg.norm(unit)

    # Cross product between vector and point
    vecCrossPoint = np.cross(unit, point)

    # Dot product between vector and point
    vecDotPoint = np.dot(unit, point)

    # Transform degrees to radians
    theta = np.radians(theta)

    vnew = np.zeros(shape=(3,), dtype=np.float64)
    for i in range(3):
        vnew[i] = point[i] * np.cos(theta)
        vnew[i] += vecCrossPoint[i] * np.sin(theta)
        vnew[i] += unit[i] * vecDotPoint * (1 - np.cos(theta))

    return vnew


def getRotationMatrix(
    origin: np.ndarray, partner: np.ndarray, theta: float
) -> np.ndarray:
    """Calculate a rotation matrix for rotation along
    vector pointing from partner -> origin."""

    # get vector and normalize
    n = np.zeros(shape=(3,), dtype=np.float64)
    n = partner - origin
    n /= np.linalg.norm(n)

    # Transform degrees to radians
    theta = np.radians(theta)

    # setup rotation matrix
    omcost = 1 - np.cos(theta)
    sint = np.sin(theta)
    cost = np.cos(theta)

    u = np.zeros(shape=(3, 3), dtype=np.float64)
    u[0][0] = cost + np.power(n[0], 2) * omcost
    u[0][1] = n[0] * n[1] * omcost - n[2] * sint
    u[0][2] = n[0] * n[2] * omcost + n[1] * sint

    u[1][0] = n[1] * n[0] * omcost + n[2] * sint
    u[1][1] = cost + np.power(n[1], 2) * omcost
    u[1][2] = n[1] * n[2] * omcost - n[2] * sint

    u[2][0] = n[2] * n[0] * omcost - n[1] * sint
    u[2][1] = n[2] * n[1] * omcost + n[0] * sint
    u[2][2] = cost + np.power(n[2], 2) * omcost

    return u


def getSubstructureFromPath(refmol: Molecule, path: np.ndarray) -> Molecule:
    """Create a molecules object from a substructure of given molecule."""

    # initialize
    refat = refmol.get_atomic_numbers()
    refxyz = refmol.get_positions()

    # atoms from complex excluding old substrate
    atoms = []
    for elem in path:
        atom = Atom(symbol=refat[elem], position=refxyz[elem, :])
        atoms.append(atom)

    # create molecule from atoms
    return Molecule(symbols=atoms)


def writeTransitionMetalConstrains(
    shift: int, n: int, bonds: np.ndarray, exclude: bool
):
    """Write out constrain file for GFN-xTB."""

    import os

    f = open("constrain.inp", "w")
    s = os.linesep
    shiftExclude = 1
    if exclude:
        shiftExclude = 0
    f.write("$fix" + s)
    f.write(" atoms: 1-{}".format(shift + shiftExclude) + s)
    f.write("$constrain" + s)
    for i in range(n):
        for partner in bonds[i]:
            f.write(
                " distance: {}, {}, auto".format(i + 1 + shift, partner + 1 + shift) + s
            )
    f.write("$end" + s)
    f.close()


def matchSubstrates(
    bonds: np.ndarray,
    new: Molecule,
    bondsNewSub: np.ndarray,
    old: Molecule,
    bondsOldSub,
    center: np.ndarray,
) -> np.ndarray:
    """Match substrates by root mean squared deviation measure."""

    # check is central atom given
    centered = False
    count = 1
    if center is not None:
        centered = True
        count = 2

    newnat = new.get_number_of_atoms()
    newxyz = np.zeros(shape=(newnat, 3), dtype=np.float64)
    newxyz = new.get_positions()

    oldnat = old.get_number_of_atoms()
    oldxyz = np.zeros(shape=(oldnat, 3), dtype=np.float64)
    oldxyz = old.get_positions()

    # shift to first atom of old substrate
    shift = np.zeros(shape=(3,), dtype=np.float64)
    shift = oldxyz[0, :]
    oldShift = shift
    oldxyz2 = np.zeros(shape=(oldnat, 3), dtype=np.float64)
    for i in range(oldnat):
        oldxyz2[i][:] = oldxyz[i][:] - shift

    # shift to first atom of new substrate
    shift = newxyz[0, :]
    newxyz2 = np.zeros(shape=(newnat, 3), dtype=np.float64)
    for i in range(newnat):
        newxyz2[i][:] = newxyz[i][:] - shift

    # covalent bonding partner in old substrate
    covOld = len(bondsOldSub[0][:])

    # covalent bonding partner in new substrate
    covNew = len(bondsNewSub[0][:])

    # Check for linear moleule case
    if covNew != 1:
        isNotLinear = True
    else:
        isNotLinear = False

    covOld += count
    covNew += count

    # get structures for RMSD alignment
    shiftOldSub = np.zeros(shape=(covOld, 3), dtype=np.float64)
    shiftOldSub[0, :] = oldxyz2[0, :]
    shiftNewSub = np.zeros(shape=(covNew, 3), dtype=np.float64)
    shiftNewSub[0, :] = newxyz2[0, :]

    # shift old substrate
    i = 0
    for j in range(len(bondsOldSub[0][:])):
        i += 1
        k = bondsOldSub[0][j]
        shiftOldSub[i, :] = oldxyz2[k, :]

    # shift new substrate
    i = 0
    for j in range(len(bondsNewSub[0][:])):
        i += 1
        k = bondsNewSub[0][j]
        shiftNewSub[i, :] = newxyz2[k, :]

    if centered:
        indexOld = covOld - 1
        indexNew = covNew - 1
        shiftOldSub[indexOld, :] = center - oldShift
        distRef = getDist(shiftOldSub[0, :], shiftOldSub[indexOld, :])
        shiftNewSub[indexNew][:] = 0
        shift = getNewSubstrateCenter(indexNew, shiftNewSub, distRef)

    bdim = np.minimum(covOld, covNew)
    print(bdim)

    if bdim >= 3:
        bdim = 3
        b1xyz = np.zeros(shape=(bdim, 3), dtype=np.ndarray)
        b2xyz = np.zeros(shape=(bdim, 3), dtype=np.ndarray)
        indexOld = covOld - 1
        indexNew = covNew - 1
        b1xyz[0][:] = shiftOldSub[0][:]
        b1xyz[1][:] = shiftOldSub[1][:]
        b1xyz[2][:] = shiftOldSub[indexOld][:]
        b2xyz[0][:] = shiftNewSub[0][:]
        b2xyz[1][:] = shiftNewSub[1][:]
        b2xyz[2][:] = shiftNewSub[indexNew][:]
    else:
        bdim = 2
        b1xyz = np.zeros(shape=(bdim, 3), dtype=np.ndarray)
        b2xyz = np.zeros(shape=(bdim, 3), dtype=np.ndarray)
        indexOld = covOld - 1
        indexNew = covNew - 1
        b1xyz[0][:] = shiftOldSub[0][:]
        b1xyz[1][:] = shiftOldSub[indexOld][:]
        b2xyz[0][:] = shiftNewSub[0][:]
        b2xyz[1][:] = shiftNewSub[indexNew][:]

    tmpxyz = np.zeros(shape=(newnat, 3), dtype=np.float64)
    outxyz = np.zeros(shape=(newnat, 3), dtype=np.float64)

    if isNotLinear:
        u = np.zeros(shape=(3, 3), dtype=np.float64)
        error = 0.0
        # get RMSD value and rotation matrix
        error, u = rmsd(bdim, b2xyz, b1xyz)
        tmpxyz = np.matmul(u.T, newxyz2[:][0:newnat].T)
        # shift
        for i in range(newnat):
            outxyz[i, :] = tmpxyz.T[i, :] + oldShift

        return outxyz

    tmpxyz = newxyz2[:][0:newnat]
    # shift
    for i in range(newnat):
        outxyz[i, :] = tmpxyz[i, :] + oldShift

    return outxyz


def getNewSubstrateCenter(index: int, shift: np.ndarray, distRef: float):
    """Get the position of the new substrate."""

    n = index - 2
    center = np.zeros(shape=(3,), dtype=np.ndarray)

    if n > 0:
        center = np.sum(shift, axis=0)
        distCenter = getDist(shift[0][:], center)
        distRel = distRef / distCenter
        center *= -distRel
        shift[index][:] = center
    else:
        # atom case: place ligand along the axis
        shift[index][n] = distRef
        shift[index][2:] = 0.0

    return shift[index][:]


def getDist(a: np.ndarray, b: np.ndarray) -> float:
    """Calculate euclidean distance between a and b."""

    dist = 0.0

    # sanity check for dimensions
    if len(a) == len(b):
        for i in range(3):
            dist += np.power(a[i] - b[i], 2)

        dist = np.sqrt(dist)
    else:
        raise RuntimeError

    return dist
