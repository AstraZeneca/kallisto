# src/kallisto/rmsd.py

from collections import Counter
from typing import Tuple

import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.spatial import distance
from scipy.spatial.transform import Rotation as R

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
    # use fac = -1 instead of u.T below
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
    eigenval, eigenvec = eigsh(smat, 1, which="LA")
    eigenvec = eigenvec.squeeze()

    # convert quaternion eigenvec to rotation matrix U
    u = rotationMatrix(eigenvec)

    # root mean squared deviation
    error = np.sqrt(np.maximum(0.0, ((x_norm + y_norm) - 2 * eigenval) / float(n)))

    return error, u


def rotationMatrix(q: np.ndarray) -> np.ndarray:
    """Constructs rotation matrix U from quaternion q."""

    q = q.squeeze()
    u = R.from_quat(q).as_matrix()

    # anti-transpose, transpose, and fix signs
    u = u[::-1, ::-1].T
    u[0, 0] *= -1
    u[1, 1] *= -1
    u[1, 2] *= -1
    u[2, 0] *= -1

    return u


def recursiveGetSubstructures(n: int, bonds, center: int):
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
    bonds,
    ref: Molecule,
    newsub: Molecule,
    newSubBonds,
    name: str,
    rotate: int,
    exclude: bool,
) -> Molecule:
    """Exchange substructure (subnr) from ref with substrate."""

    # setup initial complex
    refxyz = ref.get_positions()
    refat = ref.get_atomic_numbers()
    refnat = ref.get_number_of_atoms()

    # match new and old substrate
    newat = newsub.get_atomic_numbers()
    newnat = newsub.get_number_of_atoms()

    # extract coordinates of central atom
    centralAtom = refxyz[center, :]

    mol = Molecule()
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
            outxyz = matchSubstrates(
                bonds,
                newsub,
                newSubBonds,
                oldsub,
                oldSubBonds,
                centralAtom,
            )

            # atoms from complex excluding old substrate
            atoms = []
            for j in range(refnat):
                if j in path:
                    continue
                atom = Atom(symbol=refat[j], position=refxyz[j, :])
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
                for j in range(newnat):
                    outxyz2[j, :] = getRodriguezRotation(
                        outxyz[j, :],
                        origin,
                        partner,
                        theta,
                    )
                    outxyz2[j, :] += shift

                    atom = Atom(symbol=newat[j], position=outxyz2[j, :])
                    atoms.append(atom)
            else:
                for j in range(newnat):
                    atom = Atom(symbol=newat[j], position=outxyz[j, :])
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
            return mol

    return mol


def getRodriguezRotation(
    point: np.ndarray, origin: np.ndarray, partner: np.ndarray, theta: float
):
    """Let 'unit' be a unit vector defining a rotation axis, and let 'point' be
    any vector to rotate about 'unit' by angle 'theta' (right hand rule).
    Using the dot and cross products, the vector 'point' can be decomposed into
    components parallel (vector projection) and perpendicular
    (vector rejection) to the axis 'unit'."""

    # Define vector of covalent bond and normalize
    unit = origin - partner
    unit /= np.linalg.norm(unit)

    # Transform degrees to radians
    theta = np.radians(theta)

    # define rotation and apply
    r = R.from_rotvec(theta * unit)
    return r.apply(point)


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
        shiftNewSub[indexNew][:] = 0
        distRef = distance.euclidean(shiftOldSub[0, :], shiftOldSub[indexOld, :])
        getNewSubstrateCenter(indexNew, shiftNewSub, distRef)

    bdim = np.minimum(covOld, covNew)

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
        # get RMSD value and rotation matrix
        _, u = rmsd(bdim, b2xyz, b1xyz)
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

    if n > 0:
        center = np.sum(shift, axis=0)
        distCenter = distance.euclidean(shift[0][:], center)
        distRel = distRef / distCenter
        center *= -distRel
        shift[index][:] = center
    else:
        # atom case: place ligand along the axis
        shift[index][n] = distRef
        shift[index][2:] = 0.0

    return shift[index][:]
