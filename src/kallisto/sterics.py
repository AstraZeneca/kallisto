# src/kallisto/sterics.py

import numpy as np
from scipy.spatial.transform import Rotation as R  # type: ignore

from kallisto.molecule import Molecule


def getClassicalSterimol(mol: Molecule, origin: int, partner: int):
    """Verloop definitions for L, B1, and B5 steric parameter.
    We apply kallisto van der Waals radii."""

    # initialize coordinates
    coords = mol.get_positions()

    # get number of atoms
    nat = mol.get_number_of_atoms()

    # get van der Waals radii in Bohr
    vdw = np.zeros(shape=(nat,), dtype=np.float64)
    vdw = mol.get_vdw(charge=0, vdwtype="rahm", scale=1)

    # shift all atoms wrt origin
    coords -= coords[origin]

    # extract vector origin -> attachted and normalize
    vector = np.zeros(shape=(3,), dtype=np.float64)
    vector = coords[partner] - coords[origin]
    vector /= np.linalg.norm(vector)

    # align vector and x-axis
    xaxis = np.array([1, 0, 0])
    dot = np.dot(vector, xaxis).reshape(1) + 1

    # define zaxis and cover antiparallel case
    zaxis = np.array([0, 0, 1])
    epsilon = 1e-06
    if dot < epsilon:
        w = np.cross(vector, zaxis)
        if np.linalg.norm(w) < epsilon:
            w = np.cross(vector, xaxis)
    else:
        w = np.cross(vector, xaxis)

    # define quaternion q and normalize
    q = np.concatenate((w, dot))
    q /= np.linalg.norm(q)

    # rotate coordinates according to q
    u = R.from_quat(q)
    coords = u.apply(coords)

    # project coordinates on vector
    vector = coords[partner] - coords[origin]
    unitVector = vector / np.linalg.norm(vector)

    # project coordinates on vectors
    vdw = np.vstack(vdw).reshape(-1)
    cvalues = np.dot(unitVector.reshape(1, -1), coords.T)
    projected = cvalues + vdw

    lval = np.max(projected)
    L = unitVector * lval
    L = L.reshape(-1)

    # B min and max values
    r = 1
    slices = 360
    theta = np.linspace(0, 2 * np.pi, slices)
    x = np.zeros(shape=(len(theta),), dtype=np.float64)
    y = r * np.cos(theta)
    z = r * np.sin(theta)
    vectors = np.column_stack((x, y, z))

    cvalues = np.dot(vectors, coords.T)
    projected = cvalues + vdw

    maxValues = np.max(projected, axis=1)
    bminVal = np.min(maxValues)
    bmaxVal = np.max(maxValues)

    return lval, bminVal, bmaxVal
