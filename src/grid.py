# src/grid.py

from collections import Counter

import numpy as np

# Angular Lebedev-Laikov grids (6-48 equivalent points)
#
# Original (C): Dmitri N. Laikov
# Adapted (Fortran): Christoph van Wuellen (CvW)
# Translated into Python by permission of CvW
#
# Reference:
# V.I. Lebedev, and D.N. Laikov, "A quadrature formula for the sphere of the
# 131st algebraic order of accuracy", Doklady Mathematics, Vol. 59, No. 3,
# 1999, pp. 477-481.

# Define grid sizes
gridSize = dict(
    [
        (0, 6),
        (1, 14),
        (2, 26),
        (3, 38),
        (4, 50),
        (5, 74),
        (6, 86),
        (7, 110),
        (8, 146),
        (9, 170),
        (10, 194),
    ]
)

# keep octahedral symmetry
z = 0.0


def getLebedevLaikovGrid(n: int):
    """Generate angular Lebedev-Laikov grid.

    Parameter:
     n - size of grid
     grid - grid points
     weights - point weights"""

    # initialize count
    count = Counter()  # type: ignore
    count["points"] = 0

    size = gridSize[n]
    grid = np.zeros(shape=(size, 3), dtype=np.float64)
    weights = np.zeros(shape=(size,), dtype=np.float64)

    if n == 0:
        grid, weights, count = ld0006(grid, weights, count)
    elif n == 1:
        grid, weights, count = ld0014(grid, weights, count)
    elif n == 2:
        grid, weights, count = ld0026(grid, weights, count)
    elif n == 3:
        grid, weights, count = ld0038(grid, weights, count)
    elif n == 4:
        grid, weights, count = ld0050(grid, weights, count)
    elif n == 5:
        grid, weights, count = ld0074(grid, weights, count)
    elif n == 6:
        grid, weights, count = ld0086(grid, weights, count)
    elif n == 7:
        grid, weights, count = ld0110(grid, weights, count)
    elif n == 8:
        grid, weights, count = ld0146(grid, weights, count)
    elif n == 9:
        grid, weights, count = ld0170(grid, weights, count)
    else:
        grid, weights, count = ld0194(grid, weights, count)

    if count["points"] != gridSize[n]:
        raise RuntimeError(
            "Number of grid points is not matching the referencecount. Stop here."
        )

    return grid, weights


# Given a point on a sphere (specified by a and b), generate all equivalent
# points with octahedral symmetry, making grid points each having a weight v.


def genOh1(count: Counter, coord: np.ndarray, weights: np.ndarray, v: float):
    """6 points with octahedral symmetry from (0,0,1)"""

    a = 1.0
    i = count["points"]
    coord[i + 0, :] = [a, z, z]
    coord[i + 1, :] = [-a, z, z]
    coord[i + 2, :] = [z, a, z]
    coord[i + 3, :] = [z, -a, z]
    coord[i + 4, :] = [z, z, a]
    coord[i + 5, :] = [z, z, -a]
    weights[i : i + 6] = v
    count["points"] += 6


def genOh2(count: Counter, coord: np.ndarray, weights: np.ndarray, v: float):
    """12 points from (0,a,a) with a=1/sqrt(2)."""

    a = np.sqrt(0.5)
    i = count["points"]
    coord[i + 0, :] = [z, a, a]
    coord[i + 1, :] = [z, -a, a]
    coord[i + 2, :] = [z, a, -a]
    coord[i + 3, :] = [z, -a, -a]
    coord[i + 4, :] = [a, z, a]
    coord[i + 5, :] = [-a, z, a]
    coord[i + 6, :] = [a, z, -a]
    coord[i + 7, :] = [-a, z, -a]
    coord[i + 8, :] = [a, a, z]
    coord[i + 9, :] = [-a, a, z]
    coord[i + 10, :] = [a, -a, z]
    coord[i + 11, :] = [-a, -a, z]
    weights[i : i + 12] = v
    count["points"] += 12


def genOh3(count: Counter, coord: np.ndarray, weights: np.ndarray, v: float):
    """8 points from (a,a,a), a=1/sqrt(3)"""

    a = np.sqrt(1.0 / 3.0)
    i = count["points"]
    coord[i + 0, :] = [a, a, a]
    coord[i + 1, :] = [-a, a, a]
    coord[i + 2, :] = [a, -a, a]
    coord[i + 3, :] = [-a, -a, a]
    coord[i + 4, :] = [a, a, -a]
    coord[i + 5, :] = [-a, a, -a]
    coord[i + 6, :] = [a, -a, -a]
    coord[i + 7, :] = [-a, -a, -a]
    weights[i : i + 8] = v
    count["points"] += 8


def genOh4(count: Counter, coord: np.ndarray, weights: np.ndarray, a: float, v: float):
    """24 points from (a,a,b), b=sqrt(1-2 a^2)"""

    b = np.sqrt(1.0 - 2.0 * a * a)
    i = count["points"]
    coord[i + 0, :] = [a, a, b]
    coord[i + 1, :] = [-a, a, b]
    coord[i + 2, :] = [a, -a, b]
    coord[i + 3, :] = [-a, -a, b]
    coord[i + 4, :] = [a, a, -b]
    coord[i + 5, :] = [-a, a, -b]
    coord[i + 6, :] = [a, -a, -b]
    coord[i + 7, :] = [-a, -a, -b]
    coord[i + 8, :] = [a, b, a]
    coord[i + 9, :] = [-a, b, a]
    coord[i + 10, :] = [a, -b, a]
    coord[i + 11, :] = [-a, -b, a]
    coord[i + 12, :] = [a, b, -a]
    coord[i + 13, :] = [-a, b, -a]
    coord[i + 14, :] = [a, -b, -a]
    coord[i + 15, :] = [-a, -b, -a]
    coord[i + 16, :] = [b, a, a]
    coord[i + 17, :] = [-b, a, a]
    coord[i + 18, :] = [b, -a, a]
    coord[i + 19, :] = [-b, -a, a]
    coord[i + 20, :] = [b, a, -a]
    coord[i + 21, :] = [-b, a, -a]
    coord[i + 22, :] = [b, -a, -a]
    coord[i + 23, :] = [-b, -a, -a]
    weights[i : i + 24] = v
    count["points"] += 24


def genOh5(count: Counter, coord: np.ndarray, weights: np.ndarray, a: float, v: float):
    """24 points from (a,b,0), b=sqrt(1-a^2)."""

    i = count["points"]
    b = np.sqrt(1.0 - a * a)
    coord[i + 0, :] = [a, b, z]
    coord[i + 1, :] = [-a, b, z]
    coord[i + 2, :] = [a, -b, z]
    coord[i + 3, :] = [-a, -b, z]
    coord[i + 4, :] = [b, a, z]
    coord[i + 5, :] = [-b, a, z]
    coord[i + 6, :] = [b, -a, z]
    coord[i + 7, :] = [-b, -a, z]
    coord[i + 8, :] = [a, z, b]
    coord[i + 9, :] = [-a, z, b]
    coord[i + 10, :] = [a, z, -b]
    coord[i + 11, :] = [-a, z, -b]
    coord[i + 12, :] = [b, z, a]
    coord[i + 13, :] = [-b, z, a]
    coord[i + 14, :] = [b, z, -a]
    coord[i + 15, :] = [-b, z, -a]
    coord[i + 16, :] = [z, a, b]
    coord[i + 17, :] = [z, -a, b]
    coord[i + 18, :] = [z, a, -b]
    coord[i + 19, :] = [z, -a, -b]
    coord[i + 20, :] = [z, b, a]
    coord[i + 21, :] = [z, -b, a]
    coord[i + 22, :] = [z, b, -a]
    coord[i + 23, :] = [z, -b, -a]
    weights[i : i + 24] = v
    count["points"] += 24


def genOh6(
    count: Counter, coord: np.ndarray, weights: np.ndarray, a: float, b: float, v: float
):
    """48 points from (a,b,c), c=sqrt(1-a^2-b^2)"""

    i = count["points"]
    c = np.sqrt(1.0 - a * a - b * b)
    coord[i + 0, :] = [a, b, c]
    coord[i + 1, :] = [-a, b, c]
    coord[i + 2, :] = [a, -b, c]
    coord[i + 3, :] = [-a, -b, c]
    coord[i + 4, :] = [a, b, -c]
    coord[i + 5, :] = [-a, b, -c]
    coord[i + 6, :] = [a, -b, -c]
    coord[i + 7, :] = [-a, -b, -c]
    coord[i + 8, :] = [a, c, b]
    coord[i + 9, :] = [-a, c, b]
    coord[i + 10, :] = [a, -c, b]
    coord[i + 11, :] = [-a, -c, b]
    coord[i + 12, :] = [a, c, -b]
    coord[i + 13, :] = [-a, c, -b]
    coord[i + 14, :] = [a, -c, -b]
    coord[i + 15, :] = [-a, -c, -b]
    coord[i + 16, :] = [b, a, c]
    coord[i + 17, :] = [-b, a, c]
    coord[i + 18, :] = [b, -a, c]
    coord[i + 19, :] = [-b, -a, c]
    coord[i + 20, :] = [b, a, -c]
    coord[i + 21, :] = [-b, a, -c]
    coord[i + 22, :] = [b, -a, -c]
    coord[i + 23, :] = [-b, -a, -c]
    coord[i + 24, :] = [b, c, a]
    coord[i + 25, :] = [-b, c, a]
    coord[i + 26, :] = [b, -c, a]
    coord[i + 27, :] = [-b, -c, a]
    coord[i + 28, :] = [b, c, -a]
    coord[i + 29, :] = [-b, c, -a]
    coord[i + 30, :] = [b, -c, -a]
    coord[i + 31, :] = [-b, -c, -a]
    coord[i + 32, :] = [c, a, b]
    coord[i + 33, :] = [-c, a, b]
    coord[i + 34, :] = [c, -a, b]
    coord[i + 35, :] = [-c, -a, b]
    coord[i + 36, :] = [c, a, -b]
    coord[i + 37, :] = [-c, a, -b]
    coord[i + 38, :] = [c, -a, -b]
    coord[i + 39, :] = [-c, -a, -b]
    coord[i + 40, :] = [c, b, a]
    coord[i + 41, :] = [-c, b, a]
    coord[i + 42, :] = [c, -b, a]
    coord[i + 43, :] = [-c, -b, a]
    coord[i + 44, :] = [c, b, -a]
    coord[i + 45, :] = [-c, b, -a]
    coord[i + 46, :] = [c, -b, -a]
    coord[i + 47, :] = [-c, -b, -a]
    weights[i : i + 48] = v
    count["points"] += 48


def ld0006(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 6 point angular Lebedev-Laikov grid."""

    v = 0.1666666666666667
    genOh1(count, coord, weights, v)
    return coord, weights, count


def ld0014(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 14 point angular Lebedev-Laikov grid."""

    v = 0.6666666666666667
    genOh1(count, coord, weights, v)
    v = 0.7500000000000000
    genOh3(count, coord, weights, v)
    return coord, weights, count


def ld0026(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 26 point angular Lebedev-Laikov grid."""

    v = 0.4761904761904762e-1
    genOh1(count, coord, weights, v)
    v = 0.3809523809523810e-1
    genOh2(count, coord, weights, v)
    v = 0.3214285714285714e-1
    genOh3(count, coord, weights, v)
    return coord, weights, count


def ld0038(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 38 point angular Lebedev-Laikov grid."""

    v = 0.9523809523809524e-2
    genOh1(count, coord, weights, v)
    v = 0.3214285714285714e-1
    genOh3(count, coord, weights, v)
    a = 0.4597008433809831
    v = 0.2857142857142857e-1
    genOh5(count, coord, weights, a, v)
    return coord, weights, count


def ld0050(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 50 point angular Lebedev-Laikov grid."""

    v = 0.1269841269841270e-1
    genOh1(count, coord, weights, v)
    v = 0.2257495590828924e-1
    genOh2(count, coord, weights, v)
    v = 0.2109375000000000e-1
    genOh3(count, coord, weights, v)
    a = 0.3015113445777636
    v = 0.2017333553791887
    genOh4(count, coord, weights, a, v)
    return coord, weights, count


def ld0074(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 74 point angular Lebedev-Laikov grid."""

    v = 0.5130671797338464e-3
    genOh1(count, coord, weights, v)
    v = 0.1660406956574204e-1
    genOh2(count, coord, weights, v)
    v = -0.2958603896103896e-1
    genOh3(count, coord, weights, v)
    a = 0.4803844614152614
    v = 0.2657620708215946e-1
    genOh4(count, coord, weights, a, v)
    a = 0.3207726489807764
    v = 0.1652217099371571e-1
    genOh5(count, coord, weights, a, v)
    return coord, weights, count


def ld0086(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 86 point angular Lebedev-Laikov grid."""

    v = 0.1154401154401154e-1
    genOh1(count, coord, weights, v)
    v = 0.1194390908585628e-1
    genOh3(count, coord, weights, v)
    a = 0.3696028464541502
    v = 0.1111055571060340e-1
    genOh4(count, coord, weights, a, v)
    a = 0.6943540066026664
    v = 0.1187650129453714e-1
    genOh4(count, coord, weights, a, v)
    a = 0.3742430390903412
    v = 0.1181230374690448e-1
    genOh5(count, coord, weights, a, v)
    return coord, weights, count


def ld0110(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 110 point angular Lebedev-Laikov grid."""

    v = 0.3828270494937162e-2
    genOh1(count, coord, weights, v)
    v = 0.9793737512487512e-2
    genOh3(count, coord, weights, v)
    a = 0.1851156353447362
    v = 0.8211737283191111e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6904210483822922
    v = 0.9942814891178103e-2
    genOh4(count, coord, weights, a, v)
    a = 0.3956894730559419
    v = 0.9595471336070963e-2
    genOh4(count, coord, weights, a, v)
    a = 0.4783690288121502
    v = 0.9694996361663028e-2
    genOh5(count, coord, weights, a, v)
    return coord, weights, count


def ld0146(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 146 point angular Lebedev-Laikov grid."""

    v = 0.5996313688621381e-3
    genOh1(count, coord, weights, v)
    v = 0.7372999718620756e-2
    genOh2(count, coord, weights, v)
    v = 0.7210515360144488e-2
    genOh3(count, coord, weights, v)
    a = 0.6764410400114264
    v = 0.7116355493117555e-2
    genOh4(count, coord, weights, a, v)
    a = 0.4174961227965453
    v = 0.6753829486314477e-2
    genOh4(count, coord, weights, a, v)
    a = 0.1574676672039082
    v = 0.7574394159054034e-2
    genOh4(count, coord, weights, a, v)
    a = 0.1403553811713183
    b = 0.4493328323269557
    v = 0.6991087353303262e-2
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count


def ld0170(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 170 point angular Lebedev-Laikov grid."""

    v = 0.5544842902037365e-2
    genOh1(count, coord, weights, v)
    v = 0.6071332770670752e-2
    genOh2(count, coord, weights, v)
    v = 0.6383674773515093e-2
    genOh3(count, coord, weights, v)
    a = 0.2551252621114134
    v = 0.5183387587747790e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6743601460362766
    v = 0.6317929009813725e-2
    genOh4(count, coord, weights, a, v)
    a = 0.4318910696719410
    v = 0.6201670006589077e-2
    genOh4(count, coord, weights, a, v)
    a = 0.2613931360335988
    v = 0.5477143385137348e-2
    genOh5(count, coord, weights, a, v)
    a = 0.4990453161796037
    b = 0.1446630744325115
    v = 0.5968383987681156e-2
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count


def ld0194(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 194 point angular Lebedev-Laikov grid."""

    v = 0.1782340447244611e-2
    genOh1(count, coord, weights, v)
    v = 0.5716905949977102e-2
    genOh2(count, coord, weights, v)
    v = 0.5573383178848738e-2
    genOh3(count, coord, weights, v)
    a = 0.6712973442695226
    v = 0.5608704082587997e-2
    genOh4(count, coord, weights, a, v)
    a = 0.2892465627575439
    v = 0.5158237711805383e-2
    genOh4(count, coord, weights, a, v)
    a = 0.4446933178717437
    v = 0.5518771467273614e-2
    genOh4(count, coord, weights, a, v)
    a = 0.1299335447650067
    v = 0.4106777028169394e-2
    genOh4(count, coord, weights, a, v)
    a = 0.3457702197611283
    v = 0.5051846064614808e-2
    genOh5(count, coord, weights, a, v)
    a = 0.1590417105383530
    b = 0.8360360154824589
    v = 0.5530248916233094e-2
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count
