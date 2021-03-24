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
        (11, 230),
        (12, 266),
        (13, 302),
        (14, 350),
        (15, 434),
        (16, 590),
        (17, 770),
        (18, 974),
        (19, 1202),
        (20, 1454),
        (21, 1730),
        (22, 2030),
    ]
)

# keep octahedral symmetry
z = 0.0


def getGridFromInt(n: int):
    """Returns Lebedev-Laikov grid, weights, and count depending on input.

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

    switcher = {
        0: lambda: ld0006(grid, weights, count),
        1: lambda: ld0014(grid, weights, count),
        2: lambda: ld0026(grid, weights, count),
        3: lambda: ld0038(grid, weights, count),
        4: lambda: ld0050(grid, weights, count),
        5: lambda: ld0074(grid, weights, count),
        6: lambda: ld0086(grid, weights, count),
        7: lambda: ld0110(grid, weights, count),
        8: lambda: ld0146(grid, weights, count),
        9: lambda: ld0170(grid, weights, count),
        10: lambda: ld0194(grid, weights, count),
        11: lambda: ld0230(grid, weights, count),
        12: lambda: ld0266(grid, weights, count),
        13: lambda: ld0302(grid, weights, count),
        14: lambda: ld0350(grid, weights, count),
        15: lambda: ld0434(grid, weights, count),
        16: lambda: ld0590(grid, weights, count),
        17: lambda: ld0770(grid, weights, count),
        18: lambda: ld0974(grid, weights, count),
        19: lambda: ld1202(grid, weights, count),
        20: lambda: ld1454(grid, weights, count),
        21: lambda: ld1730(grid, weights, count),
        22: lambda: ld2030(grid, weights, count),
    }  # type: ignore
    func = switcher.get(n, lambda: "Invalid Lebedev-Laikov grid requested.")
    return func()


def getLebedevLaikovGrid(n: int):
    """Generate angular Lebedev-Laikov grid.

    Parameter:
     n - size of grid"""

    # input integer and get Lebedev-Laikov grid, weights, and count
    grid, weights, count = getGridFromInt(n)

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


def ld0230(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 230 point angular Lebedev-Laikov grid."""

    v = -0.5522639919727325e-1
    genOh1(count, coord, weights, v)
    v = 0.4450274607445226e-2
    genOh3(count, coord, weights, v)
    a = 0.4492044687397611e0
    v = 0.4496841067921404e-2
    genOh4(count, coord, weights, a, v)
    a = 0.2520419490210201e0
    v = 0.5049153450478750e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6981906658447242e0
    v = 0.3976408018051883e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6587405243460960e0
    v = 0.4401400650381014e-2
    genOh4(count, coord, weights, a, v)
    a = 0.4038544050097660e-1
    v = 0.1724544350544401e-1
    genOh4(count, coord, weights, a, v)
    a = 0.5823842309715585e0
    v = 0.4231083095357343e-2
    genOh5(count, coord, weights, a, v)
    a = 0.3545877390518688e0
    v = 0.5198069864064399e-2
    genOh5(count, coord, weights, a, v)
    a = 0.2272181808998187e0
    b = 0.4864661535886647e0
    v = 0.4695720972568883e-2
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count


def ld0266(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 266 point angular Lebedev-Laikov grid."""

    v = -0.1313769127326952e-2
    genOh1(count, coord, weights, v)
    v = -0.2522728704859336e-2
    genOh2(count, coord, weights, v)
    v = 0.4186853881700583e-2
    genOh3(count, coord, weights, v)
    a = 0.7039373391585475e0
    v = 0.5315167977810885e-2
    genOh4(count, coord, weights, a, v)
    a = 0.1012526248572414e0
    v = 0.4047142377086219e-2
    genOh4(count, coord, weights, a, v)
    a = 0.4647448726420539e0
    v = 0.4112482394406990e-2
    genOh4(count, coord, weights, a, v)
    a = 0.3277420654971629e0
    v = 0.3595584899758782e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6620338663699974e0
    v = 0.4256131351428158e-2
    genOh4(count, coord, weights, a, v)
    a = 0.8506508083520399e0
    v = 0.4229582700647240e-2
    genOh5(count, coord, weights, a, v)
    a = 0.3233484542692899e0
    b = 0.1153112011009701e0
    v = 0.4080914225780505e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.2314790158712601e0
    b = 0.5244939240922365e0
    v = 0.4071467593830964e-2
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count


def ld0302(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 302 point angular Lebedev-Laikov grid."""

    v = 0.8545911725128148e-3
    genOh1(count, coord, weights, v)
    v = 0.3599119285025571e-2
    genOh3(count, coord, weights, v)
    a = 0.3515640345570105e0
    v = 0.3449788424305883e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6566329410219612e0
    v = 0.3604822601419882e-2
    genOh4(count, coord, weights, a, v)
    a = 0.4729054132581005e0
    v = 0.3576729661743367e-2
    genOh4(count, coord, weights, a, v)
    a = 0.9618308522614784e-1
    v = 0.2352101413689164e-2
    genOh4(count, coord, weights, a, v)
    a = 0.2219645236294178e0
    v = 0.3108953122413675e-2
    genOh4(count, coord, weights, a, v)
    a = 0.7011766416089545e0
    v = 0.3650045807677255e-2
    genOh4(count, coord, weights, a, v)
    a = 0.2644152887060663e0
    v = 0.2982344963171804e-2
    genOh5(count, coord, weights, a, v)
    a = 0.5718955891878961e0
    v = 0.3600820932216460e-2
    genOh5(count, coord, weights, a, v)
    a = 0.2510034751770465e0
    b = 0.8000727494073952e0
    v = 0.3571540554273387e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.1233548532583327e0
    b = 0.4127724083168531e0
    v = 0.3392312205006170e-2
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count


def ld0350(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 350 point angular Lebedev-Laikov grid."""

    v = 0.3006796749453936e-2
    genOh1(count, coord, weights, v)
    v = 0.3050627745650771e-2
    genOh3(count, coord, weights, v)
    a = 0.7068965463912316e0
    v = 0.1621104600288991e-2
    genOh4(count, coord, weights, a, v)
    a = 0.4794682625712025e0
    v = 0.3005701484901752e-2
    genOh4(count, coord, weights, a, v)
    a = 0.1927533154878019e0
    v = 0.2990992529653774e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6930357961327123e0
    v = 0.2982170644107595e-2
    genOh4(count, coord, weights, a, v)
    a = 0.3608302115520091e0
    v = 0.2721564237310992e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6498486161496169e0
    v = 0.3033513795811141e-2
    genOh4(count, coord, weights, a, v)
    a = 0.1932945013230339e0
    v = 0.3007949555218533e-2
    genOh5(count, coord, weights, a, v)
    a = 0.3800494919899303e0
    v = 0.2881964603055307e-2
    genOh5(count, coord, weights, a, v)
    a = 0.2899558825499574e0
    b = 0.7934537856582316e0
    v = 0.2958357626535696e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.9684121455103957e-1
    b = 0.8280801506686862e0
    v = 0.3036020026407088e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.1833434647041659e0
    b = 0.9074658265305127e0
    v = 0.2832187403926303e-2
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count


def ld0434(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 434 point angular Lebedev-Laikov grid."""

    v = 0.5265897968224436e-3
    genOh1(count, coord, weights, v)
    v = 0.2548219972002607e-2
    genOh2(count, coord, weights, v)
    v = 0.2512317418927307e-2
    genOh3(count, coord, weights, v)
    a = 0.6909346307509111e0
    v = 0.2530403801186355e-2
    genOh4(count, coord, weights, a, v)
    a = 0.1774836054609158e0
    v = 0.2014279020918528e-2
    genOh4(count, coord, weights, a, v)
    a = 0.4914342637784746e0
    v = 0.2501725168402936e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6456664707424256e0
    v = 0.2513267174597564e-2
    genOh4(count, coord, weights, a, v)
    a = 0.2861289010307638e0
    v = 0.2302694782227416e-2
    genOh4(count, coord, weights, a, v)
    a = 0.7568084367178018e-1
    v = 0.1462495621594614e-2
    genOh4(count, coord, weights, a, v)
    a = 0.3927259763368002e0
    v = 0.2445373437312980e-2
    genOh4(count, coord, weights, a, v)
    a = 0.8818132877794288e0
    v = 0.2417442375638981e-2
    genOh5(count, coord, weights, a, v)
    a = 0.9776428111182649e0
    v = 0.1910951282179532e-2
    genOh5(count, coord, weights, a, v)
    a = 0.2054823696403044e0
    b = 0.8689460322872412e0
    v = 0.2416930044324775e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.5905157048925271e0
    b = 0.7999278543857286e0
    v = 0.2512236854563495e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.5550152361076807e0
    b = 0.7717462626915901e0
    v = 0.2496644054553086e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.9371809858553722e0
    b = 0.3344363145343455e0
    v = 0.2236607760437849e-2
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count


def ld0590(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 590 point angular Lebedev-Laikov grid."""

    v = 0.3095121295306187e-3
    genOh1(count, coord, weights, v)
    v = 0.1852379698597489e-2
    genOh3(count, coord, weights, v)
    a = 0.7040954938227469e0
    v = 0.1871790639277744e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6807744066455243e0
    v = 0.1858812585438317e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6372546939258752e0
    v = 0.1852028828296213e-2
    genOh4(count, coord, weights, a, v)
    a = 0.5044419707800358e0
    v = 0.1846715956151242e-2
    genOh4(count, coord, weights, a, v)
    a = 0.4215761784010967e0
    v = 0.1818471778162769e-2
    genOh4(count, coord, weights, a, v)
    a = 0.3317920736472123e0
    v = 0.1749564657281154e-2
    genOh4(count, coord, weights, a, v)
    a = 0.2384736701421887e0
    v = 0.1617210647254411e-2
    genOh4(count, coord, weights, a, v)
    a = 0.1459036449157763e0
    v = 0.1384737234851692e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6095034115507196e-1
    v = 0.9764331165051050e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6116843442009876e0
    v = 0.1857161196774078e-2
    genOh5(count, coord, weights, a, v)
    a = 0.3964755348199858e0
    v = 0.1705153996395864e-2
    genOh5(count, coord, weights, a, v)
    a = 0.1724782009907724e0
    v = 0.1300321685886048e-2
    genOh5(count, coord, weights, a, v)
    a = 0.5610263808622060e0
    b = 0.3518280927733519e0
    v = 0.1842866472905286e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.4742392842551980e0
    b = 0.2634716655937950e0
    v = 0.1802658934377451e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.5984126497885380e0
    b = 0.1816640840360209e0
    v = 0.1849830560443660e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.3791035407695563e0
    b = 0.1720795225656878e0
    v = 0.1713904507106709e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.2778673190586244e0
    b = 0.8213021581932511e-1
    v = 0.1555213603396808e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.5033564271075117e0
    b = 0.8999205842074875e-1
    v = 0.1802239128008525e-2
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count


def ld0770(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 770 point angular Lebedev-Laikov grid."""

    v = 0.2192942088181184e-3
    genOh1(count, coord, weights, v)
    v = 0.1436433617319080e-2
    genOh2(count, coord, weights, v)
    v = 0.1421940344335877e-2
    genOh3(count, coord, weights, v)
    a = 0.5087204410502360e-1
    v = 0.6798123511050502e-3
    genOh4(count, coord, weights, a, v)
    a = 0.1228198790178831e0
    v = 0.9913184235294912e-3
    genOh4(count, coord, weights, a, v)
    a = 0.2026890814408786e0
    v = 0.1180207833238949e-2
    genOh4(count, coord, weights, a, v)
    a = 0.2847745156464294e0
    v = 0.1296599602080921e-2
    genOh4(count, coord, weights, a, v)
    a = 0.3656719078978026e0
    v = 0.1365871427428316e-2
    genOh4(count, coord, weights, a, v)
    a = 0.4428264886713469e0
    v = 0.1402988604775325e-2
    genOh4(count, coord, weights, a, v)
    a = 0.5140619627249735e0
    v = 0.1418645563595609e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6306401219166803e0
    v = 0.1421376741851662e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6716883332022612e0
    v = 0.1423996475490962e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6979792685336881e0
    v = 0.1431554042178567e-2
    genOh4(count, coord, weights, a, v)
    a = 0.1446865674195309e0
    v = 0.9254401499865368e-3
    genOh5(count, coord, weights, a, v)
    a = 0.3390263475411216e0
    v = 0.1250239995053509e-2
    genOh5(count, coord, weights, a, v)
    a = 0.5335804651263506e0
    v = 0.1394365843329230e-2
    genOh5(count, coord, weights, a, v)
    a = 0.6944024393349413e-1
    b = 0.2355187894242326e0
    v = 0.1127089094671749e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.2269004109529460e0
    b = 0.4102182474045730e0
    v = 0.1345753760910670e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.8025574607775339e-1
    b = 0.6214302417481605e0
    v = 0.1424957283316783e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.1467999527896572e0
    b = 0.3245284345717394e0
    v = 0.1261523341237750e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.1571507769824727e0
    b = 0.5224482189696630e0
    v = 0.1392547106052696e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.2365702993157246e0
    b = 0.6017546634089558e0
    v = 0.1418761677877656e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.7714815866765732e-1
    b = 0.4346575516141163e0
    v = 0.1338366684479554e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.3062936666210730e0
    b = 0.4908826589037616e0
    v = 0.1393700862676131e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.3822477379524787e0
    b = 0.5648768149099500e0
    v = 0.1415914757466932e-2
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count


def ld0974(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 974 point angular Lebedev-Laikov grid."""

    v = 0.1438294190527431e-3
    genOh1(count, coord, weights, v)
    v = 0.1125772288287004e-2
    genOh3(count, coord, weights, v)
    a = 0.4292963545341347e-1
    v = 0.4948029341949241e-3
    genOh4(count, coord, weights, a, v)
    a = 0.1051426854086404e0
    v = 0.7357990109125470e-3
    genOh4(count, coord, weights, a, v)
    a = 0.1750024867623087e0
    v = 0.8889132771304384e-3
    genOh4(count, coord, weights, a, v)
    a = 0.2477653379650257e0
    v = 0.9888347838921435e-3
    genOh4(count, coord, weights, a, v)
    a = 0.3206567123955957e0
    v = 0.1053299681709471e-2
    genOh4(count, coord, weights, a, v)
    a = 0.3916520749849983e0
    v = 0.1092778807014578e-2
    genOh4(count, coord, weights, a, v)
    a = 0.4590825874187624e0
    v = 0.1114389394063227e-2
    genOh4(count, coord, weights, a, v)
    a = 0.5214563888415861e0
    v = 0.1123724788051555e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6253170244654199e0
    v = 0.1125239325243814e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6637926744523170e0
    v = 0.1126153271815905e-2
    genOh4(count, coord, weights, a, v)
    a = 0.6910410398498301e0
    v = 0.1130286931123841e-2
    genOh4(count, coord, weights, a, v)
    a = 0.7052907007457760e0
    v = 0.1134986534363955e-2
    genOh4(count, coord, weights, a, v)
    a = 0.1236686762657990e0
    v = 0.6823367927109931e-3
    genOh5(count, coord, weights, a, v)
    a = 0.2940777114468387e0
    v = 0.9454158160447096e-3
    genOh5(count, coord, weights, a, v)
    a = 0.4697753849207649e0
    v = 0.1074429975385679e-2
    genOh5(count, coord, weights, a, v)
    a = 0.6334563241139567e0
    v = 0.1129300086569132e-2
    genOh5(count, coord, weights, a, v)
    a = 0.5974048614181342e-1
    b = 0.2029128752777523e0
    v = 0.8436884500901954e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.1375760408473636e0
    b = 0.4602621942484054e0
    v = 0.1075255720448885e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.3391016526336286e0
    b = 0.5030673999662036e0
    v = 0.1108577236864462e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.1271675191439820e0
    b = 0.2817606422442134e0
    v = 0.9566475323783357e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.2693120740413512e0
    b = 0.4331561291720157e0
    v = 0.1080663250717391e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.1419786452601918e0
    b = 0.6256167358580814e0
    v = 0.1126797131196295e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.6709284600738255e-1
    b = 0.3798395216859157e0
    v = 0.1022568715358061e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.7057738183256172e-1
    b = 0.5517505421423520e0
    v = 0.1108960267713108e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.2783888477882155e0
    b = 0.6029619156159187e0
    v = 0.1122790653435766e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.1979578938917407e0
    b = 0.3589606329589096e0
    v = 0.1032401847117460e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.2087307061103274e0
    b = 0.5348666438135476e0
    v = 0.1107249382283854e-2
    genOh6(count, coord, weights, a, b, v)
    a = 0.4055122137872836e0
    b = 0.5674997546074373e0
    v = 0.1121780048519972e-2
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count


def ld1202(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 1202 point angular Lebedev-Laikov grid."""

    v = 0.1105189233267572e-3
    genOh1(count, coord, weights, v)
    v = 0.9205232738090741e-3
    genOh2(count, coord, weights, v)
    v = 0.9133159786443561e-3
    genOh3(count, coord, weights, v)
    a = 0.3712636449657089e-1
    v = 0.3690421898017899e-3
    genOh4(count, coord, weights, a, v)
    a = 0.9140060412262223e-1
    v = 0.5603990928680660e-3
    genOh4(count, coord, weights, a, v)
    a = 0.1531077852469906e0
    v = 0.6865297629282609e-3
    genOh4(count, coord, weights, a, v)
    a = 0.2180928891660612e0
    v = 0.7720338551145630e-3
    genOh4(count, coord, weights, a, v)
    a = 0.2839874532200175e0
    v = 0.8301545958894795e-3
    genOh4(count, coord, weights, a, v)
    a = 0.3491177600963764e0
    v = 0.8686692550179628e-3
    genOh4(count, coord, weights, a, v)
    a = 0.4121431461444309e0
    v = 0.8927076285846890e-3
    genOh4(count, coord, weights, a, v)
    a = 0.4718993627149127e0
    v = 0.9060820238568219e-3
    genOh4(count, coord, weights, a, v)
    a = 0.5273145452842337e0
    v = 0.9119777254940867e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6209475332444019e0
    v = 0.9128720138604181e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6569722711857291e0
    v = 0.9130714935691735e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6841788309070143e0
    v = 0.9152873784554116e-3
    genOh4(count, coord, weights, a, v)
    a = 0.7012604330123631e0
    v = 0.9187436274321654e-3
    genOh4(count, coord, weights, a, v)
    a = 0.1072382215478166e0
    v = 0.5176977312965694e-3
    genOh5(count, coord, weights, a, v)
    a = 0.2582068959496968e0
    v = 0.7331143682101417e-3
    genOh5(count, coord, weights, a, v)
    a = 0.4172752955306717e0
    v = 0.8463232836379928e-3
    genOh5(count, coord, weights, a, v)
    a = 0.5700366911792503e0
    v = 0.9031122694253992e-3
    genOh5(count, coord, weights, a, v)
    a = 0.9827986018263947e0
    b = 0.1771774022615325e0
    v = 0.6485778453163257e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.9624249230326228e0
    b = 0.2475716463426288e0
    v = 0.7435030910982369e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.9402007994128811e0
    b = 0.3354616289066489e0
    v = 0.7998527891839054e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.9320822040143202e0
    b = 0.3173615246611977e0
    v = 0.8101731497468018e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.9043674199393299e0
    b = 0.4090268427085357e0
    v = 0.8483389574594331e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.8912407560074747e0
    b = 0.3854291150669224e0
    v = 0.8556299257311812e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.8676435628462708e0
    b = 0.4932221184851285e0
    v = 0.8803208679738260e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.8581979986041619e0
    b = 0.4785320675922435e0
    v = 0.8811048182425720e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.8396753624049856e0
    b = 0.4507422593157064e0
    v = 0.8850282341265444e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.8165288564022188e0
    b = 0.5632123020762100e0
    v = 0.9021342299040653e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.8015469370783529e0
    b = 0.5434303569693900e0
    v = 0.9010091677105086e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.7773563069070351e0
    b = 0.5123518486419871e0
    v = 0.9022692938426915e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.7661621213900394e0
    b = 0.6394279634749102e0
    v = 0.9158016174693465e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.7553584143533510e0
    b = 0.6269805509024392e0
    v = 0.9131578003189435e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.7344305757559503e0
    b = 0.6031161693096310e0
    v = 0.9107813579482705e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.7043837184021765e0
    b = 0.5693702498468441e0
    v = 0.9105760258970126e-3
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count


def ld1454(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 1454 point angular Lebedev-Laikov grid."""

    v = 0.7777160743261247e-4
    genOh1(count, coord, weights, v)
    v = 0.7557646413004701e-3
    genOh3(count, coord, weights, v)
    a = 0.3229290663413854e-1
    v = 0.2841633806090617e-3
    genOh4(count, coord, weights, a, v)
    a = 0.8036733271462222e-1
    v = 0.4374419127053555e-3
    genOh4(count, coord, weights, a, v)
    a = 0.1354289960531653e0
    v = 0.5417174740872172e-3
    genOh4(count, coord, weights, a, v)
    a = 0.1938963861114426e0
    v = 0.6148000891358593e-3
    genOh4(count, coord, weights, a, v)
    a = 0.2537343715011275e0
    v = 0.6664394485800705e-3
    genOh4(count, coord, weights, a, v)
    a = 0.3135251434752570e0
    v = 0.7025039356923220e-3
    genOh4(count, coord, weights, a, v)
    a = 0.3721558339375338e0
    v = 0.7268511789249627e-3
    genOh4(count, coord, weights, a, v)
    a = 0.4286809575195696e0
    v = 0.7422637534208629e-3
    genOh4(count, coord, weights, a, v)
    a = 0.4822510128282994e0
    v = 0.7509545035841214e-3
    genOh4(count, coord, weights, a, v)
    a = 0.5320679333566263e0
    v = 0.7548535057718401e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6172998195394274e0
    v = 0.7554088969774001e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6510679849127481e0
    v = 0.7553147174442808e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6777315251687360e0
    v = 0.7564767653292297e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6963109410648741e0
    v = 0.7587991808518730e-3
    genOh4(count, coord, weights, a, v)
    a = 0.7058935009831749e0
    v = 0.7608261832033027e-3
    genOh4(count, coord, weights, a, v)
    a = 0.9955546194091857e0
    v = 0.4021680447874916e-3
    genOh5(count, coord, weights, a, v)
    a = 0.9734115901794209e0
    v = 0.5804871793945964e-3
    genOh5(count, coord, weights, a, v)
    a = 0.9275693732388626e0
    v = 0.6792151955945159e-3
    genOh5(count, coord, weights, a, v)
    a = 0.8568022422795103e0
    v = 0.7336741211286294e-3
    genOh5(count, coord, weights, a, v)
    a = 0.7623495553719372e0
    v = 0.7581866300989608e-3
    genOh5(count, coord, weights, a, v)
    a = 0.5707522908892223e0
    b = 0.4387028039889501e0
    v = 0.7538257859800743e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5196463388403083e0
    b = 0.3858908414762617e0
    v = 0.7483517247053123e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4646337531215351e0
    b = 0.3301937372343854e0
    v = 0.7371763661112059e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4063901697557691e0
    b = 0.2725423573563777e0
    v = 0.7183448895756934e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.3456329466643087e0
    b = 0.2139510237495250e0
    v = 0.6895815529822191e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.2831395121050332e0
    b = 0.1555922309786647e0
    v = 0.6480105801792886e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.2197682022925330e0
    b = 0.9892878979686097e-1
    v = 0.5897558896594636e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.1564696098650355e0
    b = 0.4598642910675510e-1
    v = 0.5095708849247346e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.6027356673721295e0
    b = 0.3376625140173426e0
    v = 0.7536906428909755e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5496032320255096e0
    b = 0.2822301309727988e0
    v = 0.7472505965575118e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4921707755234567e0
    b = 0.2248632342592540e0
    v = 0.7343017132279698e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4309422998598483e0
    b = 0.1666224723456479e0
    v = 0.7130871582177445e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.3664108182313672e0
    b = 0.1086964901822169e0
    v = 0.6817022032112776e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.2990189057758436e0
    b = 0.5251989784120085e-1
    v = 0.6380941145604121e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.6268724013144998e0
    b = 0.2297523657550023e0
    v = 0.7550381377920310e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5707324144834607e0
    b = 0.1723080607093800e0
    v = 0.7478646640144802e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5096360901960365e0
    b = 0.1140238465390513e0
    v = 0.7335918720601220e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4438729938312456e0
    b = 0.5611522095882537e-1
    v = 0.7110120527658118e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.6419978471082389e0
    b = 0.1164174423140873e0
    v = 0.7571363978689501e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5817218061802611e0
    b = 0.5797589531445219e-1
    v = 0.7489908329079234e-3
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count


def ld1730(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 1730 point angular Lebedev-Laikov grid."""

    v = 0.6309049437420976e-4
    genOh1(count, coord, weights, v)
    v = 0.6398287705571748e-3
    genOh2(count, coord, weights, v)
    v = 0.6357185073530720e-3
    genOh3(count, coord, weights, v)
    a = 0.2860923126194662e-1
    v = 0.2221207162188168e-3
    genOh4(count, coord, weights, a, v)
    a = 0.7142556767711522e-1
    v = 0.3475784022286848e-3
    genOh4(count, coord, weights, a, v)
    a = 0.1209199540995559e0
    v = 0.4350742443589804e-3
    genOh4(count, coord, weights, a, v)
    a = 0.1738673106594379e0
    v = 0.4978569136522127e-3
    genOh4(count, coord, weights, a, v)
    a = 0.2284645438467734e0
    v = 0.5435036221998053e-3
    genOh4(count, coord, weights, a, v)
    a = 0.2834807671701512e0
    v = 0.5765913388219542e-3
    genOh4(count, coord, weights, a, v)
    a = 0.3379680145467339e0
    v = 0.6001200359226003e-3
    genOh4(count, coord, weights, a, v)
    a = 0.3911355454819537e0
    v = 0.6162178172717512e-3
    genOh4(count, coord, weights, a, v)
    a = 0.4422860353001403e0
    v = 0.6265218152438485e-3
    genOh4(count, coord, weights, a, v)
    a = 0.4907781568726057e0
    v = 0.6323987160974212e-3
    genOh4(count, coord, weights, a, v)
    a = 0.5360006153211468e0
    v = 0.6350767851540569e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6142105973596603e0
    v = 0.6354362775297107e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6459300387977504e0
    v = 0.6352302462706235e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6718056125089225e0
    v = 0.6358117881417972e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6910888533186254e0
    v = 0.6373101590310117e-3
    genOh4(count, coord, weights, a, v)
    a = 0.7030467416823252e0
    v = 0.6390428961368665e-3
    genOh4(count, coord, weights, a, v)
    a = 0.8354951166354646e-1
    v = 0.3186913449946576e-3
    genOh5(count, coord, weights, a, v)
    a = 0.2050143009099486e0
    v = 0.4678028558591711e-3
    genOh5(count, coord, weights, a, v)
    a = 0.3370208290706637e0
    v = 0.5538829697598626e-3
    genOh5(count, coord, weights, a, v)
    a = 0.4689051484233963e0
    v = 0.6044475907190476e-3
    genOh5(count, coord, weights, a, v)
    a = 0.5939400424557334e0
    v = 0.6313575103509012e-3
    genOh5(count, coord, weights, a, v)
    a = 0.1394983311832261e0
    b = 0.4097581162050343e-1
    v = 0.4078626431855630e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.1967999180485014e0
    b = 0.8851987391293348e-1
    v = 0.4759933057812725e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.2546183732548967e0
    b = 0.1397680182969819e0
    v = 0.5268151186413440e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.3121281074713875e0
    b = 0.1929452542226526e0
    v = 0.5643048560507316e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.3685981078502492e0
    b = 0.2467898337061562e0
    v = 0.5914501076613073e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4233760321547856e0
    b = 0.3003104124785409e0
    v = 0.6104561257874195e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4758671236059246e0
    b = 0.3526684328175033e0
    v = 0.6230252860707806e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5255178579796463e0
    b = 0.4031134861145713e0
    v = 0.6305618761760796e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5718025633734589e0
    b = 0.4509426448342351e0
    v = 0.6343092767597889e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.2686927772723415e0
    b = 0.4711322502423248e-1
    v = 0.5176268945737826e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.3306006819904809e0
    b = 0.9784487303942695e-1
    v = 0.5564840313313692e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.3904906850594983e0
    b = 0.1505395810025273e0
    v = 0.5856426671038980e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4479957951904390e0
    b = 0.2039728156296050e0
    v = 0.6066386925777091e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5027076848919780e0
    b = 0.2571529941121107e0
    v = 0.6208824962234458e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5542087392260217e0
    b = 0.3092191375815670e0
    v = 0.6296314297822907e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.6020850887375187e0
    b = 0.3593807506130276e0
    v = 0.6340423756791859e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4019851409179594e0
    b = 0.5063389934378671e-1
    v = 0.5829627677107342e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4635614567449800e0
    b = 0.1032422269160612e0
    v = 0.6048693376081110e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5215860931591575e0
    b = 0.1566322094006254e0
    v = 0.6202362317732461e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5758202499099271e0
    b = 0.2098082827491099e0
    v = 0.6299005328403779e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.6259893683876795e0
    b = 0.2618824114553391e0
    v = 0.6347722390609353e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5313795124811891e0
    b = 0.5263245019338556e-1
    v = 0.6203778981238834e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5893317955931995e0
    b = 0.1061059730982005e0
    v = 0.6308414671239979e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.6426246321215801e0
    b = 0.1594171564034221e0
    v = 0.6362706466959498e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.6511904367376113e0
    b = 0.5354789536565540e-1
    v = 0.6375414170333233e-3
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count


def ld2030(coord: np.ndarray, weights: np.ndarray, count: Counter):
    """Get 2030 point angular Lebedev-Laikov grid."""

    v = 0.4656031899197431e-4
    genOh1(count, coord, weights, v)
    v = 0.5421549195295507e-3
    genOh3(count, coord, weights, v)
    a = 0.2540835336814348e-1
    v = 0.1778522133346553e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6399322800504915e-1
    v = 0.2811325405682796e-3
    genOh4(count, coord, weights, a, v)
    a = 0.1088269469804125e0
    v = 0.3548896312631459e-3
    genOh4(count, coord, weights, a, v)
    a = 0.1570670798818287e0
    v = 0.4090310897173364e-3
    genOh4(count, coord, weights, a, v)
    a = 0.2071163932282514e0
    v = 0.4493286134169965e-3
    genOh4(count, coord, weights, a, v)
    a = 0.2578914044450844e0
    v = 0.4793728447962723e-3
    genOh4(count, coord, weights, a, v)
    a = 0.3085687558169623e0
    v = 0.5015415319164265e-3
    genOh4(count, coord, weights, a, v)
    a = 0.3584719706267024e0
    v = 0.5175127372677937e-3
    genOh4(count, coord, weights, a, v)
    a = 0.4070135594428709e0
    v = 0.5285522262081019e-3
    genOh4(count, coord, weights, a, v)
    a = 0.4536618626222638e0
    v = 0.5356832703713962e-3
    genOh4(count, coord, weights, a, v)
    a = 0.4979195686463577e0
    v = 0.5397914736175170e-3
    genOh4(count, coord, weights, a, v)
    a = 0.5393075111126999e0
    v = 0.5416899441599930e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6115617676843916e0
    v = 0.5419308476889938e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6414308435160159e0
    v = 0.5416936902030596e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6664099412721607e0
    v = 0.5419544338703164e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6859161771214913e0
    v = 0.5428983656630975e-3
    genOh4(count, coord, weights, a, v)
    a = 0.6993625593503890e0
    v = 0.5442286500098193e-3
    genOh4(count, coord, weights, a, v)
    a = 0.7062393387719380e0
    v = 0.5452250345057301e-3
    genOh4(count, coord, weights, a, v)
    a = 0.7479028168349763e-1
    v = 0.2568002497728530e-3
    genOh5(count, coord, weights, a, v)
    a = 0.1848951153969366e0
    v = 0.3827211700292145e-3
    genOh5(count, coord, weights, a, v)
    a = 0.3059529066581305e0
    v = 0.4579491561917824e-3
    genOh5(count, coord, weights, a, v)
    a = 0.4285556101021362e0
    v = 0.5042003969083574e-3
    genOh5(count, coord, weights, a, v)
    a = 0.5468758653496526e0
    v = 0.5312708889976025e-3
    genOh5(count, coord, weights, a, v)
    a = 0.6565821978343439e0
    v = 0.5438401790747117e-3
    genOh5(count, coord, weights, a, v)
    a = 0.1253901572367117e0
    b = 0.3681917226439641e-1
    v = 0.3316041873197344e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.1775721510383941e0
    b = 0.7982487607213301e-1
    v = 0.3899113567153771e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.2305693358216114e0
    b = 0.1264640966592335e0
    v = 0.4343343327201309e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.2836502845992063e0
    b = 0.1751585683418957e0
    v = 0.4679415262318919e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.3361794746232590e0
    b = 0.2247995907632670e0
    v = 0.4930847981631031e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.3875979172264824e0
    b = 0.2745299257422246e0
    v = 0.5115031867540091e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4374019316999074e0
    b = 0.3236373482441118e0
    v = 0.5245217148457367e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4851275843340022e0
    b = 0.3714967859436741e0
    v = 0.5332041499895321e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5303391803806868e0
    b = 0.4175353646321745e0
    v = 0.5384583126021542e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5726197380596287e0
    b = 0.4612084406355461e0
    v = 0.5411067210798852e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.2431520732564863e0
    b = 0.4258040133043952e-1
    v = 0.4259797391468714e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.3002096800895869e0
    b = 0.8869424306722721e-1
    v = 0.4604931368460021e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.3558554457457432e0
    b = 0.1368811706510655e0
    v = 0.4871814878255202e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4097782537048887e0
    b = 0.1860739985015033e0
    v = 0.5072242910074885e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4616337666067458e0
    b = 0.2354235077395853e0
    v = 0.5217069845235350e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5110707008417874e0
    b = 0.2842074921347011e0
    v = 0.5315785966280310e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5577415286163795e0
    b = 0.3317784414984102e0
    v = 0.5376833708758905e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.6013060431366950e0
    b = 0.3775299002040700e0
    v = 0.5408032092069521e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.3661596767261781e0
    b = 0.4599367887164592e-1
    v = 0.4842744917904866e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4237633153506581e0
    b = 0.9404893773654421e-1
    v = 0.5048926076188130e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4786328454658452e0
    b = 0.1431377109091971e0
    v = 0.5202607980478373e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5305702076789774e0
    b = 0.1924186388843570e0
    v = 0.5309932388325743e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5793436224231788e0
    b = 0.2411590944775190e0
    v = 0.5377419770895208e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.6247069017094747e0
    b = 0.2886871491583605e0
    v = 0.5411696331677717e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.4874315552535204e0
    b = 0.4804978774953206e-1
    v = 0.5197996293282420e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5427337322059053e0
    b = 0.9716857199366665e-1
    v = 0.5311120836622945e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.5943493747246700e0
    b = 0.1465205839795055e0
    v = 0.5384309319956951e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.6421314033564943e0
    b = 0.1953579449803574e0
    v = 0.5421859504051886e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.6020628374713980e0
    b = 0.4916375015738108e-1
    v = 0.5390948355046314e-3
    genOh6(count, coord, weights, a, b, v)
    a = 0.6529222529856881e0
    b = 0.9861621540127005e-1
    v = 0.5433312705027845e-3
    genOh6(count, coord, weights, a, b, v)
    return coord, weights, count
