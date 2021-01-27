# src/units.py

from math import pi, sqrt

__version__ = "2014"


class Units(dict):
    """Dictionary for units that support .attribute access."""

    def __init__(self, *args, **kwargs):
        super(Units, self).__init__(*args, **kwargs)
        self.__dict__ = self


codata = {
    # Source: http://arxiv.org/pdf/1507.07956.pdf
    "2014": {
        "_c": 299792458.0,
        "_mu0": 4.0e-7 * pi,
        "_Grav": 6.67408e-11,
        "_hplanck": 6.626070040e-34,
        "_e": 1.6021766208e-19,
        "_me": 9.10938356e-31,
        "_mp": 1.672621898e-27,
        "_Nav": 6.022140857e23,
        "_k": 1.38064852e-23,
        "_amu": 1.660539040e-27,
    },
}


def createUnits(version):
    """Function that creates a dictionary containing units."""

    # try to obtain codata version
    try:
        u = Units(codata[version])
    except KeyError:
        raise NotImplementedError(
            "Codata version {0} not implemented".format(__version__)
        )

    # permittivity of vacuum
    u["_eps0"] = 1 / u["_mu0"] / u["_c"] ** 2
    # Planck constant / 2pi, J s
    u["_hbar"] = u["_hplanck"] / (2 * pi)
    # Angstrom
    u["Ang"] = u["Angstrom"] = 1.0
    # Nanometer
    u["nm"] = 10.0
    # Bohr
    u["Bohr"] = 4e10 * pi * u["_eps0"] * u["_hbar"] ** 2 / u["_me"] / u["_e"] ** 2
    # Electronvolt
    u["eV"] = 1.0
    # Hartree
    u["Hartree"] = (
        u["_me"] * u["_e"] ** 3 / 16 / pi ** 2 / u["_eps0"] ** 2 / u["_hbar"] ** 2
    )
    # Define kJ/mol
    u["kJ"] = 1000.0 / u["_e"]
    # Define kcal/mol
    u["kcal"] = 4.184 * u["kJ"]
    # Define mol
    u["mol"] = u["_Nav"]
    # Define Rydberg
    u["Rydberg"] = 0.5 * u["Hartree"]
    u["Ry"] = u["Rydberg"]
    # Define Hartree
    u["Ha"] = u["Hartree"]
    # Define second
    u["second"] = 1e10 * sqrt(u["_e"] / u["_amu"])
    # Define femtosecond
    u["fs"] = 1e-15 * u["second"]
    # Define Boltzmann constant, eV/K
    u["kB"] = u["_k"] / u["_e"]
    # Define Pascal J/m^3
    u["Pascal"] = (1 / u["_e"]) / 1e30
    # Define Gigapascal
    u["GPa"] = 1e9 * u["Pascal"]
    # Define Debye
    u["Debye"] = 1.0 / 1e11 / u["_e"] / u["_c"]
    # fine structure constant
    u["alpha"] = u["_e"] ** 2 / (4 * pi * u["_eps0"]) / u["_hbar"] / u["_c"]
    # Inverse centimetre: cm^-1 energy unit
    u["invcm"] = 100 * u["_c"] * u["_hplanck"] / u["_e"]
    # atomic unit of time, s:
    u["_aut"] = u["_hbar"] / (u["alpha"] ** 2 * u["_me"] * u["_c"] ** 2)
    # atomic unit of velocity, m/s:
    u["_auv"] = u["_e"] ** 2 / u["_hbar"] / (4 * pi * u["_eps0"])
    # atomic unit of force, N:
    u["_auf"] = u["alpha"] ** 3 * u["_me"] ** 2 * u["_c"] ** 3 / u["_hbar"]
    # atomic unit of pressure, Pa:
    u["_aup"] = u["alpha"] ** 5 * u["_me"] ** 4 * u["_c"] ** 5 / u["_hbar"] ** 3
    u["AUT"] = u["second"] * u["_aut"]

    # Get SI units
    # metre
    u["m"] = 1e10 * u["Ang"]
    # kilogram
    u["kg"] = 1.0 / u["_amu"]
    # second
    u["s"] = u["second"]
    # ampere
    u["A"] = 1.0 / u["_e"] / u["s"]
    # derived
    u["J"] = u["kJ"] / 1000  # Joule = kg * m**2 / s**2
    u["C"] = 1.0 / u["_e"]  # Coulomb = A * s

    return u


# Initialize expected symbols
# pylint: disable=invalid-name
(
    _Grav,
    _Nav,
    _amu,
    _auf,
    _aup,
    _aut,
    _auv,
    _c,
    _e,
    _eps0,
    _hbar,
    _hplanck,
    _k,
    _me,
    _mp,
    _mu0,
    alpha,
    eV,
    fs,
    invcm,
    kB,
    kJ,
    kcal,
    kg,
    m,
    mol,
    nm,
    s,
    second,
    A,
    AUT,
    Ang,
    Angstrom,
    Bohr,
    C,
    Debye,
    GPa,
    Ha,
    Hartree,
    J,
    Pascal,
    Ry,
    Rydberg,
) = [0.0] * 43


# Update the module scope:
globals().update(createUnits(__version__))
