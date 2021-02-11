# src/kallisto/units.py

from math import pi, sqrt

version = "arxiv2014"


class Units(dict):
    """Dictionary for units."""

    def __init__(self, *args, **kwargs):
        super(Units, self).__init__(*args, **kwargs)
        self.__dict__ = self


codata = {
    # Reference: http://arxiv.org/pdf/1507.07956.pdf
    "arxiv2014": {
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
    """Create a dictionary with units."""

    try:
        units = Units(codata[version])
    except KeyError:
        raise NotImplementedError("Codata {0} not found".format(version))

    # Define all units
    # permittivity of vacuum
    units["_eps0"] = 1 / units["_mu0"] / units["_c"] ** 2
    # Planck constant / 2pi, J s
    units["_hbar"] = units["_hplanck"] / (2 * pi)
    # Angstrom
    units["Ang"] = units["Angstrom"] = 1.0
    # Nanometer
    units["nm"] = 10.0
    # Bohr
    units["Bohr"] = (
        4e10
        * pi
        * units["_eps0"]
        * units["_hbar"] ** 2
        / units["_me"]
        / units["_e"] ** 2
    )
    # Electronvolt
    units["eV"] = 1.0
    # Hartree
    units["Hartree"] = (
        units["_me"]
        * units["_e"] ** 3
        / 16
        / pi ** 2
        / units["_eps0"] ** 2
        / units["_hbar"] ** 2
    )
    # Define kJ/mol
    units["kJ"] = 1000.0 / units["_e"]
    # Define kcal/mol
    units["kcal"] = 4.184 * units["kJ"]
    # Define mol
    units["mol"] = units["_Nav"]
    # Define Rydberg
    units["Rydberg"] = 0.5 * units["Hartree"]
    units["Ry"] = units["Rydberg"]
    # Define Hartree
    units["Ha"] = units["Hartree"]
    # Define second
    units["second"] = 1e10 * sqrt(units["_e"] / units["_amu"])
    # Define femtosecond
    units["fs"] = 1e-15 * units["second"]
    # Define Boltzmann constant, eV/K
    units["kB"] = units["_k"] / units["_e"]
    # Define Pascal J/m^3
    units["Pascal"] = (1 / units["_e"]) / 1e30
    # Define Gigapascal
    units["GPa"] = 1e9 * units["Pascal"]
    # Define Debye
    units["Debye"] = 1.0 / 1e11 / units["_e"] / units["_c"]
    # fine structure constant
    units["alpha"] = (
        units["_e"] ** 2 / (4 * pi * units["_eps0"]) / units["_hbar"] / units["_c"]
    )
    # Inverse centimetre: cm^-1 energy unit
    units["invcm"] = 100 * units["_c"] * units["_hplanck"] / units["_e"]
    # atomic unit of time, s:
    units["_aut"] = units["_hbar"] / (
        units["alpha"] ** 2 * units["_me"] * units["_c"] ** 2
    )
    # atomic unit of velocity, m/s:
    units["_auv"] = units["_e"] ** 2 / units["_hbar"] / (4 * pi * units["_eps0"])
    # atomic unit of force, N:
    units["_auf"] = (
        units["alpha"] ** 3 * units["_me"] ** 2 * units["_c"] ** 3 / units["_hbar"]
    )
    # atomic unit of pressure, Pa:
    units["_aup"] = (
        units["alpha"] ** 5 * units["_me"] ** 4 * units["_c"] ** 5 / units["_hbar"] ** 3
    )
    units["AUT"] = units["second"] * units["_aut"]
    # Define SI units
    # metre
    units["m"] = 1e10 * units["Ang"]
    # kilogram
    units["kg"] = 1.0 / units["_amu"]
    # second
    units["s"] = units["second"]
    # ampere
    units["A"] = 1.0 / units["_e"] / units["s"]
    # derived
    units["J"] = units["kJ"] / 1000  # Joule = kg * m**2 / s**2
    units["C"] = 1.0 / units["_e"]  # Coulomb = A * s

    return units


# Initialize symbols
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
globals().update(createUnits(version))
