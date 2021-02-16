# src/kallisto/molecule.py

import numpy as np

from kallisto.atom import Atom


class Molecule(object):
    """The Molecule object.

    Class for representing an isolated molecule.

    Parameters:


    symbols: list of str
        List of symbols of Atoms objects.
    positions: list of xyz-positions
        Atomic positions.
    numbers: list of int
        Atomic nuclear charges.
    charges: list of float
        Atomic charges."""

    def __init__(
        self, symbols=None, positions=None, numbers=None, charges=None,
    ):

        molecule = None

        # Get data from another molecule
        if hasattr(symbols, "get_positions"):
            molecule = symbols
            symbols = None
        # Get data from a list or tuple of Atom objects
        elif (
            isinstance(symbols, (list, tuple))
            and len(symbols) > 0
            and isinstance(symbols[0], Atom)
        ):
            data = [
                [symbol.get_raw(name) for symbol in symbols]
                for name in ["position", "number", "charge"]
            ]
            molecule = self.__class__(None, *data)
            symbols = None

        if molecule is not None:
            # Get data from another molecule
            if symbols is None and numbers is None:
                numbers = molecule.get_atomic_numbers()
            if positions is None:
                positions = molecule.get_positions()

        self.arrays = {}

        if symbols is None:
            if numbers is None:
                if positions is not None:
                    natoms = len(positions)
                else:
                    natoms = 0
                numbers = np.zeros(natoms, int)
            self.new_array(name="numbers", a=numbers, dtype=int)

        if positions is None:
            if "numbers" in self.arrays:
                length = len(self.arrays["numbers"])
            else:
                length = 0
            positions = np.zeros(shape=(length, 3))
        self.new_array(name="positions", a=positions, dtype=float, shape=(3,))

    def new_array(self, name, a, dtype=None, shape=None):
        """Add a new array."""

        if dtype is not None:
            a = np.array(a, dtype, order="C")
            if len(a) == 0 and shape is not None:
                a.shape = (-1,) + shape
        else:
            if isinstance(a, np.ndarray):
                if not a.flags["C_CONTIGUOUS"]:
                    a = np.ascontiguousarray(a)
                else:
                    a = a.copy()
            else:
                a = np.ascontiguousarray(a)

        if name in self.arrays:
            raise RuntimeError('Array "{name}" already present'.format(name=name))

        for b in self.arrays.values():
            if len(a) != len(b):
                raise ValueError(
                    'Array "{name}" has wrong length: {a} != {b}.'.format(
                        name=name, a=len(a), b=len(b)
                    )
                )
            break

        self.arrays[name] = a

    # Setter methods
    def set_array(self, name, a, dtype=None, shape=None):
        """Set array."""

        b = self.arrays.get(name)
        if b is None:
            if a is not None:
                self.new_array(name=name, a=a, dtype=dtype, shape=shape)
        else:
            if a is None:
                del self.arrays[name]
            else:
                a = np.asarray(a)
                if a.shape != b.shape:
                    raise ValueError(
                        'Array "{name}" has wrong shape {a} != {b}'.format(
                            name=name, a=a.shape, b=b.shape
                        )
                    )
                b[:] = a

    # Getter methods
    def get_array(self, name, copy=True):
        """Get an array."""

        if copy:
            return self.arrays[name].copy()
        else:
            return self.arrays[name]

    def get_atomic_numbers(self):
        """Get integer array of atomic numbers."""
        return self.arrays["numbers"].copy()

    def get_positions(self):
        """Get positions-array."""
        return self.arrays["positions"]

    def copy(self):
        """Return a copy."""
        molecule = self.__class__()

        molecule.arrays = {}
        for name, a in self.arrays.items():
            molecule.arrays[name] = a.copy()
        return molecule

    def get_number_of_atoms(self):
        """Get integer number of atoms."""
        return len(self.get_positions())

    def get_bonds(self, partner="X", thresholdBond=0.6, thresholdCN=800.0):
        """Get an index table for covalent bonding partner.

        thresholdBond defines the treshold for a covalent bond."""

        from kallisto.methods import getCovalentBondingPartner

        at = self.get_atomic_numbers()
        coords = self.get_positions()
        return getCovalentBondingPartner(
            at, coords, partner, thresholdBond, thresholdCN
        )

    def get_cns(self, cntype: str, threshold=800.0):
        """Get coordination numbers (cns).

        CN values are calculated for a given structure and are returned as an
        array. Choose functional type by "cn" defining standard (exp), covalent (cov),
        or error (err)."""

        from kallisto.methods import getCoordinationNumbers

        at = self.get_atomic_numbers()
        coords = self.get_positions()
        return getCoordinationNumbers(at, coords, cntype, threshold)

    def get_cnspheres(self, cntype: str, threshold=800.0):
        """Get coordination number spheres (cnsp)."""

        from kallisto.methods import getCoordinationNumberSpheres

        at = self.get_atomic_numbers()
        coords = self.get_positions()
        return getCoordinationNumberSpheres(at, coords, cntype, threshold)

    def get_vdw(self, charge: int, vdwtype: str, scale: float):
        """Get atomic-charge dependent van der Waals radii (vdws).

        VDW values are calculated for a given structure and are returned as an
        array. For the charge dependency EEQ atomic partial charges are used
        in an empirical scaling function as used in the dftd4 program."""

        from kallisto.methods import getVanDerWaalsRadii

        at = self.get_atomic_numbers()
        nat = self.get_number_of_atoms()
        aiw = self.get_alp(charge=charge)
        return getVanDerWaalsRadii(nat, at, aiw, charge, vdwtype, scale)

    def get_alp(self, charge: int):
        """Get atomic-charge dependent dynamic atomic polarizabilities (alps).

        ALP values are calculated for a given structure and are returned as an
        array. For the charge dependency EEQ atomic partial charges are used
        in an empirical scaling function as used in the dftd4 program."""

        from kallisto.methods import getPolarizabilities

        at = self.get_atomic_numbers()
        covcn = self.get_cns(cntype="cov")
        qs = self.get_eeq(charge)
        return getPolarizabilities(at, covcn, qs, charge)

    def get_eeq(self, charge: int):
        """Get atomic electronegativity equilibration partial charges (eeqs).

        EEQ values are calculated for a given structure and are returned as an
        array."""

        from kallisto.methods import getAtomicPartialCharges

        at = self.get_atomic_numbers()
        coords = self.get_positions()
        cns = self.get_cns(cntype="cov")
        return getAtomicPartialCharges(at, coords, cns, charge)

    def writeMolecule(self, name: str):
        """Write molecular structure."""

        import os
        from kallisto.data import chemical_symbols
        from kallisto.units import Bohr

        coord = self.get_positions()
        at = self.get_atomic_numbers()
        nat = self.get_number_of_atoms()

        f = open(name, "w")
        s = os.linesep
        f.write("{:5}".format(nat) + s)
        f.write("Created with kallisto" + s)
        for i in range(nat):
            f.write(
                "{:3} {:9.4f} {:9.4f} {:9.4f}".format(
                    chemical_symbols[at[i]],
                    coord[i][0] * Bohr,
                    coord[i][1] * Bohr,
                    coord[i][2] * Bohr,
                )
                + s
            )
        f.close()
