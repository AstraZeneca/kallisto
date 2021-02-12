# src/kallisto/atom.py

import numpy as np

from kallisto.data import atomic_numbers
from kallisto.data import chemical_symbols
from kallisto.utils import basestring


def atomProperty(name: str, doc: str):
    """Helper function to create Atom property."""

    def getter(self):
        return self.get(name)

    def setter(self, value):
        self.set(name, value)

    def deleter(self):
        pass

    return property(getter, setter, deleter, doc)


def xyzProperty(index: int):
    """Helper function to create Atom XYZ-property."""

    def getter(self):
        return self.position[index]

    def setter(self, value):
        self.position[index] = value

    return property(getter, setter, doc="XYZ"[index] + "-coordinate")


class Atom(object):
    """The Atom object.

    Class for representing a single atom.

    Parameters:

    symbol: str or int
        Can be a chemical symbol (str) or an atomic number (int).
    position: sequence of 3 floats
        Atomic position.
    charge: float
        Atomic charge.
    """

    __slots__ = ["data", "molecule", "index"]

    def __init__(self, symbol="X", position=(0, 0, 0), charge=0, molecule=None):

        self.data = data = {}

        if molecule is None:
            # This atom is not part of a molecule
            if isinstance(symbol, basestring):
                data["number"] = atomic_numbers[symbol]
            else:
                data["number"] = symbol
            data["position"] = np.array(position, float)
            data["charge"] = charge

        self.molecule = molecule

    def get_raw(self, name):
        """Get name attribute, return None if not explicitely set."""
        if name == "symbol":
            return chemical_symbols[self.get_raw("number")]

        if self.molecule is None:
            return self.data[name]

        return None

    def get(self, name):
        """Get name attribute, return default if not explicitely set."""
        value = self.get_raw(name)
        return value

    def set(self, name, value):
        """Set name attribute to value."""
        if name == "symbol":
            name = "number"
            value = atomic_numbers[value]

        if self.molecule is None:
            self.data[name] = value

    symbol = atomProperty("symbol", "Chemical symbol")
    number = atomProperty("number", "Atomic number")
    position = atomProperty("position", "XYZ-coordinates")
    charge = atomProperty("charge", "Initial atomic charge")
    x = xyzProperty(0)
    y = xyzProperty(1)
    z = xyzProperty(2)
