# src/atom.py

from data import atomic_numbers
from data import chemical_symbols
import numpy as np
from utils import basestring


def atomproperty(name, doc):
    """Helper function to easily create Atom attribute property."""

    def getter(self):
        return self.get(name)

    def setter(self, value):
        self.set(name, value)

    def deleter(self):
        self.delete(name)

    return property(getter, setter, deleter, doc)


def xyzproperty(index):
    """Helper function to easily create Atom XYZ-property."""

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

        self.data = d = {}

        if molecule is None:
            # This atom is not part of a molecule
            if isinstance(symbol, basestring):
                d["number"] = atomic_numbers[symbol]
            else:
                d["number"] = symbol
            d["position"] = np.array(position, float)
            d["charge"] = charge

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

    def delete(self, name):
        """Delete name attribute."""
        if (self.molecule is None) and (name not in ["number", "symbol", "position"]):
            self.data[name] = None

    symbol = atomproperty("symbol", "Chemical symbol")
    position = atomproperty("position", "XYZ-coordinates")
    charge = atomproperty("charge", "Initial atomic charge")
    x = xyzproperty(0)
    y = xyzproperty(1)
    z = xyzproperty(2)
