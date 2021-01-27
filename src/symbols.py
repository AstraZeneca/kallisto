# src/symbols.py

from data import atomic_numbers
from formula import Formula


def string2symbols(s):
    """Convert string to list of chemical symbols."""
    print(s, list(Formula(s)))
    return list(Formula(s))


def symbols2numbers(symbols):
    if isinstance(symbols, str):
        symbols = string2symbols(symbols)
    numbers = []
    for s in symbols:
        if isinstance(s, str):
            numbers.append(atomic_numbers[s])
        else:
            numbers.append(int(s))
    return numbers
