# src/kallisto/console.py

from typing import Tuple

import click
import numpy as np

import kallisto.reader.strucreader as ksr
from kallisto.utils import errorbye, silentPrinter


class Config(object):
    """Define global config file for click."""

    def __init__(self):
        self.silent = False
        self.shift = 0


pass_config = click.make_pass_decorator(Config, ensure=True)


@click.group(chain=True)
@click.option("--silent", is_flag=True)
@click.option("--shift", default=0, type=int, required=False)
@pass_config
def cli(config, silent: bool, shift: int):
    """kallisto calculates quantum mechanically derived atomic features.\n

    Please check out the documentation (https://ehjc.gitbook.io/kallisto/).\n

    Always cite:\n
     E. Caldeweyher J. Open Source Softw., 2021, 6, 3050.
     (https://doi.org/10.21105/joss.03050)
    """

    config.shift = shift
    config.silent = silent


@cli.command("cns")
@pass_config
@click.option(
    "--cntype",
    default="erf",
    type=str,
    show_default=True,
    help="Coordination number type (exp, cov, err).",
)
@click.option(
    "--out",
    default="-",
    type=click.File("w"),
    show_default=True,
    required=False,
    help="Write to output file.",
)
@click.argument("inp", type=str, default="coord", required=True)
def cns(config, inp: str, out: click.File, cntype: str):
    """Atomic coordination numbers."""

    # Available CNs
    availableCN = ("erf", "cov", "exp")
    if cntype not in availableCN:
        errorbye(
            'CN definition "{}" is not implemented. Please use "erf", "cov", or "exp"'.format(
                cntype
            )
        )

    molecule = ksr.constructMolecule(geometry=inp, out=out)
    cns = molecule.get_cns(cntype)
    nat = molecule.get_number_of_atoms()
    for i in range(nat):
        silentPrinter(config.silent, cns[i], out)

    return cns


@cli.command("prox")
@pass_config
@click.option(
    "--size",
    type=(int, int),
    default=(2, 3),
    show_default=True,
    required=False,
    help="Size that defines the proximity shells.",
)
@click.option(
    "--out",
    default="-",
    type=click.File("w"),
    show_default=True,
    required=False,
    help="Write to output file.",
)
@click.argument("inp", type=str, default="coord", required=True)
def prox(config, inp: str, size: Tuple[int, int], out: click.File):
    """Atomic proximity shells."""

    # Stop if outer border is smaller than inner one
    if size[0] > size[1]:
        errorbye("Outer border is smaller than inner one. Switch them and try again!")

    molecule = ksr.constructMolecule(geometry=inp, out=out)
    nat = molecule.get_number_of_atoms()
    prox = molecule.get_prox(size)
    for i in range(nat):
        silentPrinter(config.silent, prox[i], out)

    return prox


@cli.command("bonds")
@pass_config
@click.option(
    "--partner",
    type=str,
    default="X",
    show_default=True,
    required=True,
    help="Get all partner of atom X.",
)
@click.option("--constrain", is_flag=True)
@click.option(
    "--out",
    default="-",
    type=click.File("w"),
    show_default=True,
    required=False,
    help="Write to output file.",
)
@click.argument("inp", type=str, default="coord", required=True)
def bonds(config, inp: str, partner: int, constrain: bool, out: click.File):
    """Get information about covalent bonding partner."""

    import os

    molecule = ksr.constructMolecule(geometry=inp, out=out)

    if partner == "X":
        # Get index table of covalent bonding partners
        bonds = molecule.get_bonds(partner=partner)
    else:
        # Get all covalent bonding partners of atom #partner
        bonds = molecule.get_bonds(partner=partner)

    nat = molecule.get_number_of_atoms()
    # write constrain file in xtb format
    if constrain:
        f = open("constrain.inp", "w")
        s = os.linesep
        f.write("$constrain" + s)
        for i in range(nat):
            for partner in bonds[i]:
                f.write(
                    " distance: {}, {}, auto".format(
                        i + 1 + config.shift, partner + 1 + config.shift
                    )
                    + s
                )
        f.write("$end" + s)
        f.close()

    if partner != "X":
        silentPrinter(config.silent, str(bonds), out)
    else:
        for i in range(nat):
            silentPrinter(config.silent, str(bonds[i]), out)

    return bonds


@cli.command("sort")
@pass_config
@click.option(
    "--start",
    type=str,
    default="X",
    show_default=True,
    required=True,
    help="Get all partner of atom X.",
)
@click.option(
    "--out",
    default="-",
    type=click.File("w"),
    show_default=True,
    required=False,
    help="Write to output file.",
)
@click.argument("inp", type=str, default="coord", required=True)
def sort(config, inp: str, start: str, out: click.File):
    """Sort input geoemtry according to connectivity.

    start defines on which atom we start the sorting process.
    """

    molecule = ksr.constructMolecule(geometry=inp, out=out)
    bonds = molecule.get_bonds()
    nat = molecule.get_number_of_atoms()

    # construct graph
    from kallisto.sort import Graph

    g = Graph(inp, out)

    for i in range(nat):
        partners = bonds[i]
        for j in partners:
            # add edges
            g.addEdge(j, i)

    if start == "X":
        # start at #0
        g.BFS(0)
    else:
        # start at atom #start
        g.BFS(int(start))


@cli.command("eeq")
@pass_config
@click.option(
    "--chrg",
    default=0,
    type=int,
    show_default=True,
    help="Absolute charge of the system.",
)
@click.option(
    "--out",
    default="-",
    type=click.File("w"),
    show_default=True,
    required=False,
    help="Write to output file.",
)
@click.argument("inp", type=str, default="coord", required=True)
def eeq(config, inp: str, out: click.File, chrg: int):
    """Electronegativity equilibration atomic partial charges."""

    molecule = ksr.constructMolecule(geometry=inp, out=out)
    nat = molecule.get_number_of_atoms()
    eeq = molecule.get_eeq(chrg)
    for i in range(nat):
        silentPrinter(config.silent, eeq[i], out)

    return eeq


@cli.command("alp")
@pass_config
@click.option(
    "--chrg",
    default=0,
    type=int,
    show_default=True,
    help="Absolute charge of the system.",
)
@click.option("--molecular", is_flag=True)
@click.option(
    "--out",
    default="-",
    type=click.File("w"),
    show_default=True,
    required=False,
    help="Write to output file.",
)
@click.argument("inp", type=str, default="coord", required=True)
def alp(config, inp: str, out: click.File, chrg: int, molecular: bool):
    """Static atomic polarizabilities in Bohr^3."""

    molecule = ksr.constructMolecule(geometry=inp, out=out)
    nat = molecule.get_number_of_atoms()
    alp = molecule.get_alp(charge=chrg)
    if molecular:
        silentPrinter(config.silent, np.sum(alp), out)
    else:
        for i in range(nat):
            silentPrinter(config.silent, alp[i], out)

    return alp


@cli.command("vdw")
@pass_config
@click.option(
    "--chrg",
    default=0,
    type=int,
    show_default=True,
    required=True,
    help="Absolute charge of the system.",
)
@click.option(
    "--vdwtype",
    default="rahm",
    type=str,
    show_default=True,
    required=True,
    help="Scaling of van der Waals radii.",
)
@click.option("--angstrom", is_flag=True)
@click.option(
    "--out",
    default="-",
    type=click.File("w"),
    show_default=True,
    required=False,
    help="Write to output file.",
)
@click.argument("inp", type=str, default="coord", required=True)
def vdw(config, inp: str, out: click.File, chrg: int, vdwtype: str, angstrom: bool):
    """Charge-dependent atomic van der Waals radii in Bohr."""

    # Available VDWs
    availableVDW = ("rahm", "truhlar")
    if vdwtype not in availableVDW:
        errorbye(
            'VDW definition "{}" is not implemented. Please use "rahm" or "truhlar"'.format(
                vdwtype
            )
        )

    molecule = ksr.constructMolecule(geometry=inp, out=out)
    nat = molecule.get_number_of_atoms()

    if angstrom:
        from kallisto.units import Bohr

        scale = Bohr
    else:
        scale = 1.0

    vdw = molecule.get_vdw(chrg, vdwtype, scale)
    for i in range(nat):
        silentPrinter(config.silent, vdw[i], out)

    return vdw


@cli.command("rms")
@pass_config
@click.option(
    "--out",
    default="-",
    type=click.File("w"),
    show_default=True,
    required=False,
    help="Write to output file.",
)
@click.argument("inp", type=(str, str), default=("coord1", "coord2"), required=True)
def rms(config, inp: Tuple[str, str], out: click.File):
    """Calculate the root mean squared deviation between two structures using quaternions.
    Based on a Fortran implementation by Chaok Seok, Evangelos
    Coutsias, and Ken Dill."""

    from kallisto.rmsd import rmsd

    mol1 = ksr.constructMolecule(geometry=inp[0], out=out)
    nat1 = mol1.get_number_of_atoms()
    mol2 = ksr.constructMolecule(geometry=inp[1], out=out)
    nat2 = mol2.get_number_of_atoms()

    # for RMSD comparison both coordinates need the same atom count
    if nat1 != nat2:
        errorbye(
            "Error: number of atoms do not match in {} and in {}".format(
                inp[0], inp[1]
            ),
        )

    coord1 = mol1.get_positions()
    coord2 = mol2.get_positions()

    # get RMSD error and rotation matrix u
    error, u = rmsd(nat1, coord1, coord2)

    silentPrinter(config.silent, "RMSD {} Angstrom".format(error), out)
    silentPrinter(config.silent, "Rotation Matrix", out)
    click.echo(u, file=out)  # type: ignore

    return error, u


@cli.command("lig")
@pass_config
@click.option(
    "--center",
    type=int,
    default=0,
    show_default=True,
    required=True,
    help="Central atom for which all bonding partner (ligands) are defined.",
)
@click.option(
    "--out",
    default="-",
    type=click.File("w"),
    show_default=True,
    required=False,
    help="Write to output file.",
)
@click.argument("inp", type=str, default="coord", required=True)
def lig(config, inp: str, center: int, out: click.File):
    """Get all substructures (or ligands) that are bound to the center atom."""

    # setup reference molecular structure
    ref = ksr.constructMolecule(geometry=inp, out=out)
    nat = ref.get_number_of_atoms()

    # get all covalent bonding partner in reference complex
    covbonds = ref.get_bonds()

    from kallisto.rmsd import recursiveGetSubstructures

    silentPrinter(config.silent, "Write out substructures for {}".format(center), out)

    substructures = recursiveGetSubstructures(nat, covbonds, center)  # type: ignore

    k = 0
    for path in substructures:
        silentPrinter(
            config.silent,
            "Substructure {}: {}".format(k, path),
            out,
        )
        k += 1


@cli.command("exs")
@pass_config
@click.option(
    "--center",
    type=int,
    default=0,
    show_default=True,
    required=True,
    help="Central metal atom in transition-metal complex.",
)
@click.option(
    "--subnr",
    type=int,
    default=0,
    show_default=True,
    required=True,
    help="Number of substructure to be exchanged.",
)
@click.option(
    "--name",
    type=str,
    default="newstructure",
    show_default=True,
    help="Name of the output xyz file.",
)
@click.option(
    "--rotate",
    type=int,
    default=0,
    show_default=True,
    help="Rotate new substrate around covalent bond to center around specified degree .",
)
@click.option("--exclude", is_flag=True)
@click.option(
    "--out",
    default="-",
    type=click.File("w"),
    show_default=True,
    required=False,
    help="Write to output file.",
)
@click.argument("inp", type=(str, str), default=("coord1", "coord2"), required=True)
def exs(
    config,
    inp: str,
    center: int,
    subnr: int,
    name: str,
    exclude: bool,
    rotate: int,
    out: click.File,
):
    """Exchange a substrate within a transition metal complex with another
    substrate. Use an root mean squared deviation (RMSD) measure to rotate
    and shift new substrate to the position of the old substrate.
    Calculate the RMSD between two structures using quaternions."""

    # setup reference molecular structure
    ref = ksr.constructMolecule(geometry=inp[0], out=out)
    substrate = ksr.constructMolecule(geometry=inp[1], out=out)
    nat = ref.get_number_of_atoms()

    # get all covalent bonding partner in reference complex
    covBonds = ref.get_bonds()

    # get covalent bonds in new substrate
    newSubBonds = substrate.get_bonds()
    newSubBonds = np.array(newSubBonds, dtype=object)

    from kallisto.rmsd import exchangeSubstructure

    mol = exchangeSubstructure(
        nat,
        center,
        subnr,
        covBonds,  # type: ignore
        ref,
        substrate,
        newSubBonds,
        name,
        rotate,
        exclude,
    )

    return mol


@cli.command("stm")
@pass_config
@click.option(
    "--origin",
    type=int,
    default=0,
    show_default=True,
    required=True,
    help="Origin atom.",
)
@click.option(
    "--partner",
    type=int,
    default=1,
    show_default=True,
    required=True,
    help="Partner atom.",
)
@click.option(
    "--out",
    default="-",
    type=click.File("w"),
    show_default=True,
    required=False,
    help="Write to output file.",
)
@click.argument("inp", type=str, default="coord", required=True)
def stm(config, inp: str, origin: int, partner: int, out: click.File):
    """Calculate sterimol descriptors using kallisto van der Waals radii."""

    # setup molecular structure
    mol = ksr.constructMolecule(geometry=inp, out=out)

    # calculate Sterimol descriptors: L, bmin, bmax
    from kallisto.sterics import getClassicalSterimol

    L, bmin, bmax = getClassicalSterimol(mol, origin, partner)

    # print origin and partner
    silentPrinter(
        config.silent,
        "Calculated for atom {0} (origin) and atom {1} (partner)".format(
            origin, partner
        ),
        out,
    )

    # descriptors in Bohr
    silentPrinter(
        config.silent,
        "L, Bmin, Bmax / au: {:8.6f} {:8.6f} {:8.6f}".format(L, bmin, bmax),
        out,
    )

    # descriptors in Angstrom
    from kallisto.units import Bohr

    silentPrinter(
        config.silent,
        "L, Bmin, Bmax / A: {:8.6f} {:8.6f} {:8.6f}".format(
            L * Bohr, bmin * Bohr, bmax * Bohr
        ),
        out,
    )

    return L, bmin, bmax
