# src/kallisto/console.py

from typing import Tuple

import click
import numpy as np

import kallisto.reader.strucreader as ksr
from kallisto.utils import errorbye, verbosePrinter


class Config(object):
    """Define global config file for click."""

    def __init__(self):
        self.verbose = False
        self.context = click.get_current_context()
        self.shift = 0


pass_config = click.make_pass_decorator(Config, ensure=True)


@click.group(chain=True)
@click.option("--verbose", is_flag=True)
@click.option("--shift", default=0, type=int, required=False)
@pass_config
def cli(config, verbose: bool, shift: int):
    """kallisto calculates quantum mechanically derived atomic features.\n

    Please check out the documentation.
    """

    config.shift = shift
    config.verbose = verbose
    pass


@cli.command("cns")
@pass_config
@click.option(
    "--inp",
    default="coord",
    type=str,
    show_default=True,
    required=True,
    help="Geometry file in xyz or turbomole format.",
)
@click.option(
    "--cntype",
    default="erf",
    type=str,
    show_default=True,
    help="Coordination number type (exp, cov, err).",
)
@click.argument("out", type=click.File("w"), default="-", required=False)
def cns(config, inp: str, out: click.File, cntype: str):
    """Atomic coordination numbers."""

    molecule = ksr.constructMolecule(geometry=inp, out=out)
    cns = molecule.get_cns(cntype)
    verbosePrinter(config.verbose, cns, out)

    return cns


@cli.command("cnsp")
@pass_config
@click.option(
    "--inp",
    default="coord",
    type=str,
    show_default=True,
    required=True,
    help="Geometry file in xyz or turbomole format.",
)
@click.option(
    "--cntype",
    default="erf",
    type=str,
    show_default=True,
    help="Coordination number type (exp, cov, err).",
)
@click.argument("out", type=click.File("w"), default="-", required=False)
def cnsp(config, inp: str, out: click.File, cntype: str):
    """Atomic coordination number spheres."""

    molecule = ksr.constructMolecule(geometry=inp, out=out)
    nat = molecule.get_number_of_atoms()
    cnsp = molecule.get_cnspheres(cntype)
    for i in range(nat):
        verbosePrinter(config.verbose, cnsp[i], out)

    return cnsp


@cli.command("bonds")
@pass_config
@click.option(
    "--inp",
    default="coord",
    type=str,
    show_default=True,
    required=True,
    help="Geometry file in xyz or turbomole format.",
)
@click.option(
    "--partner",
    type=str,
    default="X",
    show_default=True,
    required=True,
    help="Get all partner of atom X.",
)
@click.option("--constrain", is_flag=True)
@click.argument("out", type=click.File("w"), default="-", required=False)
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

    # write constrain file in xtb format
    if constrain:
        nat = molecule.get_number_of_atoms()
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
        pass

    return bonds


@cli.command("sort")
@pass_config
@click.option(
    "--inp",
    default="coord",
    type=str,
    show_default=True,
    required=True,
    help="Geometry file in xyz or turbomole format.",
)
@click.option(
    "--start",
    type=str,
    default="X",
    show_default=True,
    required=True,
    help="Get all partner of atom X.",
)
@click.argument("out", type=click.File("w"), default="-", required=False)
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
    "--inp",
    default="coord",
    type=str,
    show_default=True,
    required=True,
    help="Geometry file in xyz or turbomole format.",
)
@click.option(
    "--chrg",
    default=0,
    type=int,
    show_default=True,
    help="Absolute charge of the system.",
)
@click.argument("out", type=click.File("w"), default="-", required=False)
def eeq(config, inp: str, out: click.File, chrg: int):
    """Electronegativity equilibration atomic partial charges."""

    molecule = ksr.constructMolecule(geometry=inp, out=out)
    nat = molecule.get_number_of_atoms()
    eeq = molecule.get_eeq(chrg)
    for i in range(nat):
        verbosePrinter(config.verbose, eeq[i], out)

    return eeq


@cli.command("alp")
@pass_config
@click.option(
    "--inp",
    default="coord",
    type=str,
    show_default=True,
    help="Calculate atomic polarizabilities (static and dynamic) in Bohr^3.",
)
@click.option(
    "--chrg",
    default=0,
    type=int,
    show_default=True,
    help="Absolute charge of the system.",
)
@click.option("--molecular", is_flag=True)
@click.argument("out", type=click.File("w"), default="-", required=False)
def alp(config, inp: str, out: click.File, chrg: int, molecular: bool):
    """Static atomic polarizabilities in Bohr^3."""

    molecule = ksr.constructMolecule(geometry=inp, out=out)
    nat = molecule.get_number_of_atoms()
    alp = molecule.get_alp(charge=chrg)
    if molecular:
        verbosePrinter(config.verbose, np.sum(alp), out)
    else:
        for i in range(nat):
            verbosePrinter(config.verbose, alp[i], out)

    return alp


@cli.command("vdw")
@pass_config
@click.option(
    "--inp",
    default="coord",
    type=str,
    show_default=True,
    required=True,
    help="Calculate atomic van der Waals radii in Bohr.",
)
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
@click.argument("out", type=click.File("w"), default="-", required=False)
def vdw(config, inp: str, out: click.File, chrg: int, vdwtype: str, angstrom: bool):
    """Charge-dependent atomic van der Waals radii in Bohr."""

    molecule = ksr.constructMolecule(geometry=inp, out=out)
    nat = molecule.get_number_of_atoms()
    vdw = np.zeros(shape=(nat,), dtype=np.float64)

    if angstrom:
        from kallisto.units import Bohr

        scale = Bohr
    else:
        scale = 1.0

    vdw = np.zeros(shape=(nat,), dtype=np.float64)
    vdw = molecule.get_vdw(chrg, vdwtype, scale)
    for i in range(nat):
        verbosePrinter(config.verbose, vdw[i], out)

    return vdw


@cli.command("rms")
@pass_config
@click.option(
    "--compare",
    type=(str, str),
    show_default=True,
    help="Coordinates to be compared.",
    required=True,
)
@click.argument("out", type=click.File("w"), default="-", required=False)
def rms(config, compare: Tuple[str, str], out: click.File):
    """Calculate the root mean squared deviation between two structures using quaternions.
    Based on a Fortran implementation by Chaok Seok, Evangelos
    Coutsias, and Ken Dill."""

    from kallisto.rmsd import rmsd

    mol1 = ksr.constructMolecule(geometry=compare[0], out=out)
    nat1 = mol1.get_number_of_atoms()
    mol2 = ksr.constructMolecule(geometry=compare[1], out=out)
    nat2 = mol2.get_number_of_atoms()

    # for RMSD comparison both coordinates need the same atom count
    if nat1 != nat2:
        errorbye(
            "Error: number of atoms do not match in {} and in {}".format(
                compare[0], compare[1]
            ),
            out,
        )

    coord1 = mol1.get_positions()
    coord2 = mol2.get_positions()

    # get RMSD error and rotation matrix u
    error, u = rmsd(nat1, coord1, coord2)

    verbosePrinter(config.verbose, "RMSD {} Angstrom".format(error), out)
    verbosePrinter(config.verbose, "Rotation Matrix", out)
    click.echo(u, file=out)  # type: ignore

    return error, u


@cli.command("lig")
@pass_config
@click.option(
    "--inp",
    type=str,
    show_default=True,
    required=True,
    help="Input structure to determine covalent substructures for center atom.",
)
@click.option(
    "--center",
    type=int,
    default=0,
    show_default=True,
    required=True,
    help="Central atom for which all bonding partner (ligands) are defined.",
)
@click.argument("out", type=click.File("w"), default="-", required=False)
def lig(config, inp: str, center: int, out: click.File):
    """Get all substructures (or ligands) that are bound to the center atom."""

    # setup reference molecular structure
    ref = ksr.constructMolecule(geometry=inp, out=out)
    nat = ref.get_number_of_atoms()

    # get all covalent bonding partner in reference complex
    covbonds = config.context.invoke(bonds, inp=inp)

    from kallisto.rmsd import recursiveGetSubstructures

    verbosePrinter(config.verbose, "Write out substructures for {}".format(center), out)

    substructures = recursiveGetSubstructures(nat, covbonds, center)

    if config.verbose:
        k = 0
        for path in substructures:
            verbosePrinter(
                config.verbose, "Substructure {}: {}".format(k, path), out,
            )
            k += 1


@cli.command("exs")
@pass_config
@click.option(
    "--inp",
    type=(str, str),
    show_default=True,
    required=True,
    help="Reference structure and substrate.",
)
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
@click.argument("out", type=click.File("w"), default="-", required=False)
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
    covBonds = config.context.invoke(bonds, inp=inp[0])
    # get covalent bonds in new substrate
    newSubBonds = substrate.get_bonds(partner="X")
    newSubBonds = np.array(newSubBonds, dtype=object)

    from kallisto.rmsd import exchangeSubstructure

    exchangeSubstructure(
        nat,
        center,
        subnr,
        covBonds,
        ref,
        substrate,
        newSubBonds,
        name,
        rotate,
        exclude,
    )


@cli.command("stm")
@pass_config
@click.option(
    "--inp",
    type=str,
    show_default=True,
    required=True,
    help="Reference structure and substrate.",
)
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
    default=0,
    show_default=True,
    required=True,
    help="Partner atom.",
)
@click.argument("out", type=click.File("w"), default="-", required=False)
def stm(config, inp: str, origin: int, partner: int, out: click.File):
    """Calculate sterimol descriptors using kallisto van der Waals radii."""

    # setup molecular structure
    mol = ksr.constructMolecule(geometry=inp, out=out)

    # calculate Sterimol descriptors: L, bmin, bmax
    from kallisto.sterics import getClassicalSterimol

    L, bmin, bmax = getClassicalSterimol(mol, origin, partner)

    if config.verbose:
        # print values in Bohr
        verbosePrinter(
            config.verbose,
            "L, Bmin, Bmax / au: {:5.2f} {:5.2f} {:5.2f}".format(L, bmin, bmax),
            out,
        )

        # print values in Angstrom
        from kallisto.units import Bohr

        verbosePrinter(
            config.verbose,
            "L, Bmin, Bmax / A: {:5.2f} {:5.2f} {:5.2f}".format(
                L * Bohr, bmin * Bohr, bmax * Bohr
            ),
            out,
        )

    return L, bmin, bmax
