# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os

from monty.io import zopen
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.analysis.transition_state import NEBAnalysis
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pybat.core import Cathode

"""
Utility commands for the pybat package.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "May 2018"


def show_path(directory, filename):
    """
    Show the final migration for a NEB calculation.

    Returns:

    """
    # TODO This is quite inefficient, since the NEBAnalysis script also
    # parses the OUTCAR files. However, this was the fastest solution
    # implementation wise. Improve when feeling less lazy.
    neb = NEBAnalysis.from_dir(directory)

    transition_structure = neb.structures[0].copy()
    for structure in neb.structures[1:]:
        for site in structure:
            transition_structure.append(site.specie, site.frac_coords)

    transition_structure.to("cif", filename)


def conventional_structure(structure_file, fmt="json"):
    """

    Args:
        structure_file:
        fmt:

    Returns:

    """
    cathode = Cathode.from_file(structure_file)
    spg = SpacegroupAnalyzer(cathode)

    conv_structure_file = structure_file.split(".")[0] + "_conv" + "." + fmt
    spg.get_conventional_standard_structure().to(fmt, conv_structure_file)


def primitive_structure(structure_file, fmt="json"):
    """

    Args:
        structure_file:
        fmt:

    Returns:

    """
    cathode = Cathode.from_file(structure_file)
    spg = SpacegroupAnalyzer(cathode)

    prim_structure_file = structure_file.split(".")[0] + "_conv" + "." + fmt
    spg.get_primitive_standard_structure().to(fmt, prim_structure_file)


def make_supercell(structure_file, supercell, fmt="json"):
    """
    Make a supercell of the Structure in the structure file.

    #TODO Assumes the new structure must be of the same
    class as the original, i.e. a Cathode is transformed into another cathode.

    Args:
        structure_file:
        supercell:
        fmt:

    Returns:

    """

    # Turn into list TODO add checks
    supercell_list = [int(number) for number in supercell]

    # Load the structure as a Cathode
    cathode = Cathode.from_file(structure_file)
    cathode.make_supercell(supercell_list)

    super_structure_file = structure_file.split(".")[0] + "_" + supercell \
                           + "." + fmt

    cathode.to(fmt, super_structure_file)


def data(vasprun_file):
    """
    Extract the main data from a vasprun.xml file and write it as a data.json.

    Args:
        vasprun_file (str):

    Returns:

    """
    directory = os.path.dirname(os.path.abspath(vasprun_file))
    vasprun = Vasprun(vasprun_file)

    with zopen(os.path.join(directory, "data.json"), "w") as file:
        file.write(vasprun.to_json())


def print_structure(structure_file):
    print(Cathode.from_file(structure_file))
