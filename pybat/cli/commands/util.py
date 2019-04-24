# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os

from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pybat.core import Cathode

"""
Utility commands for the pybat package.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "pre-alpha"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Mar 2019"


def show_path(directory, filename):
    """
    Show the final migration for a NEB calculation.

    Returns:

    """
    initial_structure = Cathode.from_file(
        os.path.join(directory, "initial", "final_cathode.json")).as_ordered_structure()
    final_structure = Cathode.from_file(
        os.path.join(directory, "final", "final_cathode.json")).as_ordered_structure()

    image_dirs = [d for d in os.listdir(directory) if "0" in d and os.path.isdir(d)]
    image_dirs.sort()
    image_structures = [Structure.from_file(os.path.join(directory, idir, "CONTCAR"))
                        for idir in image_dirs[1:-1]]

    image_structures.append(final_structure)

    transition_structure = initial_structure.copy()
    for structure in image_structures:
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


def print_structure(structure_file):
    print(Cathode.from_file(structure_file))
