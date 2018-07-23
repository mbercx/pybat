# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os

from pathlib import Path
from ruamel.yaml import YAML

from pymatgen.core import Structure
from pymatgen.analysis.transition_state import NEBAnalysis
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

"""
Utility commands for the pybat package.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "May 2018"

# Load the workflow configuration
CONFIG_FILE = os.path.join(os.path.expanduser("~"), ".pybat_wf_config.yaml")

def workflow_config(settings="all"):
    """
    Script to set up the configuration of the workflow server and jobscripts.

    Returns:
        None

    """
    yaml = YAML()
    yaml.default_flow_style = False

    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, 'r') as config_file:
            config_dict = yaml.load(config_file.read())
    else:
        config_dict = {"SERVER": {}, "WORKFLOW": {}}

    if settings in ["server", "all"]:
        config_dict["SERVER"]["host"] = input("Please provide the server "
                                              "host: ")
        config_dict["SERVER"]["port"] = input("Please provide the port "
                                              "number: ")
        config_dict["SERVER"]["name"] = input("Please provide the server "
                                              "name: ")
        config_dict["SERVER"]["username"] = input("Please provide your "
                                                  "username: ")
        config_dict["SERVER"]["password"] = input("Please provide your "
                                                  "password: ")

    if settings in ["workflow", "all"]:
        script_path = input(
            "Please provide the full path to the workflow script: "
        )
        if not os.path.exists(script_path):
            raise FileNotFoundError("Could not find suggested path.")
        elif not os.path.isabs(script_path):
            print("Provided path is not an absolute path. Finding absolute "
                  "path for proper configuration of the workflows...")
            script_path = os.path.abspath(script_path)

        config_dict["WORKFLOW"]["script_path"] = script_path

    with Path(CONFIG_FILE) as config_file:
        yaml.dump(config_dict, config_file)

def show_path(directory, filename):
    """
    Show the final migration for a NEB calculation.

    Returns:

    """
    # TODO This is quite inefficient, since the NEBAnalysis script also parses the OUTCAR files. However, this was the fastest solution implementation wise.
    neb = NEBAnalysis.from_dir(directory)

    transition_structure = neb.structures[0].copy()
    for structure in neb.structures[1:]:
        for site in structure:
            transition_structure.append(site.specie, site.frac_coords)

    transition_structure.to("cif", filename + ".cif")


def conventional_structure(structure_file, fmt="cif"):
    """

    Args:
        structure_file:
        fmt:

    Returns:

    """
    structure = Structure.from_file(structure_file)
    spg = SpacegroupAnalyzer(structure)

    conv_structure_file = structure_file.split(".")[0] + "_conv" + "." + fmt
    spg.get_conventional_standard_structure().to(fmt, conv_structure_file)


def primitive_structure(structure_file, fmt="cif"):
    """

    Args:
        structure_file:
        fmt:

    Returns:

    """
    structure = Structure.from_file(structure_file)
    spg = SpacegroupAnalyzer(structure)

    prim_structure_file = structure_file.split(".")[0] + "_conv" + "." + fmt
    spg.get_primitive_standard_structure().to(fmt, prim_structure_file)


def make_supercell(structure_file, supercell, fmt="cif"):
    """

    Args:
        structure_file:
        supercell:
        fmt:

    Returns:

    """
    # Turn into list TODO add checks
    supercell_list = [int(number) for number in supercell]

    structure = Structure.from_file(structure_file)
    structure.make_supercell(supercell_list)

    super_structure_file = structure_file.split(".")[0] + "_" + supercell\
        + "." + fmt

    structure.to(fmt, super_structure_file)
