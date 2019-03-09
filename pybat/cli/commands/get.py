# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import pdb

from pybat.core import Cathode, DimerNEBAnalysis
from pymatgen import Structure
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.analysis.transition_state import NEBAnalysis

"""
Set of scripts used to extract information from VASP output files for analysis.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "pre-alpha"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Mar 2019"


def get_structure(directory, write_cif=False):
    """
    Construct a .json file with the structure and magnetic moment from the
    output of a VASP calculation, i.e. the CONTCAR and OUTCAR file.

    Args:
        directory (str): Directory in which the geometry optimization
                output files (i.e. CONTCAR and OUTCAR) are stored.
        write_cif (bool): Flag that indicates whether the structure should
            also be written as a .cif file.
    """
    directory = os.path.abspath(directory)
    structure = Structure.from_file(os.path.join(directory, "CONTCAR"))
    out = Outcar(os.path.join(directory, "OUTCAR"))

    magmom = [site["tot"] for site in out.magnetization]

    # Add the magnetic moments to the Structure
    try:
        structure.add_site_property("magmom", magmom)
    except ValueError:
        # If something goes wrong in assigning the magnetic moments,
        # give the user a warning and assign magnetic moment zero to all sites.
        print("WARNING: Could not assign the magnetic moments found in the "
              "OUTCAR file. They may be missing.")
        structure.add_site_property("magmom", len(structure.sites) * [0])

    structure.to("json", "structure.json")

    if write_cif:
        structure.to("cif", "structure.cif")


def get_cathode(directory, to_current_dir=False, write_cif=False,
                ignore_magmom=False):
    """
    Construct a .json file of the updated Cathode from a geometry
    optimization, based on the initial_cathode.json file and the output of a
    VASP calculation, i.e. the CONTCAR and OUTCAR files. All these files must
    be present in the directory.

    Args:
        directory (str): Directory in which the geometry optimization
            calculation was performed. Must contain the initial_cathode.json,
            OUTCAR and CONTCAR file.
        to_current_dir (bool): Write the output final_cathode files to the
            current working directory.
        write_cif (bool): Flag that determines whether a .cif file of the
            cathode structure is written to the directory.
        ignore_magmom (bool): Flag that indicates that the final magnetic
            moments of the optimized structure should be ignored. This means
            that the magnetic moments of the initial structure will be used.

    Returns:
        None

    """
    directory = os.path.abspath(directory)
    cathode = Cathode.from_file(os.path.join(directory,
                                             "initial_cathode.json"))
    cathode.update_sites(directory, ignore_magmom=ignore_magmom)

    if to_current_dir:
        filename = os.path.join(os.getcwd(), "final_cathode")
    else:
        filename = os.path.join(directory, "final_cathode")

    cathode.to("json", filename + ".json")

    if write_cif:
        cathode.to("cif", filename + ".cif")


def get_barrier(directory, method="pymatgen"):
    """
    Plot the migration barrier of a transition in a directory.
    Args:
        directory (str):
        method (str):

    Returns:

    """
    if method == "pymatgen":
        # The pymatgen.analysis.transition_state module has an object that
        # allows you to

        neb = NEBAnalysis.from_dir(directory, relaxation_dirs=('initial',
                                                               'final'))
        neb.get_plot().show()

    if method == "dimers":
        # This method makes some assumptions about the directory structure
        # for it to work:
        #
        # - The image directories are two characters long, and there are no
        # other directories which are two characters long.
        # - The directory in which the nudged elastic band was performed
        #  contains the dimer indices, delimited by '_', and with no other
        # numbers delimited in such a way present.

        if os.path.exists(os.path.join(directory, "neb_data.json")):
            neb = DimerNEBAnalysis.from_file(
                os.path.join(directory, "neb_data.json")
            )
        else:
            neb = DimerNEBAnalysis.from_dir(directory)
            neb.to("json", os.path.join(directory, "neb_data.json"))

        neb.setup_spline({"saddle_point": "zero_slope"})
        neb.get_plot(label_barrier=False).show()


def get_voltage(directory, calculation="relax", functional=None):
    """
    Calculate the voltage of a battery consisting of a cathode specified by the
    directory versus a metallic Li anode.

    Args:
        directory:
        calculation:
        functional:
    """
    raise NotImplementedError


def get_endiff(directory):
    """
    Calculate the energy difference for a transition in a directory.

    Args:
        directory:

    Returns:

    """
    initial_outcar = Outcar(os.path.join(directory, "initial", "OUTCAR"))
    final_outcar = Outcar(os.path.join(directory, "final", "OUTCAR"))

    initial_energy = initial_outcar.final_energy
    final_energy = final_outcar.final_energy

    print("The energy difference is: ", end="")
    print(str(final_energy - initial_energy) + " eV")


# SO plagiarism

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
