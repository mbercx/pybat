# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os

from pymatgen import Structure
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.analysis.transition_state import NEBAnalysis

"""
Set of scripts used to extract information from VASP output files for analysis.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "May 2018"

# Total Energy per Li of metallic lithium
LI_ENERGY = -1.89


def get_structure(directory, write_cif=False):
    """
    Construct a .json file with the structure and magnetic moment from the
    output of a VASP calculation, i.e. the CONTCAR and OUTCAR file.

    Args:
        directory:

    Returns:

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
        structure.add_site_property("magmom", len(structure.sites)*[0])

    structure.to("json", "structure.json")

    if write_cif:
        structure.to("cif", "structure.cif")


def get_barrier(directory):
    """
    Plot the migration barrier of a transition in a directory.
    Args:
        directory:

    Returns:

    """
    neb = NEBAnalysis.from_dir(directory, relaxation_dirs=('initial',
                                                           'final'))
    neb.get_plot().show()


def get_voltage(directory, calculation="relax", functional=None):
    """
    Calculate the voltage of a battery consisting of a cathode specified by the
    directory versus a metallic Li anode.

    Args:
        directory:
        calculation:
        functional:

    Returns:

    """
    pass


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
    print(str(final_energy-initial_energy) + " eV")