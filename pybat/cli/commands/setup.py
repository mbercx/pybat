import numpy as np
import os

from pybat.sets import pybatRelaxSet, pybatNEBSet

from pymatgen.core import Structure
from pymatgen.analysis.path_finder import ChgcarPotential, NEBPathfinder
from pymatgen.io.vasp.outputs import Chgcar, Outcar
from pymatgen.io.vasp.sets import MPRelaxSet, MPHSERelaxSet, MPStaticSet

"""
Setup scripts for the calculations.
"""

DFT_FUNCTIONAL = "PBE_54"

def set_up_relaxation(structure_file, calculation_dir, hse_calculation=False):
    """
    Set up a standard relaxation of a structure.

    """
    USER_INCAR_SETTINGS = {"ISMEAR": 0}

    structure_file = os.path.abspath(structure_file)
    structure = Structure.from_file(structure_file)
    calculation_dir = os.path.abspath(calculation_dir)

    # Check if a magnetic moment was not provided for the sites. If not, make
    # sure it is zero for the calculations.
    if not "magmom" in structure.site_properties.keys():
        structure.add_site_property("magmom", [0] * len(structure.sites))

    if hse_calculation:
        geo_optimization = MPHSERelaxSet(structure=structure,
                                         potcar_functional=DFT_FUNCTIONAL,
                                         user_incar_settings={"EDIFFG": -0.01})

    else:
        geo_optimization = MPRelaxSet(structure=structure,
                                      potcar_functional=DFT_FUNCTIONAL,
                                      user_incar_settings=USER_INCAR_SETTINGS)

    geo_optimization.write_input(calculation_dir)


def set_up_transition(directory, initial_structure, final_structure,
                      is_migration=False, hse_calculation=False):
    """
    This script will set up the geometry optimizations for the initial and final
    structures.

    If requested, it will also set up a charge density calculation for the
    "host structure", i.e. the structure with vacancies at the initial and
    final locations of the migrating ion.

    """

    # Check if a magnetic moment was not provided for the sites. If not, make
    # sure it is zero for the calculations.
    if not "magmom" in initial_structure.site_properties.keys():
        initial_structure.add_site_property("magmom",
                                            [0]*len(initial_structure.sites))

    if not "magmom" in final_structure.site_properties.keys():
        final_structure.add_site_property("magmom",
                                          [0] * len(initial_structure.sites))

    # Set up the initial and final optimization calculations
    initial_optimization = pybatRelaxSet(structure=initial_structure,
                                         potcar_functional=DFT_FUNCTIONAL,
                                         hse_calculation=hse_calculation)

    final_optimization = pybatRelaxSet(structure=final_structure,
                                       potcar_functional=DFT_FUNCTIONAL,
                                       hse_calculation=hse_calculation)

    # Set up the root directory for the neb calculation
    neb_dir = os.path.abspath(directory)

    # Write input files to directories
    initial_optimization.write_input(os.path.join(neb_dir, "initial"))
    final_optimization.write_input(os.path.join(neb_dir, "final"))

    # If the transition is a migration of an atom in the structure, set up the
    # calculation for the charge density, used to find a better initial pathway
    if is_migration:

        migration_site_index = initial_structure.sites.index(
            find_migrating_ion(initial_structure, final_structure)
        )

        host_structure = initial_structure.copy()
        host_structure.remove_sites([migration_site_index])
        host_scf = MPStaticSet(host_structure,
                               potcar_functional=DFT_FUNCTIONAL)
        host_scf.write_input(os.path.join(neb_dir, "host"))


def set_up_NEB(directory, nimages=8, is_migration=False,
               hse_calculation=False):
    """
    Set up the NEB calculation from the initial and final structures.

    """
    directory = os.path.abspath(directory)

    # Extract the optimized initial and final geometry
    initial_dir = os.path.join(directory,"initial")
    final_dir = os.path.join(directory,"final")

    initial_structure = Structure.from_file(os.path.join(initial_dir,
                                                         "CONTCAR"))
    final_structure = Structure.from_file(os.path.join(final_dir,
                                                       "CONTCAR"))

    # Add the magnetic configuration to the initial structure
    initial_out = Outcar(os.path.join(initial_dir,"OUTCAR"))
    initial_magmom = [site["tot"] for site in initial_out.magnetization]
    initial_structure.add_site_property("magmom", initial_magmom)

    # In case the transition is a migration
    if is_migration:
        # Set up the static potential for the Pathfinder from the host charge
        # density
        host_charge_density = Chgcar.from_file(os.path.join(directory, "host"))
        host_potential = ChgcarPotential(host_charge_density)

        migration_site_index = initial_structure.sites.index(
            find_migrating_ion(initial_structure, final_structure)
        )

        neb = NEBPathfinder(start_struct=initial_structure,
                            end_struct=final_structure,
                            relax_sites= migration_site_index,
                            v=host_potential)

        images = neb.images
        neb.plot_images("neb.vasp")

    else:
        # Linearly interpolate the initial and final structures
        images = initial_structure.interpolate(end_structure=final_structure,
                                               nimages=nimages)


    neb_calculation = pybatNEBSet(images, potcar_functional=DFT_FUNCTIONAL)

    # Set up the NEB calculation
    neb_calculation.write_input(directory)

    # Make a file to visualize the transition
    neb_calculation.visualize_transition()

###########
# UTILITY #
###########

def find_transition_structures(directory, initial_contains="init",
                               final_contains="final"):
    """
    Find the initial and final structures for a transition from the files in a
    directory.

    Args:
        directory:

    Returns:

    """
    directory = os.path.abspath(directory)

    initial_structure_file = None
    final_structure_file = None

    for item in os.listdir(directory):

        if initial_contains in item and os.path.isfile(item):
            initial_structure_file = os.path.join(directory, item)

        if final_contains in item and os.path.isfile(item):
            final_structure_file = os.path.join(directory, item)

    if initial_structure_file:
        initial_structure = Structure.from_file(initial_structure_file)
    else:
        raise FileNotFoundError("No suitably named initial structure file in "
                                "directory.")

    if final_structure_file:
        final_structure = Structure.from_file(final_structure_file)
    else:
        raise FileNotFoundError("No suitably named final structure file in "
                                "directory.")

    return (initial_structure, final_structure)

def find_migrating_ion(initial_structure, final_structure):
    """
    Tries to find the index of the ion in the structure that is migrating, by
    considering the distance that the sites have moved. The site whose
    coordinates have changed the most is considered to be the migrating ion.

    Args:
        initial_structure:
        final_structure:

    Returns:
        (int) Index of the migrating site

    """

    # Check if both structures have the same number of sites.
    if len(initial_structure.sites) != len(final_structure.sites):
        raise IOError("Provided structures do not have the same number of "
                      "atoms.")

    max_distance = 0
    migrating_site = None

    # TODO Build in some checks, i.e. make sure that the other ions have not moved significantly, and that the migrating ion has moved sufficiently.

    for site_pair in zip(initial_structure.sites, final_structure.sites):

        distance = np.linalg.norm(site_pair[0].coords - site_pair[1].coords)

        if distance > max_distance:
            max_distance = distance
            migrating_site = site_pair[0]

    return migrating_site







