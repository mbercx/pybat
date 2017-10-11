import numpy as np
import os

from pymatgen.core import Structure
from pymatgen.analysis.path_finder import ChgcarPotential, NEBPathfinder
from pymatgen.io.vasp.sets import MPRelaxSet, MPHSERelaxSet, MITNEBSet, \
    MPStaticSet

"""
Setup scripts for the calculations.
"""

DFT_FUNCTIONAL = "PBE_54"
PBE_RELAX_INCAR = {"ISMEAR":0, "EDIFF":1e-4}

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


def define_migration():
    """
    This script allows the user to define a migration in a structure.

    """
    pass

def set_up_transition(directory, initial_structure, final_structure,
                      is_migration=False):
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
    initial_optimization = MPRelaxSet(structure=initial_structure,
                                      potcar_functional=DFT_FUNCTIONAL,
                                      user_incar_settings=PBE_RELAX_INCAR)

    final_optimization = MPRelaxSet(structure=final_structure,
                                    potcar_functional=DFT_FUNCTIONAL,
                                    user_incar_settings=PBE_RELAX_INCAR)

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

        host_structure = initial_structure.remove_sites([migration_site_index])
        host_scf = MPStaticSet(host_structure)
        host_scf.write_input(os.path.join(neb_dir, "host"))


def set_up_NEB(initial_structure, final_structure, charge_density):
    """
    Set up the NEB calculation from the initial and final structure.

    """
    # Set up the static potential for the Pathfinder from the host charge
    # density
    host_potential = ChgcarPotential(charge_density)

    migration_site_index = initial_structure.sites.index(
        find_migrating_ion(initial_structure, final_structure)
    )

    neb = NEBPathfinder(start_struct=initial_structure,
                        end_struct=final_structure,
                        relax_sites= migration_site_index,
                        v=host_potential)

    neb.plot_images("neb.vasp")

    #neb_calculation = MITNEBSet()

# Utility scripts

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







