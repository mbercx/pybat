import numpy as np
import os
import shutil

from pybat.core import LiRichCathode
from pybat.sets import PybatRelaxSet, PybatNEBSet

from pymatgen.core import Structure
from pymatgen.analysis.path_finder import ChgcarPotential, NEBPathfinder
from pymatgen.io.vasp.outputs import Chgcar, Outcar
from pymatgen.io.vasp.sets import MPRelaxSet, MPHSERelaxSet, MPStaticSet

"""
Setup scripts for the calculations.
"""

DFT_FUNCTIONAL = "PBE_54"


def relax(structure_file, calculation_dir, hse_calculation=False):
    """
    Set up a standard relaxation of a structure.

    """
    user_incar_settings = {"ISMEAR": 0}

    structure_file = os.path.abspath(structure_file)
    structure = Structure.from_file(structure_file)
    calculation_dir = os.path.abspath(calculation_dir)

    # Check if a magnetic moment was not provided for the sites. If not, make
    # sure it is zero for the calculations.
    if "magmom" not in structure.site_properties.keys():
        structure.add_site_property("magmom", [0] * len(structure.sites))

    if hse_calculation:
        geo_optimization = MPHSERelaxSet(structure=structure,
                                         potcar_functional=DFT_FUNCTIONAL,
                                         user_incar_settings={"EDIFFG": -0.01})

    else:
        geo_optimization = MPRelaxSet(structure=structure,
                                      potcar_functional=DFT_FUNCTIONAL,
                                      user_incar_settings=user_incar_settings)

    geo_optimization.write_input(calculation_dir)


def transition(directory, initial_structure, final_structure,
               is_migration=False, hse_calculation=False):
    """
    This script will set up the geometry optimizations for the initial and
    final structures.

    If requested, it will also set up a charge density calculation for the
    "host structure", i.e. the structure with vacancies at the initial and
    final locations of the migrating ion.

    """

    # Check if a magnetic moment was not provided for the sites. If not, make
    # sure it is zero for the calculations.
    if "magmom" not in initial_structure.site_properties.keys():
        initial_structure.add_site_property("magmom",
                                            [0]*len(initial_structure.sites))

    if "magmom" not in final_structure.site_properties.keys():
        final_structure.add_site_property("magmom",
                                          [0] * len(initial_structure.sites))

    # Set up the initial and final optimization calculations
    initial_optimization = PybatRelaxSet(structure=initial_structure,
                                         potcar_functional=DFT_FUNCTIONAL,
                                         hse_calculation=hse_calculation)

    final_optimization = PybatRelaxSet(structure=final_structure,
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

        migration_site_index = find_migrating_ion(initial_structure,
                                                  final_structure)

        host_structure = initial_structure.copy()
        host_structure.remove_sites([migration_site_index])
        host_scf = MPStaticSet(host_structure,
                               potcar_functional=DFT_FUNCTIONAL)
        host_scf.write_input(os.path.join(neb_dir, "host"))


def dimers(structure_file, dimer_distance=1.4, hse_calculation=False):
    """

    Args:
        structure_file:

    Returns:

    """

    # Load the cathode structure
    cathode = LiRichCathode.from_file(structure_file)

    # Find the non-equivalent dimers
    dimers = cathode.find_noneq_dimers()

    # Set up the geometry optimization calculations for the various dimers
    for dimer in dimers:

        # Set up the dimer directory
        dimer_directory = "".join(str(dimer.indices[0]) + "_"
                                  + str(dimer.indices[1]))

        os.mkdir(dimer_directory)
        initial_dir = os.path.join(dimer_directory, "initial")
        os.mkdir(initial_dir)
        final_dir = os.path.join(dimer_directory, "final")
        os.mkdir(final_dir)

        # Copy the structure into the 'initial' directory
        shutil.copy(structure_file, os.path.join(initial_dir,
                                                 "structure.json"))

        # Set up the dimer structure, i.e. move the oxygen pair closer together
        dimer_structure = cathode.copy()
        dimer_structure.change_site_distance(site_indices=dimer.indices,
                                             distance=dimer_distance)

        # Set up the geometry optimization for the dimer structure
        dimer_optimization = PybatRelaxSet(structure=dimer_structure,
                                           potcar_functional=DFT_FUNCTIONAL,
                                           hse_calculation=hse_calculation)

        # Write the calculation files to the 'final' directory
        dimer_optimization.write_input(final_dir)


def neb(directory, nimages=8, is_migration=False,
        hse_calculation=False):
    """
    Set up the NEB calculation from the initial and final structures.

    """
    directory = os.path.abspath(directory)

    # Extract the optimized initial and final geometries
    initial_dir = os.path.join(directory, "initial")
    final_dir = os.path.join(directory, "final")

    try:
        initial_structure = Structure.from_file(os.path.join(initial_dir,
                                                             "CONTCAR"))

        # Add the magnetic configuration to the initial structure
        initial_out = Outcar(os.path.join(initial_dir, "OUTCAR"))
        initial_magmom = [site["tot"] for site in initial_out.magnetization]
        initial_structure.add_site_property("magmom", initial_magmom)

    except FileNotFoundError:
        # In case the required output files are not present, check to see if
        # the structure is present in a json format
        initial_structure = Structure.from_file(os.path.join(initial_dir,
                                                             "structure.json"))

    final_structure = Structure.from_file(os.path.join(final_dir,
                                                       "CONTCAR"))

    # In case the transition is a migration
    if is_migration:
        # Set up the static potential for the Pathfinder from the host charge
        # density
        host_charge_density = Chgcar.from_file(os.path.join(directory, "host"))
        host_potential = ChgcarPotential(host_charge_density)

        migration_site_index = find_migrating_ion(initial_structure,
                                                  final_structure)

        neb = NEBPathfinder(start_struct=initial_structure,
                            end_struct=final_structure,
                            relax_sites=migration_site_index,
                            v=host_potential)

        images = neb.images
        neb.plot_images("neb.vasp")

    else:
        # Linearly interpolate the initial and final structures
        images = initial_structure.interpolate(end_structure=final_structure,
                                               nimages=nimages)

    neb_calculation = PybatNEBSet(images, potcar_functional=DFT_FUNCTIONAL)

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
        initial_contains:
        final_contains:

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

    return initial_structure, final_structure


def find_migrating_ion(initial_structure, final_structure):
    """
    Tries to find the index of the ion in the structure that is migrating, by
    considering the distance that the sites have moved. The site whose
    coordinates have changed the most is considered to be the migrating ion.

    Note that the site indices of the initial and final structure have to
    correspond for the algorithm to work.

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

        # TODO This definition of distance is inadequate, i.e. if the site has crossed into another unit cell, the distance will be very big. This can be fixed by using the nearest image!
        distance = np.linalg.norm(site_pair[0].coords - site_pair[1].coords)

        if distance > max_distance:
            max_distance = distance
            migrating_site = site_pair[0]

    return initial_structure.sites.index(migrating_site)
