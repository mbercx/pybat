# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

from pybat.core import Cathode
from pymatgen.core import Structure

"""
Set of scripts used to define structural changes easily using the command line
interface.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "May 2018"


def define_migration(structure_file, write_cif=False):
    """
    This script allows the user to define a migration in a structure.

    The user has to identify the site that is migrating, as well as provide the
    coordinates to which the site will migrate.

    """
    initial_structure = Structure.from_file(structure_file)
    final_structure = initial_structure.copy()

    # Ask the user for the migration site
    print("")
    print(initial_structure)
    print("")
    migration_site_index = int(input("Please provide the index of the "
                                     "migrating cation:\n"))
    migration_species = initial_structure.sites[migration_site_index].specie

    # Ask the user for the final coordinates of the migrating ion
    final_coords = input("Please provide the final fractional coordinates "
                         "of the migrating site:\n")
    final_coords = [float(number) for number
                    in list(final_coords.split(" "))]

    # Remove the original site from the final structure
    final_structure.remove_sites([migration_site_index])

    # Add the final position of the migrating ion
    final_structure.insert(i=migration_site_index,
                           species=migration_species,
                           coords=final_coords,
                           properties=final_structure.sites[
                                migration_site_index].properties)

    # Write out the initial and final structures
    initial_structure.to("json", "init.json")
    final_structure.to("json", "final.json")

    # Write the structures to a cif format if requested
    if write_cif:
        initial_structure.to("cif", "init.cif")
        final_structure.to("cif", "final.cif")


def define_dimer(structure_file, dimer_indices=(0, 0), distance=0,
                 remove_cations=False):
    """
    Define a dimerization of oxygen in a structure.

    Returns:

    """

    cathode = Cathode.from_file(structure_file)

    if dimer_indices == (0, 0):
        print("")
        print("No site indices were given for structure:")
        print("")
        print(cathode)
        print("")
        dimer_indices = input("Please provide the two indices of the elements "
                              "that need to form a dimer, separated by a "
                              "space: \n")

        dimer_indices = tuple([int(number) for number in
                               list(dimer_indices.split(" "))])

    if distance == 0:
        print("")
        distance = input("Please provide the final distance between the atoms "
                         "in the dimer: \n")
        distance = float(distance)

    if remove_cations:
        # Remove the cations around the oxygen dimer
        cathode.remove_dimer_cations(dimer_indices)

    dimer_structure_file = structure_file.split(".")[0] + "_dimer_init" \
        + ".json"
    cathode.to("json", dimer_structure_file)

    cathode.change_site_distance(dimer_indices, distance)

    dimer_structure_file = structure_file.split(".")[0] + "_dimer_final" \
        + ".json"
    cathode.to("json", dimer_structure_file)
