# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import _warnings as warnings

from string import ascii_letters
from pybat.core import Cathode
from pymatgen.core import Composition

"""
Set of scripts used to define structural changes easily using the command line
interface.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "pre-alpha"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Mar 2019"


def define_migration(structure_file, migration_indices=(0, 0),
                     write_cif=False):
    """
    This script allows the user to define a migration of a cation in a
    Cathode structure.

    The user has to identify the site that is migrating, as well as provide the
    fractional coordinates to which the site will migrate, or a vacant site
    index.

    Args:
        structure_file (str): Path to the structure file.
        migration_indices (tuple): Tuple of the indices which designate the
        migrating site and the vacant site to which the cation will migrate.
        write_cif (bool): Flag that determines if the initial and final
        structures should also be written in a cif format.

    Returns:
        migration_dir (str): The absolute path to the migration directory.

    """
    cathode = Cathode.from_file(structure_file)
    final_structure = cathode.copy()

    if migration_indices == (0, 0):
        # Prompt the user for the migration site
        print("")
        print(cathode)
        print("")
        migration_site_index = int(input("Please provide the index of the "
                                         "migrating cation (Note: Not the "
                                         "VESTA index!): "))
        print("")

        migration_site = cathode.sites[migration_site_index]
        migration_species = migration_site.species_and_occu

        # Check if the site to which the ion should migrate is actually
        # occupied
        if migration_species == Composition():
            raise ValueError("Chosen site is vacant.")

        # Ask the user for the final coordinates of the migrating ion
        final_coords = input("Please provide the index of the site the cation "
                             "is migrating to, or the final fractional "
                             "coordinates of the migration site: ")
        print("")
        final_coords = [float(number) for number
                        in list(final_coords.split(" "))]
    else:
        migration_site_index = migration_indices[0]
        migration_site = cathode.sites[migration_site_index]
        migration_species = migration_site.species_and_occu

        # Check if the site to which the ion should migrate is actually
        # occupied
        if migration_species == Composition():
            raise ValueError("Chosen site is vacant.")
        final_coords = [migration_indices[1]]

    # In case of a site index as input
    if len(final_coords) == 1:

        # Grab the required information about the final site
        final_site_index = int(final_coords[0])
        final_site = cathode.sites[final_site_index]
        final_coords = final_site.frac_coords
        final_species = final_site.species_and_occu

        # Check if site is occupied
        if final_species != Composition():
            raise ValueError("Chosen final site is not vacant.")

        # Change the coordinates of the migration site with the ones of
        # the final site.
        final_structure.replace(i=migration_site_index,
                                species=migration_species,
                                coords=final_coords,
                                properties=migration_site.properties)

        # Do the opposite for the final site
        final_structure.replace(i=final_site_index,
                                species=final_species,
                                coords=migration_site.frac_coords,
                                properties=final_site.properties)

        migration_id = str(migration_site_index) + "_" + str(final_site_index)

    # In case of a set of fractional coordinates as input
    elif len(final_coords) == 3:

        # Replace the site with the site of the new coordinates
        final_structure.replace(i=migration_site_index,
                                species=migration_species,
                                coords=final_coords,
                                properties=final_structure.sites[
                                    migration_site_index].properties)

        migration_id = str(migration_site_index) + "_a"

        letter_index = 0

        while "migration_" + migration_id in os.listdir(os.getcwd()) and \
                        letter_index < len(ascii_letters):
            letter_index += 1
            migration_id = str(migration_site_index) + "_" + \
                           ascii_letters[letter_index]

    else:
        raise IOError("Provided input is incorrect.")

    # Replace the
    final_structure.remove_sites([migration_site_index])

    # Add the final position of the migrating ion
    final_structure.insert(i=migration_site_index,
                           species=migration_species,
                           coords=final_coords,
                           properties=final_structure.sites[
                               migration_site_index].properties)

    # Set up the migration directory
    migration_dir = os.path.join(os.getcwd(), "migration_" + migration_id)
    try:
        os.mkdir(migration_dir)
    except FileExistsError:
        print("WARNING: Migration directory already exists.")

    # Set up the filenames
    initial_structure_file = ".".join(structure_file.split("/")[-1].split(".")[
                                      0:-1]) + "_m_" + migration_id + "_init"
    final_structure_file = ".".join(structure_file.split("/")[-1].split(".")[
                                    0:-1]) + "_m_" + migration_id + "_final"

    # Write out the initial and final structures
    cathode.to("json", migration_dir + "/" + initial_structure_file + ".json")
    final_structure.to("json",
                       migration_dir + "/" + final_structure_file + ".json")

    # Write the structures to a cif format if requested
    if write_cif:
        cathode.to("cif",
                   migration_dir + "/" + initial_structure_file + ".cif")
        final_structure.to("cif",
                           migration_dir + "/" + final_structure_file + ".cif")

    return migration_dir


def define_dimer(structure_file, dimer_indices=(0, 0), distance=0,
                 remove_cations=False, write_cif=False):
    """
    Define the formation of an oxygen dimer in a Cathode structure.

    The user has to provide the indices of the oxygen pair that will form
    the dimer, as well as the final distance between the oxygen atoms.

    Args:
        structure_file (str): Path to the Cathode structure file.
        dimer_indices (tuple): Indices of the oxygen sites which are to form a
        dimer.
        distance (float): Final distance between the oxygen atoms.
        remove_cations (bool): Flag that allows the user to remove the
        cations (Li, Na, ...) around the chosen oxygen pair.
        write_cif (bool): Flag that indicates that the initial and final
        structure files should also be written in a cif format.

    Returns:
        dimer_dir (str): Absolute path to the dimer directory.

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
                              "space (Note: Not the VESTA indices!): ")

        dimer_indices = tuple([int(number) for number in
                               list(dimer_indices.split(" "))])

    if distance == 0:
        distance = input("Please provide the final distance between the atoms "
                         "in the dimer: ")
        print("")
        distance = float(distance)

    if remove_cations:
        # Remove the cations around the oxygen dimer
        cathode.remove_dimer_cations(dimer_indices)

    # Change the distance between the oxygen atoms for the dimer structure
    dimer_structure = cathode.copy()
    dimer_structure.change_site_distance(dimer_indices, distance)

    # Create the dimer directory
    dimer_dir = os.path.join(
        os.getcwd(), "dimer_" + "_".join([str(el) for el in dimer_indices])
    )
    try:
        os.mkdir(dimer_dir)
    except FileExistsError:
        warnings.warn("Warning: " + dimer_dir + " already exists, "
                                                "overwriting...")

    # Set up the filenames
    initial_structure_file = ".".join(
        structure_file.split("/")[-1].split(".")[0:-1]) + "_d_" \
                             + "_".join([str(el) for el in dimer_indices]) \
                             + "_init"
    dimer_structure_file = ".".join(
        structure_file.split("/")[-1].split(".")[0:-1]) + "_d_" \
                           + "_".join([str(el) for el in dimer_indices]) \
                           + "_final"

    # Write out the initial and final structures
    cathode.to("json", dimer_dir + "/" + initial_structure_file + ".json")
    dimer_structure.to("json",
                       dimer_dir + "/" + dimer_structure_file + ".json")

    # Write the structures to a cif format if requested
    if write_cif:
        cathode.to("cif", dimer_dir + "/" + initial_structure_file + ".cif")
        dimer_structure.to("cif",
                           dimer_dir + "/" + dimer_structure_file + ".cif")

    return dimer_dir
