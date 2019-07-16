# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os

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


def define_migration(structure, site, final_site, write_cif=False):
    """
    This script allows the user to define a migration of an element in a
    Cathode structure. The user has to identify the sites of the migrating element,
    as well as the site the element is migrating to. Note that the final site must be
    empty!

    Args:
        structure (pybat.core.Cathode): Cathode for which to define a migration.
        site: Site or site index of the migrating element.
        final_site: Site or site index of the site the element is migrating to.
        write_cif (bool): Flag that determines if the initial and final
        structures should also be written in a cif format.

    Returns:
        migration_dir (str): The absolute path to the migration directory.

    """
    final_structure = structure.copy()
    final_structure.migrate_element(site=site, final_site=final_site)

    # Set up the migration directory
    site_index = site if isinstance(site, int) else final_structure.index(site)
    final_site_index = final_site if isinstance(final_site, int) \
        else final_structure.index(final_site)
    migration_id = str(site_index) + "_" + str(final_site_index)

    migration_dir = os.path.join(os.getcwd(), "migration_" + migration_id)
    try:
        os.makedirs(migration_dir)
    except FileExistsError:
        print("WARNING: Migration directory already exists.")

    # Set up the filenames
    struc_name = str(structure.composition.reduced_composition).replace(" ", "")
    migration_name = "_m_" + migration_id

    initial_structure_file = struc_name + migration_name + "_init"
    dimer_structure_file = struc_name + migration_name + "_final"

    # Write out the initial and final structures
    structure.to("json", os.path.join(migration_dir, initial_structure_file + ".json"))
    final_structure.to(
        "json", os.path.join(migration_dir, dimer_structure_file + ".json")
    )

    # Write the structures to a cif format if requested
    if write_cif:
        structure.to("cif", os.path.join(migration_dir, initial_structure_file + ".cif"))
        final_structure.to(
            "cif", os.path.join(migration_dir, dimer_structure_file + ".cif")
        )

    return migration_dir


def define_dimer(structure, directory, dimer_indices=(0, 0), distance=0, write_cif=False):
    """
    Define the formation of an oxygen dimer in a Cathode structure.

    The user has to provide the indices of the oxygen pair that will form
    the dimer, as well as the final distance between the oxygen atoms.

    Args:
        structure (pybat.core.LiRichCathode): LiRichCathode for which to define a dimer.
        directory (str): Path to the directory in which the dimer structure files
            should be written.
        dimer_indices (tuple): Indices of the oxygen sites which are to form a
            dimer.
        distance (float): Final distance between the oxygen atoms.
        write_cif (bool): Flag that indicates that the initial and final
        structure files should also be written in a cif format.

    Returns:
        dimer_dir (str): Absolute path to the dimer directory.

    """

    if dimer_indices == (0, 0):
        print("")
        print("No site indices were given for structure:")
        print("")
        print(structure)
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

    # Change the distance between the oxygen atoms for the dimer structure
    dimer_structure = structure.copy()
    dimer_structure.change_site_distance(dimer_indices, distance)

    # Set up the dimer directory
    try:
        os.makedirs(directory)
    except FileExistsError:
        pass

    # Set up the filenames
    struc_name = str(structure.composition.reduced_composition).replace(" ", "")
    dimer_name = "_d_" + "_".join([str(el) for el in dimer_indices])

    initial_structure_file = struc_name + dimer_name + "_init"
    dimer_structure_file = struc_name + dimer_name + "_final"

    # Write out the initial and final structures
    structure.to("json", os.path.join(directory, initial_structure_file + ".json"))
    dimer_structure.to("json", os.path.join(directory, dimer_structure_file + ".json"))

    # Write the structures to a cif format if requested
    if write_cif:
        structure.to("cif", os.path.join(directory, initial_structure_file + ".cif"))
        dimer_structure.to("cif", os.path.join(directory, dimer_structure_file + ".cif"))
