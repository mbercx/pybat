
from pybat.core import change_site_distance
from pymatgen.core import Structure

"""
Set of scripts used to define state transitions easily from the command line
interface.

"""



def define_migration(structure_file, write_cif=False):
    """
    This script allows the user to define a migration in a structure.

    As a first implementation, the script will simply start from the fully
    lithiated cathode structure and then the user can specify with cation
    site should migrate to which other one. The script then removes the atom
    which occupies the site to which the cation migrates. After this, the user
    can still remove other atoms as desired.

    """
    final_structure = Structure.from_file(structure_file)

    # Ask the user for the migration sites and the other sites which are to be
    # removed
    print("")
    print(final_structure)
    print("")
    migration_site_index = int(input("Please provide the index of the "
                                     "migrating cation:\n"))
    print("")
    final_site_index = int(input("Please provide the index of the site to "
                                 "which the cation is migrating:\n"))
    print("")
    remove_site_indices = input("Provide addition site indices to remove:\n")

    remove_site_indices = tuple([int(number) for number in
                                 list(remove_site_indices.split(" "))])

    initial_structure = final_structure.copy()
    # Replace the final migration site by the migrating species
    final_structure.replace(final_site_index,
                            final_structure.sites[migration_site_index].specie,
                            properties=final_structure.sites[
                                migration_site_index].properties)

    # Remove the original site from the final structure, as well as the other
    # sites which are to be removed.
    final_structure.remove_sites([migration_site_index] + remove_site_indices)

    # Remove the other sites from the initial structure as well
    initial_structure.remove_sites(remove_site_indices)

    # Write out the initial and final structures
    initial_structure.to("json", "init.json")
    final_structure.to("json", "final.json")

    # Write the structures to a cif format if requested
    if write_cif:
        initial_structure.to("cif", "init.cif")
        final_structure.to("cif", "final.cif")


def define_dimer(structure_file, site_indices=(0,0), distance=0):
    """
    Define a dimerization of oxygen in a structure.

    Returns:

    """

    structure = Structure.from_file(structure_file)

    if site_indices == (0,0):
        print("")
        print("No site indices were given for structure:")
        print("")
        print(structure)
        print("")
        site_indices = input("Please provide the two indices of the elements "
                             "that need to form a dimer, separated by a "
                             "space: \n")

        site_indices = tuple([int(number) for number in
                              list(site_indices.split(" "))])

    if distance == 0:
        print("")
        distance = input("Please provide the final distance between the atoms "
                         "in the dimer: \n")
        distance = float(distance)

    change_site_distance(structure, site_indices, distance)

    dimer_structure_file = structure_file.split(".")[0] + "_dimer" + ".json"
    structure.to("json", dimer_structure_file)

