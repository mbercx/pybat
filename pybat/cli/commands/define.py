
from pybat.core import Cathode, change_site_distance
from pymatgen.core import Structure

"""
Set of scripts used to define state transitions easily from the command line
interface.

"""



def define_migration(structure_file, provide_coords=False, write_cif=False):
    """
    This script allows the user to define a migration in a structure.

    As a first implementation, the script will simply start from the fully
    lithiated cathode structure and then the user can specify with cation
    site should migrate to which other one. The script then removes the atom
    which occupies the site to which the cation migrates. After this, the user
    can still remove other atoms as desired.

    A second addition is the possibility of choosing the final migration
    coordinates, in case the final site is not among the sites of the
    structure.

    """
    final_structure = Structure.from_file(structure_file)

    # Ask the user for the migration sites and the other sites which are to be
    # removed
    print("")
    print(final_structure)
    print("")
    migration_site_index = int(input("Please provide the index of the "
                                     "migrating cation:\n"))
    if provide_coords:
        final_coords = input("Please provide the final fractional coordinates "
                             "of the migrating site:\n")
        final_coords = [float(number) for number
                        in list(final_coords.split(" "))]
    else:
        print("")
        final_site_index = int(input("Please provide the index of the site to "
                                     "which the cation is migrating:\n"))
    print("")
    remove_site_indices = input("Provide addition site indices to remove:\n")

    # A bit of processing
    if remove_site_indices == "":
        remove_site_indices = []
    else:
        remove_site_indices = [int(number) for number in
                               list(remove_site_indices.split(" "))]

    migration_species = final_structure.sites[migration_site_index].specie

    initial_structure = final_structure.copy()
    if provide_coords:
        # Add the final position of the migrating ion
        final_structure.append(migration_species,
                               final_coords,
                               properties=final_structure.sites[
                                    migration_site_index].properties)
    else:
        # Replace the final migration site by the migrating species
        final_structure.replace(final_site_index,
                                migration_species,
                                properties=final_structure.sites[
                                    migration_site_index].properties)

    # Remove the original site from the final structure, as well as the other
    # sites which are to be removed.
    final_structure.remove_sites([migration_site_index] + remove_site_indices)

    if provide_coords:
        # Remove the other sites as determined by the user
        initial_structure.remove_sites(remove_site_indices)
    else:
        # Remove the final site from the initial structure, as well as the
        # other sites which are to be removed.
        initial_structure.remove_sites([final_site_index] +
                                       remove_site_indices)

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
    cathode.to(dimer_structure_file)

    cathode.change_site_distance(dimer_indices, distance)

    dimer_structure_file = structure_file.split(".")[0] + "_dimer_final" \
                           + ".json"
    cathode.to("json", dimer_structure_file)

