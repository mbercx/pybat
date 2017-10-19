
from pybat.core import change_site_distance
from pymatgen.core import Structure

"""
Set of scripts used to define state transitions easily from the command line
interface.

"""



def define_migration():
    """
    This script allows the user to define a migration in a structure.

    """
    pass


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

