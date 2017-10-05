# Encoding: utf-8

import pymatgen as pmg
import numpy as np

from pymatgen.core import Structure

"""
Module that contains tools to represent and calculate the properties of
battery cathodes

"""

class Cathode(Structure):
    """
    A class representing a cathode material in a battery.
    """
    def __init__(self, lattice, species, coords, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None):

        super(Cathode,
              self).__init__(lattice=lattice, species=species,
                             coords=coords,
                             validate_proximity=validate_proximity,
                             to_unit_cell=to_unit_cell,
                             coords_are_cartesian=coords_are_cartesian,
                             site_properties=site_properties)

    def remove_cations(self):
        """
        Remove the cations from the cathode, i.e. delithiate the structure in
        case Li is the cation of the cathode.

        :return:
        """
        pass

    def set_to_high_spin(self):
        """

        :return:
        """
        pass

    def set_to_low_spin(self):
        """

        :return:
        """
        pass

class LiRich(Cathode):
    """
    Class representing a Li-rich cathode material.

    """
    pass


def change_site_distance(structure, site_indices, distance):
    """
    Change the coordinates of two sites in a structure in order to adjust
    their distance.

    Args:
        structure: (pymatgen.core.structure.Structure)
        site_indices:

    Returns:
        None
    """

    # TODO Carefully consider the place of this method, i.e. as an independent method or a class method

    site_A = structure.sites[site_indices[0]]
    site_B = structure.sites[site_indices[1]]

    # print("Sites being considered:")
    # print(site_A)
    # print(site_B)

    # Find the distance between the sites, as well as the image of site B that
    # is closest to site A
    (original_distance, closest_image_B) = site_A.distance_and_image(site_B)

    # print("")
    # print("Original Distance = " + str(original_distance))
    # print("Closest image:")
    # print(closest_image_B)

    image_cart_coords = structure.lattice.get_cartesian_coords(
        site_B.frac_coords + closest_image_B
    )

    # print("")
    # print("Image cartesian coordinates:")
    # print(image_cart_coords)

    # Calculate the vector that connects site A with site B
    connection_vector = image_cart_coords - site_A.coords

    # print("")
    # print("Connection vector:")
    # print(connection_vector)
    # print("Connection vector length = " + str(np.linalg.norm(connection_vector)))

    # Make it a unit vector
    connection_vector /= np.linalg.norm(connection_vector)

    # Calculate the distance the sites need to be moved.
    site_move_distance = (original_distance - distance)/2

    # Calculate the new cartesian coordinates of the sites
    new_site_A_coords = site_A.coords + site_move_distance * connection_vector
    new_site_B_coords = site_B.coords - site_move_distance * connection_vector

    # Change the sites in the structure
    structure.replace(i=site_indices[0], species=site_A.species_string,
                      coords=new_site_A_coords,
                      coords_are_cartesian=True, properties=site_A.properties)

    structure.replace(i=site_indices[1], species=site_B.species_string,
                      coords=new_site_B_coords,
                      coords_are_cartesian=True, properties=site_B.properties)

