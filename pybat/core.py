# Encoding: utf-8

import itertools

import numpy as np
import math

from monty.json import MSONable

from pymatgen.core import Structure, Element, Molecule, Site, PeriodicSite
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.analysis.chemenv.coordination_environments.voronoi \
    import DetailedVoronoiContainer

"""
Module that contains tools to represent and calculate the properties of
battery cathodes

"""

# TODO Currently, the dimers are defined by their indices. This is a consequence of the fact that the DetailedVoronoiContainer expects indices for its neighbour method. Frankly, I would prefer sites as the basis of the dimer definition as well as it's environment.

# Values for determining the neighbors of a site in a voronoi decomposition
VORONOI_DIST_FACTOR = 1.3
VORONOI_ANG_FACTOR = 0.7

# Tuple of possible cations. This idea should work rather well, considering
# the fact that these cation elements rarely serve another purpose than being
# a cation.
CATIONS = (Element("Li"), Element("Na"), Element("Mg"))

# Tolerance for the template determination. This can be pretty big, since
# the dimer structure is known.
TEMPLATE_DIST_TOL = 5e-1
TEMPLATE_ANGLE_TOL = 2e-1

# Dimer template symmetry permutations
SYMMETRY_PERMUTATIONS = [[1, 2, 4, 3, 6, 5, 8, 7, 9, 10, 11, 12],
                         [2, 1, 3, 4, 7, 8, 5, 6, 11, 12, 9, 10],
                         [1, 2, 3, 4, 5, 6, 7, 8, 10, 9, 12, 11],
                         [2, 1, 4, 3, 8, 7, 6, 5, 11, 12, 9, 10],
                         [1, 2, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11],
                         [2, 1, 4, 3, 8, 7, 6, 5, 12, 11, 10, 9],
                         [2, 1, 3, 4, 7, 8, 5, 6, 12, 11, 10, 9],
                         [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]]

class Cathode(Structure):
    """
    A class representing a cathode material in a battery.

    The Cathode starts in a completely discharged state, i.e. with all cations
    present.
    """
    def __init__(self, lattice, species, coords, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None):

        super(Cathode, self).__init__(
            lattice=lattice, species=species, coords=coords,
            validate_proximity=validate_proximity, to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties
        )

        self._cation_configuration = None
        self._voronoi = None

    @property
    def cation_sites(self):
        """
        A list of all sites which correspond to cations in the Cathode.

        Returns:

        """
        return [site for site in self.sites if site.specie in CATIONS]

    @property
    def cation_configuration(self):
        """
        The configuration of the cations present in the Cathode, i.e. the
        sites of all the present cations.

        Returns:

        """
        if self._cation_configuration is None:
            self._cation_configuration = self.cation_sites

        return self._cation_configuration

    @property
    def voronoi(self):
        """
        ChemEnv voronoi decomposition of the cathode structure.

        Returns:
            pymatgen.analysis.chemenv.coordination_environments.voronoi.\
            DetailedVoronoiContainer

        """
        if self._voronoi is None:
            self._voronoi = DetailedVoronoiContainer(self)
            return self._voronoi
        else:
            return self._voronoi

    def remove_cations(self, sites=None):
        """
        Remove the cations from the cathode, i.e. delithiate the structure in
        case Li is the cation of the cathode.

        :return:
        """

        # If no indices are given
        if sites is None:
            # Remove all the cations
            self._cation_configuration = []

        # Else remove the requested sites from the cation configuration
        else:
            for site in sites:
                # Check if the site provided corresponds to a cation site
                if not site.specie in CATIONS:
                    raise IOError("Provided indices do not all correspond to "
                                  "a cation site!")
                else:
                    # Remove the cation site
                    try:
                        self._cation_configuration.remove(site)
                    except ValueError:
                        raise Warning("Requested site not found in cation "
                                      "configuration.")

    def find_cation_configurations(self):
        """
        Plan is to find all non-equivalent cation configurations. Is probably
        already implemented elsewhere.

        Returns:

        """
        raise NotImplementedError

    def change_site_distance(self, site_indices, distance):
        """
        Change the coordinates of two sites in a structure in order to adjust
        their distance.

        Args:
            self: (pymatgen.core.structure.Structure)
            site_indices:

        Returns:
            None
        """

        site_A = self.sites[site_indices[0]]
        site_B = self.sites[site_indices[1]]

        # Find the distance between the sites, as well as the image of site B
        # closest to site A
        (original_distance, closest_image_B) = site_A.distance_and_image(
            site_B)

        image_cart_coords = self.lattice.get_cartesian_coords(
            site_B.frac_coords + closest_image_B
        )

        # Calculate the vector that connects site A with site B
        connection_vector = image_cart_coords - site_A.coords

        # Make it a unit vector
        connection_vector /= np.linalg.norm(connection_vector)

        # Calculate the distance the sites need to be moved.
        site_move_distance = (original_distance - distance) / 2

        # Calculate the new cartesian coordinates of the sites
        new_site_A_coords = site_A.coords + site_move_distance * connection_vector
        new_site_B_coords = site_B.coords - site_move_distance * connection_vector

        # Change the sites in the structure
        self.replace(i=site_indices[0], species=site_A.species_string,
                     coords=new_site_A_coords,
                     coords_are_cartesian=True,
                     properties=site_A.properties)

        self.replace(i=site_indices[1], species=site_B.species_string,
                     coords=new_site_B_coords,
                     coords_are_cartesian=True,
                     properties=site_B.properties)

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


class LiRichCathode(Cathode):
    """
    A class representing a Li-rich cathode material.

    """
    def __init__(self, lattice, species, coords, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None):

        super(LiRichCathode, self).__init__(
            lattice=lattice, species=species, coords=coords,
            validate_proximity=validate_proximity, to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties
        )

    def find_oxygen_dimers(self, site_index):
        """
        Returns a list of index pairs corresponding to the oxygen dimers that
        can be formed around the site provided by the user, i.e. with oxygens
        that are neighbours of the provided site.

        Args:
            site_index:

        Returns:

        """

        # Determine the oxygen neighbors for the provided site
        oxygen_neighbors_indices = [
            neighbor["index"] for neighbor
            in self.voronoi.neighbors(site_index, VORONOI_DIST_FACTOR,
                                      VORONOI_ANG_FACTOR)
            if self.sites[neighbor["index"]].specie == Element["O"]
        ]

        if len(oxygen_neighbors_indices) <= 1:
            raise ValueError("Provided site does not have two oxygen "
                             "neighbours.\n")

        # Find all oxygen neighbour combinations that share another neighbour
        oxygen_dimers = []
        oxygen_combinations = itertools.combinations(
            oxygen_neighbors_indices, 2
        )

        for oxygen_pair in oxygen_combinations:

            oxygen_A_neighbors = [
                neighbor["index"] for neighbor
                in self.voronoi.neighbors(oxygen_pair[0], VORONOI_DIST_FACTOR,
                                          VORONOI_ANG_FACTOR)
            ]
            oxygen_B_neighbors = [
                neighbor["index"] for neighbor
                in self.voronoi.neighbors(oxygen_pair[1], VORONOI_DIST_FACTOR,
                                          VORONOI_ANG_FACTOR)
            ]

            # Check how many neighbours the oxygen have in common
            shared_neighbours = list(
                set(oxygen_A_neighbors).intersection(oxygen_B_neighbors)
            )

            if len(shared_neighbours) == 2:
                oxygen_dimers.append(oxygen_pair)

        return oxygen_dimers

    def remove_dimer_cations(self, dimer_indices):
        """

        Args:
            dimer_indices:

        Returns:

        """

        # Find the indices of the dimer environment
        dimer_environment = self.get_dimer_environment(dimer_indices)

        remove_sites = [site for site in dimer_environment
                       if site.specie in CATIONS]

        self.remove_cations(remove_sites)

    def find_noneq_dimers(self, site_index=None):
        """
        A script that distills the non-equivalent oxygen dimers around a site
        index.

        Returns:

        """

        noneq_dimers = []

        if site_index is None:

            ignore = set(CATIONS).union((Element("O"),))

            transition_metal_indices = [index for index in
                                        range(len(self.sites))
                                        if self.sites[index].specie
                                        not in ignore]

            for index in transition_metal_indices:

                dimers = [Dimer(self, dimer_indices) for dimer_indices
                          in self.find_oxygen_dimers(index)]

                for dimer in dimers:
                    if dimer not in noneq_dimers:
                        noneq_dimers.append(dimer)

        else:

            dimers = [Dimer(self, dimer_indices) for dimer_indices
                      in self.find_oxygen_dimers(site_index)]

            for dimer in dimers:
                if dimer not in noneq_dimers:
                    noneq_dimers.append(dimer)

        return noneq_dimers


class Dimer(MSONable):
    """
    Class representing an oxygen dimer in a Li-rich cathode structure.

    """

    # TODO Give the definition of this class a good hard thinking over.

    def __init__(self, cathode, dimer_indices):
        """
        Initialize

        Args:
            cathode
            dimer_indices

        Returns:
            pybat.core.Dimer

        """
        self._cathode = cathode
        self._indices = tuple(dimer_indices)
        self._sites = list
        self._center = None
        self._template = list

    def __eq__(self, other):
        """
        Checks if the dimer environments of two dimers are the same.


        Args:
            other:

        Returns:

        """
        is_equal = False

        for permutation in SYMMETRY_PERMUTATIONS:

            if [self.template[index] for index in range(1,13)] == \
                    [other.template[key] for key in permutation]:

                is_equal = True

        return is_equal

    @property
    def cathode(self):
        return self._cathode

    @property
    def indices(self):
        return self._indices

    @property
    def sites(self):

        if self._sites is list:
            # Find the oxygen neighbours
            oxygen_A_neighbors = [
                neighbor["index"] for neighbor
                in self.cathode.voronoi.neighbors(self.indices[0],
                                                  VORONOI_DIST_FACTOR,
                                                  VORONOI_ANG_FACTOR)
            ]
            oxygen_B_neighbors = [
                neighbor["index"] for neighbor
                in self.cathode.voronoi.neighbors(self.indices[1],
                                                  VORONOI_DIST_FACTOR,
                                                  VORONOI_ANG_FACTOR)
            ]

            # Determine the indices of the oxygen environment. The indices are
            # sorted in such a way that the oxygen indices come first,
            # followed by the indices of the shared neighbors.
            # TODO This can be done better. Really.
            shared_neighbors = tuple(set(oxygen_A_neighbors).intersection(oxygen_B_neighbors))
            other_neighbors = set(oxygen_A_neighbors).union(oxygen_B_neighbors)
            other_neighbors.remove(shared_neighbors[0])
            other_neighbors.remove(shared_neighbors[1])
            other_neighbors = tuple(other_neighbors)

            dimer_environment_indices = self._indices + shared_neighbors \
                                        + other_neighbors

            # Recover the corresponding sites
            self._sites = [self.cathode.sites[index] for index
                                 in dimer_environment_indices]

        return self._sites

    @property
    def center(self):

        if self._center is None:

            # Find the center of the oxygen sites
            oxygen_sites = []
            for site in self.sites:
                if site.specie == Element("O"):
                    oxygen_sites.append(site)

            (distance, oxygen_image) = oxygen_sites[0].distance_and_image(oxygen_sites[1])

            image_cart_coords = oxygen_sites[1].coords \
                                + np.dot(oxygen_image,
                                         self.cathode.lattice.matrix)

            self._center = ( oxygen_sites[0].coords + image_cart_coords ) / 2

            # print("Oxygen 1 coordinates:")
            # print(oxygen_sites[0].coords)
            #
            # print("Oxygen 2 coordinates:")
            # print(oxygen_sites[1].coords)
            #
            # print("Dimer center:")
            # print(self._center)
            #
            # print("Found distance using pmg")
            # print(distance)
            #
            # print("Found distance using cart coords")
            # print(np.linalg.norm(oxygen_sites[0].coords - image_cart_coords))
            #
            # print()

        return self._center

    @property
    def template(self):
        if self._template is list:

            # A template # TODO Add definition

            # The template idea seems like the fastest way of being able to
            # compare two dimers. By reducing them to a template, we are able
            # compare dimers by applying permutations that represent symmetry
            # transformations. This is not a very general approach, clearly.

            # Note that the sites property is built in such a way that the
            # first two sites correspond to the oxygen atoms and the next two
            # correspond to the shared neighbours. This convention is made to
            # save us some work here.

            dimer_molecule = self.get_dimer_molecule()

            oxy_1 = dimer_molecule.sites[0]
            oxy_2 = dimer_molecule.sites[1]

            shared_neighbor_3 = dimer_molecule.sites[2]
            shared_neighbor_4 = dimer_molecule.sites[3]

            # The template is defined as a dictionary between site numbers and
            # dimer environment sites
            template = {1:oxy_1.specie,
                        2:oxy_2.specie,
                        3:shared_neighbor_3.specie,
                        4:shared_neighbor_4.specie}

            # Loop over the remaining sites to find their template positions
            for site in dimer_molecule.sites[4:]:

                # Find the sites which are in the plane of the oxygens and
                # their shared neighbors.
                if np.linalg.norm(oxy_1.coords
                        - (shared_neighbor_4.coords - oxy_1.coords)
                        - site.coords) < TEMPLATE_DIST_TOL:

                    template[5] = site.specie

                if np.linalg.norm(oxy_1.coords
                        - (shared_neighbor_3.coords - oxy_1.coords)
                        - site.coords) < TEMPLATE_DIST_TOL:

                    template[6] = site.specie

                if np.linalg.norm(oxy_2.coords
                        - (shared_neighbor_4.coords - oxy_2.coords)
                        - site.coords) < TEMPLATE_DIST_TOL:

                    template[7] = site.specie

                if np.linalg.norm(oxy_2.coords
                        - (shared_neighbor_3.coords - oxy_2.coords)\
                        - site.coords) < TEMPLATE_DIST_TOL:

                    template[8] = site.specie

                # Find the sites which are out of plane
                oxy_1_oop = np.cross(
                    shared_neighbor_4.coords - oxy_1.coords,
                    shared_neighbor_3.coords - oxy_1.coords
                )
                oxy_2_oop = np.cross(
                    shared_neighbor_3.coords - oxy_2.coords,
                    shared_neighbor_4.coords - oxy_2.coords
                )

                if angle_between(oxy_1_oop, site.coords - oxy_1.coords) < TEMPLATE_ANGLE_TOL:

                    template[9] = site.specie

                if angle_between(oxy_1_oop, site.coords - oxy_1.coords) > math.pi - TEMPLATE_ANGLE_TOL:

                    template[10] = site.specie

                if angle_between(oxy_2_oop, site.coords - oxy_2.coords) < TEMPLATE_ANGLE_TOL:

                    template[11] = site.specie

                if angle_between(oxy_2_oop, site.coords - oxy_2.coords) > math.pi - TEMPLATE_ANGLE_TOL:

                    template[12] = site.specie

            self._template = template

            if len(template) != len(self.sites):
                raise ValueError("Template creation failed!")

        return self._template

    def get_dimer_molecule(self):

        molecule_sites = []

        for site in self.sites:

            # print()
            # print("Distance between image and oxygen center:")

            (distance, jimage) = site.distance_and_image_from_frac_coords(
                self.cathode.lattice.get_fractional_coords(self.center))

            image_cart_coords = site.coords \
                                - np.dot(jimage, self.cathode.lattice.matrix)

            # print(distance)
            #
            # print("Distance found using cart coords")
            # print(np.linalg.norm(image_cart_coords - self.center))

            molecule_sites.append(Site(site.specie, image_cart_coords))

        return Molecule.from_sites(molecule_sites)

    def visualize_dimer_environment(self, filename=None):
        """
        Creates a .xyz file of the oxygen dimer environment in order to
        visualize it in e.g. VESTA.

        Returns:

        """

        # Turn the structure into a molecule
        dimer_environment_molecule = self.get_dimer_molecule()

        if filename == None:
            filename = str(self.cathode.composition).replace(" ", "") + "_" \
                       + str(self.indices[0]) + "_" + str(self.indices[1])

        dimer_environment_molecule.to("xyz", filename + ".xyz")


    # TODO Check if MSONable methods need to be implemented

    def as_dict(self):
        pass

    @classmethod
    def from_dict(cls, d):
        pass


def test_script(structure_file):

    cat = LiRichCathode.from_file(structure_file)

    print(cat)
    site_index = int(input("Please give the site around which you would like to "
                       "study O-O dimers: "))

    dimers = cat.find_oxygen_dimers(site_index)

    print("Oxygen dimers of site")
    print(dimers)

    for dimer in dimers:
        print("")
        print("Dimer")
        print(dimer)
        print("")
        Dimer(cat, dimer).visualize_dimer_environment()


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """
    Returns the angle in radians between vectors 'v1' and 'v2'::
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

# def rotation_matrix(axis, theta):
#     """
#     Return the rotation matrix associated with clockwise rotation about
#     the given axis by theta radians.
#     """
#     axis = np.asarray(axis)
#     axis = axis/math.sqrt(np.dot(axis, axis))
#     a = math.cos(theta/2.0)
#     b, c, d = axis*math.sin(theta/2.0)
#     aa, bb, cc, dd = a*a, b*b, c*c, d*d
#     bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
#     return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
#                      [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
#                      [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def permute(iterable, permutation):

    if len(iterable) != len(permutation):
        raise ValueError("Length of list does not match permutation length!")

    if len(set([i for i in range(1, len(iterable)+1)]).intersection(permutation)) \
        != len(permutation):

        raise ValueError("Permutation is ill-defined.")

    return [iterable[index - 1] for index in permutation]
