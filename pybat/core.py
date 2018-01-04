# Encoding: utf-8

import itertools
import math
import json

import numpy as np

from monty.io import zopen
from monty.json import MSONable
from fnmatch import fnmatch

from pymatgen.core import Structure, Element, Molecule, Site, PeriodicSite
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

# Tolerance for the representation determination. This can be pretty big, since
# the dimer structure is known.
REPRESENTATION_DIST_TOL = 5e-1
REPRESENTATION_ANGLE_TOL = 2e-1

# Dimer representation symmetry permutations
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

    @cation_configuration.setter
    def cation_configuration(self, configuration):
        self._cation_configuration = configuration

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

    @voronoi.setter
    def voronoi(self, voronoi_container):
        self._voronoi = voronoi_container

    def remove_cations(self, sites=None):
        """
        Remove the cations from the cathode, i.e. delithiate the structure in
        case Li is the cation of the cathode.

        Note that this does not remove the sites from the pymatgen Structure.
        The cation_configuration of the Cathode is simply adjusted by removing
        the requested cations from this List.

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
                if site.specie not in CATIONS:
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
            site_indices:
            distance:

        Returns:

        """

        site_a = self.sites[site_indices[0]]
        site_b = self.sites[site_indices[1]]

        # Find the distance between the sites, as well as the image of site B
        # closest to site A
        (original_distance, closest_image_b) = site_a.distance_and_image(
            site_b)

        image_cart_coords = self.lattice.get_cartesian_coords(
            site_b.frac_coords + closest_image_b
        )

        # Calculate the vector that connects site A with site B
        connection_vector = image_cart_coords - site_a.coords

        # Make it a unit vector
        connection_vector /= np.linalg.norm(connection_vector)

        # Calculate the distance the sites need to be moved.
        site_move_distance = (original_distance - distance) / 2

        # Calculate the new cartesian coordinates of the sites
        new_site_a_coords = site_a.coords + site_move_distance\
            * connection_vector
        new_site_b_coords = site_b.coords - site_move_distance\
            * connection_vector

        # Change the sites in the structure
        self.replace(i=site_indices[0], species=site_a.species_string,
                     coords=new_site_a_coords,
                     coords_are_cartesian=True,
                     properties=site_a.properties)

        self.replace(i=site_indices[1], species=site_b.species_string,
                     coords=new_site_b_coords,
                     coords_are_cartesian=True,
                     properties=site_b.properties)

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

    def to(self, fmt=None, filename=None, **kwargs):
        """
        Outputs the structure to a file or string.

        Overwritten from IStructure in order to consider the cation
        configuration for .cif and VASP POSCAR files. For .json formats the
        as_dict method successfully stores the cation configuration. Importing
        the Cathode object from .cif or POSCAR files

        Args:
            fmt (str): Format to output to. Defaults to JSON unless filename
                is provided. If fmt is specifies, it overrides whatever the
                filename is. Options include "cif", "poscar", "cssr", "json".
                Non-case sensitive.
            filename (str): If provided, output will be written to a file. If
                fmt is not specified, the format is determined from the
                filename. Defaults is None, i.e. string output.

        Returns:
            (str) if filename is None. None otherwise.

        """

        if fmt in ["cif", "poscar"] or fnmatch(filename, "*.cif*") \
            or fnmatch(filename, "POSCAR"):

            structure = self.copy()
            for cation_site in self.cation_sites:
                if cation_site not in self.cation_configuration:
                    structure.remove(cation_site)

            super(Cathode, structure).to(fmt=fmt, filename=filename, **kwargs)

        elif fmt == "json" or fnmatch(filename, "*.json*"):

            super(Cathode, self).to(fmt=fmt, filename=filename, **kwargs)

        else:
            raise NotImplementedError("Only json, cif or VASP POSCAR formats "
                                      "are currently supported.")

    @classmethod
    def from_structure(cls, structure):

        return cls.from_sites(structure.sites)

    @classmethod
    def from_str(cls, input_string, fmt, primitive=False, sort=False,
                 merge_tol=0.0):
        if fmt is not "json":
            return super(Cathode, cls).from_str(input_string, fmt, primitive,
                                                sort, merge_tol)
        else:
            d = json.loads(input_string)
            return cls.from_dict(d)

    def as_dict(self, verbosity=1, fmt=None, **kwargs):

        d = super(Cathode, self).as_dict(verbosity=verbosity,
                                         fmt=fmt, **kwargs)

        d["cation_configuration"] = [cation.as_dict() for cation
                                     in self.cation_configuration]

        # Note The voronoi decomposition can unfortunately not be MSONabled,
        # because of a recursion issue where the as_dict method of the
        # DetailedVoronoiContainer uses the as_dict() method from the structure
        # it is based on (infinite recursion).

        # TODO Add voronoi decomposition, by adding the dictionary components manually

        return d

    @classmethod
    def from_dict(cls, d, fmt=None):

        structure = super(Cathode, cls).from_dict(d)
        cathode = cls.from_structure(structure)
        cation_configuration = d.get("cation_configuration", None)

        if cation_configuration is not None:
            cathode.cation_configuration = [PeriodicSite.from_dict(cation) for
                                            cation in cation_configuration]

        return cathode


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

            oxygen_a_neighbors = [
                neighbor["index"] for neighbor
                in self.voronoi.neighbors(oxygen_pair[0], VORONOI_DIST_FACTOR,
                                          VORONOI_ANG_FACTOR)
            ]
            oxygen_b_neighbors = [
                neighbor["index"] for neighbor
                in self.voronoi.neighbors(oxygen_pair[1], VORONOI_DIST_FACTOR,
                                          VORONOI_ANG_FACTOR)
            ]

            # Check how many neighbours the oxygen have in common
            shared_neighbours = list(
                set(oxygen_a_neighbors).intersection(oxygen_b_neighbors)
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

        # Find the sites of the dimer environment
        dimer_environment = Dimer(self, dimer_indices).sites

        remove_sites = [site for site in dimer_environment
                        if site.specie in CATIONS]

        self.remove_cations(remove_sites)

    def find_noneq_dimers(self, site_index=None):
        """
        A script that distills the non-equivalent oxygen dimers around a site
        index.

        In case no site index is provided, the method will loop over all sites
        which do not contain oxygen or a cation.

        Returns:
            List of pybat.Dimer objects

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
    Class definition of an oxygen dimer in a Li-rich cathode structure.

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
        self._representation = list

    def __eq__(self, other):
        """
        Checks if the dimer environments of two dimers are the same.


        Args:
            other:

        Returns:

        """
        is_equal = False

        for permutation in SYMMETRY_PERMUTATIONS:

            if [self.representation[index] for index in range(1, 13)] == \
                    [other.representation[key] for key in permutation]:

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
            oxygen_a_neighbors = [
                neighbor["index"] for neighbor
                in self.cathode.voronoi.neighbors(self.indices[0],
                                                  VORONOI_DIST_FACTOR,
                                                  VORONOI_ANG_FACTOR)
            ]
            oxygen_b_neighbors = [
                neighbor["index"] for neighbor
                in self.cathode.voronoi.neighbors(self.indices[1],
                                                  VORONOI_DIST_FACTOR,
                                                  VORONOI_ANG_FACTOR)
            ]

            # Determine the indices of the oxygen environment. The indices are
            # sorted in such a way that the oxygen indices come first,
            # followed by the indices of the shared neighbors.
            # TODO This can be done better. Really.
            shared_neighbors = tuple(
                set(oxygen_a_neighbors).intersection(oxygen_b_neighbors)
            )
            other_neighbors = set(oxygen_a_neighbors).union(oxygen_b_neighbors)
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

            (distance, oxygen_image) = oxygen_sites[0].distance_and_image(
                oxygen_sites[1])

            image_cart_coords = oxygen_sites[1].coords \
                + np.dot(oxygen_image, self.cathode.lattice.matrix)

            self._center = (oxygen_sites[0].coords + image_cart_coords) / 2

        return self._center

    @property
    def representation(self):
        if self._representation is list:

            # A representation # TODO Add definition

            # The representation idea seems like the fastest way of being able
            # to compare two dimers. By reducing them to a representation, we
            # are able compare dimers by applying permutations that represent
            # symmetry transformations. This is not a very general approach,
            # clearly.

            # Note that the sites property is built in such a way that the
            # first two sites correspond to the oxygen atoms and the next two
            # correspond to the shared neighbours. This convention is made to
            # save us some work here.

            dimer_molecule = self.get_dimer_molecule()

            oxy_1 = dimer_molecule.sites[0]
            oxy_2 = dimer_molecule.sites[1]

            shared_neighbor_3 = dimer_molecule.sites[2]
            shared_neighbor_4 = dimer_molecule.sites[3]

            # The representation is defined as a dictionary between site
            # numbers and dimer environment sites
            representation = {1:oxy_1.specie,
                              2:oxy_2.specie,
                              3:shared_neighbor_3.specie,
                              4:shared_neighbor_4.specie}

            # TODO Find a cleaner way of assigning representation positions

            # Loop over the remaining sites to find their representation
            # positions
            for site in dimer_molecule.sites[4:]:

                # Find the sites which are in the plane of the oxygens and
                # their shared neighbors.
                if np.linalg.norm(oxy_1.coords
                        - (shared_neighbor_4.coords - oxy_1.coords)
                        - site.coords) < REPRESENTATION_DIST_TOL:

                    representation[5] = site.specie

                if np.linalg.norm(oxy_1.coords
                        - (shared_neighbor_3.coords - oxy_1.coords)
                        - site.coords) < REPRESENTATION_DIST_TOL:

                    representation[6] = site.specie

                if np.linalg.norm(oxy_2.coords
                        - (shared_neighbor_4.coords - oxy_2.coords)
                        - site.coords) < REPRESENTATION_DIST_TOL:

                    representation[7] = site.specie

                if np.linalg.norm(oxy_2.coords
                        - (shared_neighbor_3.coords - oxy_2.coords)\
                        - site.coords) < REPRESENTATION_DIST_TOL:

                    representation[8] = site.specie

                # Find the sites which are out of plane
                oxy_1_oop = np.cross(
                    shared_neighbor_4.coords - oxy_1.coords,
                    shared_neighbor_3.coords - oxy_1.coords
                )
                oxy_2_oop = np.cross(
                    shared_neighbor_3.coords - oxy_2.coords,
                    shared_neighbor_4.coords - oxy_2.coords
                )

                if angle_between(oxy_1_oop, site.coords - oxy_1.coords) \
                        < REPRESENTATION_ANGLE_TOL:

                    representation[9] = site.specie

                if angle_between(oxy_1_oop, site.coords - oxy_1.coords) \
                        > math.pi - REPRESENTATION_ANGLE_TOL:

                    representation[10] = site.specie

                if angle_between(oxy_2_oop, site.coords - oxy_2.coords) \
                        < REPRESENTATION_ANGLE_TOL:

                    representation[11] = site.specie

                if angle_between(oxy_2_oop, site.coords - oxy_2.coords) \
                        > math.pi - REPRESENTATION_ANGLE_TOL:

                    representation[12] = site.specie

            self._representation = representation

            if len(representation) != len(self.sites):
                raise ValueError("Failed to create dimer representation.")

        return self._representation

    def get_dimer_molecule(self):

        molecule_sites = []

        for site in self.sites:

            (distance, jimage) = site.distance_and_image_from_frac_coords(
                self.cathode.lattice.get_fractional_coords(self.center))

            image_cart_coords = \
                site.coords - np.dot(jimage, self.cathode.lattice.matrix)

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

        if filename is None:
            filename = str(self.cathode.composition).replace(" ", "") + "_" \
                       + str(self.indices[0]) + "_" + str(self.indices[1])

        dimer_environment_molecule.to("xyz", filename + ".xyz")

    # TODO Check if MSONable methods need to be implemented

    def to(self, fmt="json", filename=None):

        if fmt == "json":
            if filename:
                with zopen(filename, "wt", encoding='utf8') as file:
                    return json.dump(self.as_dict(), file)
            else:
                return json.dumps(self.as_dict())
        else:
            raise NotImplementedError("Currently only json format is "
                                      "supported.")
    @classmethod
    def from_str(cls, input_string, fmt="json"):
        """
        Initialize a Facet from a string.

        Currently only supports 'json' formats.

        Args:
            input_string (str): String from which the Facet is initialized.
            fmt (str): Format of the string representation.

        Returns:
            (*cage.Facet*)
        """
        if fmt == "json":
            d = json.loads(input_string)
            return cls.from_dict(d)
        else:
            raise NotImplementedError('Only json format has been '
                                      'implemented.')

    @classmethod
    def from_file(cls, filename):

        with zopen(filename) as file:
            contents = file.read()

        return cls.from_str(contents)


    def as_dict(self):

        d = {}

        d["cathode"] = self.cathode.as_dict()
        d["dimer_indices"] = self.indices

        return d

    @classmethod
    def from_dict(cls, d):

        return cls(cathode=LiRichCathode.from_dict(d["cathode"]),
                   dimer_indices=d["dimer_indices"])


def test_script(structure_file):

    cat = LiRichCathode.from_file(structure_file)

    print(cat)
    site_index = int(input("Please give the site around which you would like "
                           "to study O-O dimers: "))

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

    permutation_numbers = set([i for i in range(1, len(iterable)+1)])

    if len(permutation_numbers.intersection(permutation)) != len(permutation):
        raise ValueError("Permutation is ill-defined.")

    return [iterable[index - 1] for index in permutation]
