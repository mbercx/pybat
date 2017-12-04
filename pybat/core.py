# Encoding: utf-8

import itertools

import numpy as np

from monty.json import MSONable

from pymatgen.core import Structure, Element, Molecule, Site
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

    def find_noneq_dimers(self, dimers):
        """
        A script that distills the non-equivalent oxygen dimers from a list of
        dimers.

        Returns:

        """

        noneq_dimers = []
        pass

    def compare_dimers(self, dimers):
        """
        A script that checks if two oxygen dimers in the structure are
        equivalent.

        Args:
            dimers:

        Returns:

        """

        # Obtain the dimer environments
        dimer_environment_A = self.get_dimer_environment(dimers[0])
        dimer_environment_B = self.get_dimer_environment(dimers[1])

        # Transform the dimer environments into molecules
        dimer_A = self.get_dimer_molecule(dimer_environment_A)
        dimer_B = self.get_dimer_molecule(dimer_environment_B)

        # Center the molecules
        pointgroup_analyzer_A = PointGroupAnalyzer(dimer_A)

        print(pointgroup_analyzer_A.symmops)

        # TODO FINISH

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
        self._indices = dimer_indices
        self._sites = list
        self._center = None

    def __eq__(self, other):
        """
        Checks if the dimer environments of two dimers are the same.


        Args:
            other:

        Returns:

        """
        pass

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

            # Determine the indices of the oxygen environment
            dimer_environment_indices = set(self._indices) \
                .union(oxygen_A_neighbors).union(oxygen_B_neighbors)

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

            self._center = sum(
                [site.frac_coords for site in oxygen_sites]) / 2

        return self._center

    def get_dimer_molecule(self):

        molecule_sites = []

        for site in self.sites:

            jimage = site.distance_and_image_from_frac_coords(self.center)[1]
            image_cart_coords = self.cathode.lattice.get_cartesian_coords(
                site.frac_coords - jimage
            )
            molecule_sites.append(Site(site.specie, image_cart_coords))

        return Molecule.from_sites(molecule_sites)

    def visualize_dimer_environment(self, filename=None):
        """
        Remove all the atoms which are not part of the environment of a dimer
        and creates a .xyz file of the oxygen dimer environment in order to
        visualize it more clearly.

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

    cat = Cathode.from_file(structure_file)

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
        cat.visualize_dimer_environment(dimer)
        cat.compare_dimers([dimer, dimers[0]])



