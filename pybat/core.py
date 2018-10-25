# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import itertools
import math
import json
import os

import numpy as np

from monty.io import zopen
from monty.json import MSONable
from pymatgen.core import Structure, Composition, Molecule, Site, Element
from pymatgen.analysis.chemenv.coordination_environments.voronoi \
    import DetailedVoronoiContainer
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

"""
Module that contains tools to represent and calculate the properties of
battery cathodes.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "May 2018"

# TODO Currently, the dimers are defined by their indices.
# This is a consequence of the fact that the DetailedVoronoiContainer
# expects indices for its neighbor method. Frankly, I would prefer sites as
# the basis of the dimer definition as well as it's environment.

# Values for determining the neighbors of a site in a voronoi decomposition
VORONOI_DIST_FACTOR = 1.3
VORONOI_ANG_FACTOR = 0.7

# Tuple of possible cations. This idea should work fine, considering the
# fact that these cation elements rarely serve another purpose than being
# a cation.
CATIONS = ("Li", "Na", "Mg")

# Tolerance for determining whether two oxygens are on opposite sides of an
# octahedron. If the angle between the two vectors connecting the site and
# the corresponding oxygens is larger than this value, the oxygens are
# considered to be opposites.
OXYGEN_ANGLE_TOL = math.pi * 9 / 10

# Tolerance for the representation determination. This can be pretty big, since
# the dimer environment structure is known.
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

    The main idea of this class is to keep track of the original sites by
    using sites with empty Compositions. This is important to make sure that
    the voronoi decomposition is successful, and hence interesting if we
    want to look at coordinations and neighbors. Another advantage is that
    we can consider the empty cation sites for final positions of transition
    metal migrations.

    """

    def __init__(self, lattice, species, coords, charge=None,
                 validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None):

        super(Cathode, self).__init__(
            lattice=lattice, species=species, coords=coords,
            validate_proximity=validate_proximity, to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties
        )

        self._voronoi = None

    @property
    def cation_configuration(self):
        """
        A list of all sites which correspond to cations in the Cathode.

        Returns:

        """
        return [site for site in self.sites if site.species_string in CATIONS]

    @cation_configuration.setter
    def cation_configuration(self, configuration):

        # TODO Add checks

        # Remove all cations
        for cation in [cat for cat in CATIONS if Element[cat] in set(
                self.composition.keys())]:
            self.replace_species({cation: {cation: 0}})

        # Add the cation sites
        if isinstance(configuration, dict):
            for cation in configuration.keys():
                for index in configuration[cation]:
                    self.replace(index, cation, properties={"magmom": 0})

        elif all([isinstance(item, Site) for item in configuration]):
            for catsite in configuration:
                for i, site in enumerate(self):
                    if np.linalg.norm(site.distance(catsite)) < 0.05:
                        self.replace(i, catsite.specie,
                                     properties={"magmom": 0})

        else:
            raise TypeError("Cation configurations should be a dictionary "
                            "mapping cations to site indices or a list of "
                            "sites.")

    def __str__(self):
        """
        Overwritten string representation, in order to provide information
        about the cation configuration, as well as the VESTA index, which
        is useful when defining structural changes.

        """
        outs = ["Full Formula ({s})".format(s=self.composition.formula),
                "Reduced Formula: {}".format(self.composition.reduced_formula)]
        to_s = lambda x: "%0.6f" % x
        outs.append("abc   : " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.abc]))
        outs.append("angles: " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.angles]))
        if self._charge:
            if self._charge >= 0:
                outs.append("Overall Charge: +{}".format(self._charge))
            else:
                outs.append("Overall Charge: -{}".format(self._charge))
        outs.append("Sites ({i})".format(i=len(self)))
        data = []
        props = self.site_properties
        keys = sorted(props.keys())
        vesta_index = 1
        for i, site in enumerate(self):
            if site.species_and_occu.num_atoms == 0:
                row = [str(i), "-", "Vac"]

            else:
                row = [str(i), vesta_index, site.species_string]
                vesta_index += 1

            row.extend([to_s(j) for j in site.frac_coords])
            for k in keys:
                row.append(props[k][i])
            data.append(row)

        from tabulate import tabulate
        outs.append(
            tabulate(data,
                     headers=["#", "#VESTA", "SP", "a", "b", "c"] + keys,
                     ))
        return "\n".join(outs)

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

    def add_cations(self, sites=None):
        """
        Args:
            sites:

        Returns:

        """

        # Add the cation sites
        if isinstance(sites, dict):
            for cation in sites.keys():
                for index in sites[cation]:
                    self.replace(index, cation, properties={"magmom": 0})

        elif all([isinstance(item, Site) for item in sites]):
            for catsite in sites:
                for i, site in enumerate(self):
                    if np.linalg.norm(site.distance(catsite)) < 0.05:
                        self.replace(i, catsite.specie,
                                     properties={"magmom": 0})

        else:
            raise TypeError("Cation configurations should be a dictionary "
                            "mapping cations to site indices or a list of "
                            "sites.")

    def remove_cations(self, sites=None):
        """
        Remove the cations from the cathode, i.e. delithiate the structure in
        case Li is the cation of the cathode.

        Note that this does not remove the sites from the pymatgen Structure.
        The occupancy is simply adjusted to an empty Composition object.

        Args:
            sites: List of indices
            List of pymatgen.core.Sites which are to be removed.

        """

        # TODO add checks

        # If no sites are given
        if sites is None:
            # Remove all the cations
            self.cation_configuration = []

        # If a List of integers is given
        elif all([isinstance(item, int) for item in sites]):
            for index in sites:
                self.replace(index, Composition(), properties={"magmom": 0})

        # If a List of sites is given
        elif all([isinstance(item, Site) for item in sites]):
            for site in sites:
                # Check if the provided site corresponds to a cation site
                if site in self.cation_configuration:
                    cat_conf = self.cation_configuration.copy()
                    cat_conf.remove(site)
                    self.cation_configuration = cat_conf
                else:
                    raise Warning("Requested site not found in cation "
                                  "configuration.")
        else:
            raise IOError("Incorrect site input.")

    def change_site_distance(self, site_indices, distance):
        """
        Change the coordinates of two sites in a structure in order to adjust
        their distance.

        Args:
            site_indices:
            distance:
        """

        # TODO Add possibility of site_indices simply being the sites

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
        new_site_a_coords = site_a.coords \
                            + site_move_distance * connection_vector
        new_site_b_coords = site_b.coords \
                            - site_move_distance * connection_vector

        # Change the sites in the structure
        self.replace(i=site_indices[0], species=site_a.species_string,
                     coords=new_site_a_coords,
                     coords_are_cartesian=True,
                     properties=site_a.properties)

        self.replace(i=site_indices[1], species=site_b.species_string,
                     coords=new_site_b_coords,
                     coords_are_cartesian=True,
                     properties=site_b.properties)

    def update_sites(self, directory, ignore_magmom=False):
        """
        Based on the CONTCAR and OUTCAR of a geometry optimization, update the
        site coordinates and magnetic moments that were optimized. Note that
        this method relies on the cation configuration of the not having
        changed.

        Args:
            directory (str): Directory in which the geometry optimization
                output files (i.e. CONTCAR and OUTCAR) are stores.
        """

        new_cathode = Cathode.from_file(os.path.join(directory, "CONTCAR"))

        out = Outcar(os.path.join(directory, "OUTCAR"))

        if ignore_magmom:
            new_cathode.add_site_property("magmom",
                                          self.site_properties["magmom"])
        else:
            magmom = [site["tot"] for site in out.magnetization]
            new_cathode.add_site_property("magmom", magmom)

        # Update the lattice
        self.modify_lattice(new_cathode.lattice)

        # Update the coordinates of the occupied sites.
        new_index = 0
        for i, site in enumerate(self):

            # If the site is not empty
            if site.species_and_occu != Composition():
                new_site = new_cathode.sites[new_index]
                # Update the site coordinates
                self.replace(i, species=new_site.species_and_occu,
                             coords=new_site.frac_coords,
                             properties=new_site.properties)
                new_index += 1

    def find_noneq_cations(self):
        """
        Find a list of the site indices of all non-equivalent cations.

        Returns:
            List of site indices

        """
        symmops = SpacegroupAnalyzer(self).get_space_group_operations()

        cation_indices = [
            index for index in range(len(self.sites))
            if not self.sites[index].species_string == "O"
        ]

        # Start with adding the first cation
        inequiv_cations = [cation_indices[0],]

        for index in cation_indices[1:]:

            s1 = [self.sites[index], ]

            # Check if the site is equivalent with one of the sites in the
            # inequivalent list.
            inequivalent = True

            for inequive_index in inequiv_cations:

                s2 = [self.sites[inequive_index], ]

                if symmops.are_symmetrically_equivalent(s1, s2):
                    inequivalent = False

            if inequivalent:
                inequiv_cations.append(index)

        return inequiv_cations

    def find_cation_configurations(self):
        """
        Plan is to find all non-equivalent cation configurations. Is probably
        already implemented elsewhere.

        Returns:

        """
        raise NotImplementedError

    def set_to_high_spin(self):
        """

        :return:
        """
        raise NotImplementedError

    def set_to_low_spin(self):
        """

        :return:
        """
        raise NotImplementedError

    def as_ordered_structure(self):
        """
        Return the structure as a pymatgen.core.Structure, removing the
        unoccupied sites. This is because many of the IO methods of pymatgen
        run into issues when empty occupancies are present.

        Returns:
            pymatgen.core.Structure

        """

        return Structure.from_sites(
            [site for site in self.sites
             if site.species_and_occu != Composition()]
        )

    def to(self, fmt=None, filename=None, **kwargs):
        """
        Structure method override to solve issue with writing the Cathode to a
        POSCAR file

        # TODO Figure out what exactly was the problem here again... Should
        have written this down immediately! I think it had something to do
        with the order of the Sites changing...

        Args:
            fmt:
            filename:
            **kwargs:

        Returns:

        """

        if fmt == "poscar":
            structure = self.as_ordered_structure()
            structure.to(fmt, filename, **kwargs)

        else:
            super(Cathode, self).to(fmt, filename, **kwargs)

    @classmethod
    def from_structure(cls, structure):
        """
        Initializes a Cathode from a pymatgen.core.Structure.

        Args:
            structure (pymatgen.core.Structure): Structure from which to
            initialize the Cathode.

        Returns:
            pybat.core.Structure

        """

        return cls.from_sites(structure.sites)


class LiRichCathode(Cathode):
    """
    A class representing a Li-rich cathode material.

    """

    def __init__(self, lattice, species, coords, charge=None,
                 validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None):

        super(LiRichCathode, self).__init__(
            lattice=lattice, species=species, coords=coords,
            validate_proximity=validate_proximity, to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties
        )

    def find_oxygen_dimers(self, site_index=None,
                           oxygen_angle_tol=OXYGEN_ANGLE_TOL):
        """
        Returns a list of index pairs corresponding to the oxygen dimers that
        can be formed around the site provided by the user, i.e. with oxygens
        that are neighbours of the provided site.

        Args:
            site_index:

        Returns:

        """

        if site_index is None:

            # TODO This is a lazy way of dealing with this...
            # Just using the method for single indices and then filtering
            # duplicates is probably not the fastest implementation. But it
            # sure is the fastest to implement :P

            oxygen_dimers = set()

            cation_indices = [
                index for index in range(len(self.sites))
                if not self.sites[index].species_string == "O"
            ]

            for index in cation_indices:
                for dimer in self.find_oxygen_dimers(index):
                    oxygen_dimers.add(dimer)

            return list(oxygen_dimers)

        else:

            # Determine the oxygen neighbors for the provided site
            oxygen_neighbors_indices = [
                neighbor["index"] for neighbor
                in self.voronoi.neighbors(site_index, VORONOI_DIST_FACTOR,
                                          VORONOI_ANG_FACTOR)
                if self.sites[neighbor["index"]].species_string == "O"
            ]

            # pdb.set_trace()

            if len(oxygen_neighbors_indices) <= 1:
                raise ValueError("Provided site does not have two oxygen "
                                 "neighbours.\n")

            # Find all oxygen neighbour combinations that can form dimers. This
            # means they are not opposites on the site's octahedron environment.
            oxygen_dimers = []
            oxygen_combinations = itertools.combinations(
                oxygen_neighbors_indices, 2
            )

            for oxygen_pair in oxygen_combinations:

                site = self.sites[site_index]
                oxygen_site_A = self.sites[oxygen_pair[0]]
                oxygen_site_B = self.sites[oxygen_pair[1]]

                oxygen_image_A = site.distance_and_image(oxygen_site_A)[1]
                oxygen_image_B = site.distance_and_image(oxygen_site_B)[1]

                image_A_cart_coords = oxygen_site_A.coords \
                                      + np.dot(oxygen_image_A,
                                               self.lattice.matrix)
                image_B_cart_coords = oxygen_site_B.coords \
                                      + np.dot(oxygen_image_B,
                                               self.lattice.matrix)

                oxygen_vector_A = image_A_cart_coords - site.coords
                oxygen_vector_B = image_B_cart_coords - site.coords

                if angle_between(oxygen_vector_A, oxygen_vector_B) < \
                        oxygen_angle_tol:
                    oxygen_dimers.append(oxygen_pair)

            if len(oxygen_dimers) > 12:
                raise Warning("Found more than 12 oxygen pairs around a single "
                              "site. This can be caused by the use of a small "
                              "unit cell. Results may not be useful.")

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
                        if site.species_string in CATIONS]

        self.remove_cations(remove_sites)

    def find_noneq_dimers(self, site_index=None, method="symmops"):
        """
        A script that distills the non-equivalent oxygen dimers around a site
        index.

        In case no site index is provided, the method will loop over all sites
        which do not correspond to an oxygen.

        Args:
            site_index (int): Index of the site around which the dimers
            should be considered.
            method (str): Method for determining equivalency:

            "symmops" - Two dimers are considered when they are
            symmetrically equivalent, using the symmetry operations of the
            cathode unit cell.

            "representation" - Two dimers are equivalent if their
            environments are the same, i.e. when they have the same
            representation.

        Returns:
            List of Tuples with the dimer indices

        """

        if method == "symmops":

            symmops = SpacegroupAnalyzer(self).get_space_group_operations()

            if site_index is None:

                ineq_cations = self.find_noneq_cations()

                ineq_dimers = []

                for index in ineq_cations:

                    site_dimers = self.find_oxygen_dimers(index)

                    for dimer in site_dimers:

                        d1 = [self.sites[dimer[0]], self.sites[dimer[1]]]

                        inequivalent = True

                        for ineq_dimer in ineq_dimers:

                            d2 = [self.sites[ineq_dimer[0]],
                                  self.sites[ineq_dimer[1]]]

                            if symmops.are_symmetrically_equivalent(d1, d2):

                                inequivalent = False

                        if inequivalent:

                            ineq_dimers.append(dimer)
            else:
                raise NotImplementedError()

            return ineq_dimers

        elif method == "representation":
            # TODO update this method

            raise NotImplementedError()

            # if site_index is None:
            #
            #     ignore = set(CATIONS).union(("O",))
            #
            #     transition_metal_indices = [index for index in
            #                                 range(len(self.sites))
            #                                 if self.sites[index].species_string
            #                                 not in ignore]
            #
            #     for index in transition_metal_indices:
            #
            #         dimers = [Dimer(self, dimer_indices) for dimer_indices
            #                   in self.find_oxygen_dimers(index)]
            #
            #         for dimer in dimers:
            #             if dimer not in noneq_dimers:
            #                 noneq_dimers.append(dimer)
            #
            # else:
            #
            #     dimers = [Dimer(self, dimer_indices) for dimer_indices
            #               in self.find_oxygen_dimers(site_index)]
            #
            #     for dimer in dimers:
            #         if dimer not in noneq_dimers:
            #             noneq_dimers.append(dimer)
            #
            # return noneq_dimers

        else:
            raise IOError("Method for finding non-equivalent dimers is not "
                          "recognized.")

    def list_noneq_dimers(self):
        """
        Create a list of lists of equivalent dimers of the various
        non-equivalent dimers, i.e. group all dimers in the structure in
        lists of dimers that are equivalent to each other.

        Returns:

        """

        symmops = SpacegroupAnalyzer(self).get_space_group_operations()

        dimers = self.find_oxygen_dimers()

        noneq_dimer_lists = [[dimer,] for dimer in self.find_noneq_dimers()]

        for dimer in dimers:

            d1 = [self.sites[dimer[0]], self.sites[dimer[1]]]

            for noneq_dimer_list in noneq_dimer_lists:

                d2 = [self.sites[noneq_dimer_list[0][0]],
                      self.sites[noneq_dimer_list[0][1]]]

                if symmops.are_symmetrically_equivalent(d1, d2):

                    noneq_dimer_list.append(dimer)

        return noneq_dimer_lists


# TODO Currently the whole dimer representation only works for the O-O
# dimers in the O3 stacking. Allowing for different oxygen frameworks will
# require some more possible representations. One way is to figure out the
# structure of the molecule and set up the representation for each
# structure, but I think it would pay off to think if there isn't some
# better way of checking the equivalency of dimer environments. Moreover,
# we should probably first test to see if considering only the immediate
# environment is sufficient.

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
            # TODO Fix issue for small unit cells
            # The issue for small unit cells is that the oxygen atoms have
            # more shared neighbors than two according to the voronoi
            # decomposition, which messes up the assignment of the
            # environment atoms in the representation. This needs to be fixed.
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
                if site.species_string == "O":
                    oxygen_sites.append(site)

            (distance, oxygen_image) = oxygen_sites[0].distance_and_image(
                oxygen_sites[1])

            image_cart_coords = oxygen_sites[1].coords \
                                + np.dot(oxygen_image,
                                         self.cathode.lattice.matrix)

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
            representation = {1: oxy_1.species_and_occu,
                              2: oxy_2.species_and_occu,
                              3: shared_neighbor_3.species_and_occu,
                              4: shared_neighbor_4.species_and_occu}

            # TODO Find a cleaner way of assigning representation positions

            # Loop over the remaining sites to find their representation
            # positions
            for site in dimer_molecule.sites[4:]:

                # Find the sites which are in the plane of the oxygens and
                # their shared neighbors.
                if np.linalg.norm(oxy_1.coords
                                          - (
                                    shared_neighbor_4.coords - oxy_1.coords)
                                          - site.coords) < REPRESENTATION_DIST_TOL:
                    representation[5] = site.species_and_occu

                if np.linalg.norm(oxy_1.coords
                                          - (
                                    shared_neighbor_3.coords - oxy_1.coords)
                                          - site.coords) < REPRESENTATION_DIST_TOL:
                    representation[6] = site.species_and_occu

                if np.linalg.norm(oxy_2.coords
                                          - (
                                    shared_neighbor_4.coords - oxy_2.coords)
                                          - site.coords) < REPRESENTATION_DIST_TOL:
                    representation[7] = site.species_and_occu

                if np.linalg.norm(oxy_2.coords
                                          - (
                                    shared_neighbor_3.coords - oxy_2.coords) \
                                          - site.coords) < REPRESENTATION_DIST_TOL:
                    representation[8] = site.species_and_occu

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
                    representation[9] = site.species_and_occu

                if angle_between(oxy_1_oop, site.coords - oxy_1.coords) \
                        > math.pi - REPRESENTATION_ANGLE_TOL:
                    representation[10] = site.species_and_occu

                if angle_between(oxy_2_oop, site.coords - oxy_2.coords) \
                        < REPRESENTATION_ANGLE_TOL:
                    representation[11] = site.species_and_occu

                if angle_between(oxy_2_oop, site.coords - oxy_2.coords) \
                        > math.pi - REPRESENTATION_ANGLE_TOL:
                    representation[12] = site.species_and_occu

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

            molecule_sites.append(Site(site.species_and_occu,
                                       image_cart_coords))

        return Molecule.from_sites(molecule_sites)

    def visualize_dimer_environment(self, filename=None):
        """
        Creates a .xyz file of the oxygen dimer environment in order to
        visualize it in e.g. VESTA.

        Returns:

        """

        # Turn the structure into a molecule, ignoring the sites with zero
        # occupancy
        dimer_environment_molecule = Molecule.from_sites(
            [site for site in self.get_dimer_molecule().sites
             if not site.species_and_occu == Composition({"": 0})]
        )

        if filename is None:
            filename = str(self.cathode.composition.reduced_composition
                           ).replace(" ", "") + "_" \
                       + str(self.indices[0]) + "_" + str(self.indices[1]) \
                       + ".xyz"

        dimer_environment_molecule.to("xyz", filename)

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

        d = {"cathode": self.cathode.as_dict(),
             "dimer_indices": self.indices}

        return d

    @classmethod
    def from_dict(cls, d):

        return cls(cathode=LiRichCathode.from_dict(d["cathode"]),
                   dimer_indices=d["dimer_indices"])


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

    permutation_numbers = set([i for i in range(1, len(iterable) + 1)])

    if len(permutation_numbers.intersection(permutation)) != len(permutation):
        raise ValueError("Permutation is ill-defined.")

    return [iterable[index - 1] for index in permutation]
