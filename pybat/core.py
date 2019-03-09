# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import itertools
import math
import json
import os
import pdb

import numpy as np

from monty.io import zopen
from monty.json import jsanitize, MSONable
from pymatgen.core import Structure, Composition, Molecule, Site, Element
from pymatgen.analysis.chemenv.coordination_environments.voronoi \
    import DetailedVoronoiContainer
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.transition_state import NEBAnalysis
from pymatgen.util.plotting import pretty_plot
from tabulate import tabulate
from icet.tools.structure_enumeration import enumerate_structures

scipy_old_piecewisepolynomial = True
try:
    from scipy.interpolate import PiecewisePolynomial
except ImportError:
    from scipy.interpolate import CubicSpline

    scipy_old_piecewisepolynomial = False

"""
Module that contains basic classes to represent and calculate the properties of
battery cathodes.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "pre-alpha"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "May 2018"

# TODO Currently, the dimers are defined by their indices.
# This is a consequence of the fact that the DetailedVoronoiContainer
# expects indices for its neighbor method. Frankly, I would prefer sites as
# the basis of the dimer definition as well as it's environment.

# Values for determining the neighbors of a site in a voronoi decomposition
VORONOI_DIST_FACTOR = 1.4
VORONOI_ANG_FACTOR = 0.6

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

    The main idea of this class is to keep track of the original sites of removed
    working ions (e.g. Li, Na, ...) by using sites with empty Compositions. This is
    important to make sure that the voronoi decomposition is successful, and hence
    essential if we want to look at coordinations and neighbors. Another advantage
    is that we can consider the empty cation sites for final positions of transition
    metal migrations.

    By designing a new class, we can update the Structure I/O methods that do not deal
    with sites that have an empty composition well. Moreover, we can design and bundle
    new methods which are useful in the context of battery cathode research.

    However, the class has currently not been fully tested yet for all the methods it
    inherits from Structure, so some of these may produce some unintented results.

    """

    # Tuple of standard working ions for typical battery insertion cathodes.
    # Lawrencium is also in there for the enumerate workaround.
    standard_working_ions = ("Li", "Na", "Lr")

    # Tuple of standard anions for typical battery insertion cathodes.
    standard_anions = ("O", "F")

    def __init__(self, lattice, species, coords, charge=None,
                 validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None):

        super(Cathode, self).__init__(
            lattice=lattice, species=species, coords=coords, charge=charge,
            validate_proximity=validate_proximity, to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties
        )

        self._voronoi = None

    def __str__(self):
        """
        Overwritten string representation, in order to provide information about the
        vacancy sites, as well as the VESTA index, which can be useful when defining
        structural changes.

        Returns:
            (str) String representation of the Cathode.

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

        outs.append(
            tabulate(data,
                     headers=["#", "#VESTA", "SP", "a", "b", "c"] + keys,
                     ))
        return "\n".join(outs)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return self == other.__hash__()

    @property
    def working_ion_configuration(self):
        """
        A list of all sites which correspond to working ions in the Cathode.

        Returns:
            (list): A list of pymatgen.Sites that correspond to the working ions
                in the Cathode.

        """
        return [site for site in self.sites
                if site.species_string in Cathode.standard_working_ions]

    @working_ion_configuration.setter
    def working_ion_configuration(self, configuration):
        """

        Args:
            configuration (dict or list): A dictionary mapping or list of pymatgen.Sites
                that describes the configuration of the working ions in the cathode.

        Returns:
            None

        """
        # TODO Add checks

        # Remove all working ions
        for working_ion in [ion for ion in Cathode.standard_working_ions
                            if Element(ion) in set(self.composition.keys())]:
            self.replace_species({working_ion: {working_ion: 0}})

        # Add the working ion sites
        if isinstance(configuration, dict):
            for working_ion in configuration.keys():
                for index in configuration[working_ion]:
                    self.replace(index, working_ion, properties={"magmom": 0})

        elif all([isinstance(item, Site) for item in configuration]):
            for ion_site in configuration:
                for i, site in enumerate(self):
                    if np.linalg.norm(site.distance(ion_site)) < 0.05:
                        self.replace(i, ion_site.specie,
                                     properties={"magmom": 0})

        else:
            raise TypeError("Working ion configurations should be a dictionary "
                            "mapping working ions to site indices, or a list of "
                            "sites.")

    @property
    def concentration(self):
        """
        The working ion concentration of the cathode, defined versus the pristine,
        i.e. fully discharged structure as a percentage.

        Returns:
            (float): The working ion concentration

        """
        working_ion_sites = [site for site in self.sites if
                             site.species_string in self.standard_working_ions
                             or site.species_and_occu == Composition()]
        return len(self.working_ion_configuration) / len(working_ion_sites)

    @property
    def voronoi(self):
        """
        Pymatgen ChemEnv voronoi decomposition of the cathode structure.

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
            sites: A dictionary mapping or list of pymatgen.Sites that describes the
                configuration of the working ions in the cathode.

        Returns:
            None

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

    def remove_working_ions(self, sites=None):
        """
        Remove working ions from the cathode, i.e. delithiate the structure in
        case Li is the working ion of the cathode.

        Note that this does not remove the sites from the pymatgen Structure.
        The occupancy is simply adjusted to an empty Composition object.

        Args:
            sites: List of indices OR
                List of pymatgen.core.Sites which are to be removed.

        Returns:
            None

        """

        # TODO add checks

        # If no sites are given
        if sites is None:
            # Remove all the working ions
            self.working_ion_configuration = []

        # If a List of integers is given
        elif all([isinstance(item, int) for item in sites]):
            for index in sites:
                self.replace(index, Composition(), properties={"magmom": 0})

        # If a List of sites is given
        elif all([isinstance(item, Site) for item in sites]):
            for site in sites:
                # Check if the provided site corresponds to a working ion site
                if site in self.working_ion_configuration:
                    ion_configuration = self.working_ion_configuration.copy()
                    ion_configuration.remove(site)
                    self.working_ion_configuration = ion_configuration
                else:
                    raise Warning("Requested site not found in working ion "
                                  "configuration.")
        else:
            raise IOError("Incorrect site input.")

    def change_site_distance(self, sites, distance):
        """
        Change the coordinates of two sites in a structure in order to adjust
        their distance.

        Args:
            sites (list): List of two site indices or pymatgen.Sites of
                elements whose distance should be changed.
            distance (float): Final distance between the two sites provided.

        Returns:
            None

        """

        if all(isinstance(el, int) for el in sites):
            site_a = self.sites[sites[0]]
            site_b = self.sites[sites[1]]
        elif all(isinstance(el, Site) for el in sites):
            site_a = sites[0]
            site_b = sites[1]
        else:
            raise IOError("Incorrect input provided.")

        # Find the distance between the sites, as well as the image of site B
        # closest to site A
        (original_distance, closest_image_b) = site_a.distance_and_image(site_b)

        image_cart_coords = self.lattice.get_cartesian_coords(
            site_b.frac_coords + closest_image_b
        )

        # Calculate the vector that connects site A with site B
        connection_vector = image_cart_coords - site_a.coords
        connection_vector /= np.linalg.norm(connection_vector)  # Unit vector

        # Calculate the distance the sites need to be moved.
        site_move_distance = (original_distance - distance) / 2

        # Calculate the new cartesian coordinates of the sites
        new_site_a_coords = site_a.coords + site_move_distance * connection_vector
        new_site_b_coords = site_b.coords - site_move_distance * connection_vector

        # Change the sites in the structure
        self.replace(i=sites[0], species=site_a.species_string,
                     coords=new_site_a_coords,
                     coords_are_cartesian=True,
                     properties=site_a.properties)

        self.replace(i=sites[1], species=site_b.species_string,
                     coords=new_site_b_coords,
                     coords_are_cartesian=True,
                     properties=site_b.properties)

    def update_sites(self, directory, ignore_magmom=False):
        """
        Based on the CONTCAR and OUTCAR of a geometry optimization, update the
        site coordinates and magnetic moments that were optimized. Note that
        this method relies on the cation configuration of the cathode not
        having changed.

        Args:
            directory (str): Directory in which the geometry optimization
                output files (i.e. CONTCAR and OUTCAR) are stores.
            ignore_magmom (bool): Flag that indicates that the final magnetic
                moments of the optimized structure should be ignored. This means
                that the magnetic moments of the Cathode structure will
                remain the same.

        Returns:
            None

        """

        new_cathode = Cathode.from_file(os.path.join(directory, "CONTCAR"))

        out = Outcar(os.path.join(directory, "OUTCAR"))

        if ignore_magmom:
            magmom = [site.properties["magmom"] for site in self.sites
                      if site.species_and_occu != Composition()]
            new_cathode.add_site_property("magmom", magmom)
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

    def find_noneq_cations(self):
        """
        Find a list of the site indices of all non-equivalent cations.

        Returns:
            (list): List of site indices

        """
        symmops = SpacegroupAnalyzer(self).get_space_group_operations()

        cation_indices = [
            index for index in range(len(self.sites))
            if not self.sites[index].species_string not in Cathode.standard_anions
        ]

        # Start with adding the first cation
        inequiv_cations = [cation_indices[0], ]

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

    def get_cation_configurations(self, substitution_sites, cation_list, sizes,
                                  concentration_restrictions=None,
                                  max_configurations=None):
        """
        Get all non-equivalent cation configurations within a specified range of unit
        cell sizes and based on certain restrictions.

        Based on the icet.tools.structure_enumeration.enumerate_structures() method.
        Because there are some issues with this method, related to the allowed
        concentrations, the method will have to be updated later. There is also the fact
        that vacancies can not be inserted in enumerate_structures, which will require
        some workaround using Lawrencium.

        Currently also returns a list of Cathodes, for easy implementation and usage. It
        might be more useful/powerful to design it as a generator later.

        Args:
            substitution_sites (list): List of site indices or pymatgen.Sites to be
                substituted.
            cation_list (list): List of string representations of the cation elements
                which have to be substituted on the substitution sites. Can also
                include "Vac" to introduce vacancy sites.
                E.g. ["Li", "Vac"]; ["Mn", "Co", "Ni"]; ...
            sizes (list): List of unit supercell sizes to be considered for the
                enumeration of the configurations.
                E.g. [1, 2]; range(1, 4); ...
            concentration_restrictions (dict): Dictionary of allowed concentration
                ranges for each element. Note that the concentration is defined
                versus the total amount of atoms in the unit cell.
                E.g. {"Li": (0.2, 0.3)}; {"Ni": (0.1, 0.2, "Mn": (0.05, 0.1)}; ...
            max_configurations (int): Maximum number of configurations to generate.

        Returns:
            (list): List of Cathodes representing different configurations.

        """
        # Check substitution_site input
        if all(isinstance(site, int) for site in substitution_sites):
            substitution_sites = [self.sites[index] for index in substitution_sites]

        # Set up the configuration space
        configuration_space = []
        cation_list = ["Lr" if cat == "Vac" else cat for cat in cation_list]

        for site in self.sites:
            if site in substitution_sites:
                configuration_space.append(cation_list)
            else:
                configuration_space.append([site.species_string, ])

        # TODO adjust this once concentration restrictions work correctly for
        # enumerate_structures
        if concentration_restrictions:
            # Get the largest concentration restriction
            enum_conc_restrictions = [
                {k: v} for k, v in concentration_restrictions.items()
                if v == max(concentration_restrictions.values())][0]
        else:
            concentration_restrictions = {}
            enum_conc_restrictions = None

        configuration_list = []
        configuration_generator = enumerate_structures(
            atoms=AseAtomsAdaptor.get_atoms(self.as_ordered_structure()),
            sizes=sizes,
            chemical_symbols=configuration_space,
            concentration_restrictions=enum_conc_restrictions
        )
        try:
            self.site_properties["magmom"]
        except KeyError:
            print("No magnetic moments found in structure, setting to zero.")
            self.add_site_property("magmom", [0] * len(self))

        for atoms in configuration_generator:

            structure = AseAtomsAdaptor.get_structure(atoms)
            structure.add_site_property(
                "magmom",
                self.site_properties["magmom"] * int(len(structure) / len(self))
            )

            frac_composition = structure.composition.fractional_composition
            elements = [str(el) for el in structure.composition.elements]
            c = concentration_restrictions

            if all([c[el][0] < frac_composition[el] < c[el][1]
                    for el in elements if c.get(el, False)]):
                cathode = Cathode.from_structure(structure.get_sorted_structure())
                cathode.remove_working_ions(
                    [i for i, site in enumerate(cathode)
                     if site.species_string == "Lr"]
                )
                configuration_list.append(cathode)

            if len(configuration_list) == max_configurations:
                break

        return configuration_list

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
        POSCAR file.

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
            return structure.to(fmt, filename, **kwargs)

        else:
            return super(Cathode, self).to(fmt, filename, **kwargs)

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

    # Tolerance for determining whether two oxygens are on opposite sides of an
    # octahedron. If the angle between the two vectors connecting the site and
    # the corresponding oxygens is larger than this value, the oxygens are
    # considered to be opposites.
    oxygen_angle_tol = math.pi * 9 / 10

    def __init__(self, lattice, species, coords, charge=None,
                 validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None):

        super(LiRichCathode, self).__init__(
            lattice=lattice, species=species, coords=coords, charge=charge,
            validate_proximity=validate_proximity, to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties
        )

    def find_oxygen_dimers(self, site_index=None, oxygen_angle_tol=None):
        """
        Returns a list of index pairs corresponding to the oxygen dimers that
        can be formed around the site provided by the user, i.e. with oxygens
        that are neighbours of the provided site.

        Args:
            site_index:
            oxygen_angle_tol:

        Returns:

        """

        if not oxygen_angle_tol:
            oxygen_angle_tol = LiRichCathode.oxygen_angle_tol

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
                oxygen_site_a = self.sites[oxygen_pair[0]]
                oxygen_site_b = self.sites[oxygen_pair[1]]

                oxygen_image_a = site.distance_and_image(oxygen_site_a)[1]
                oxygen_image_b = site.distance_and_image(oxygen_site_b)[1]

                image_a_cart_coords = oxygen_site_a.coords + np.dot(oxygen_image_a,
                                                                    self.lattice.matrix)
                image_b_cart_coords = oxygen_site_b.coords + np.dot(oxygen_image_b,
                                                                    self.lattice.matrix)

                oxygen_vector_a = image_a_cart_coords - site.coords
                oxygen_vector_b = image_b_cart_coords - site.coords

                if angle_between(oxygen_vector_a, oxygen_vector_b) < \
                        oxygen_angle_tol:
                    oxygen_dimers.append(oxygen_pair)

            if len(oxygen_dimers) > 12:
                raise Warning(
                    "Found more than 12 oxygen pairs around a single "
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
                        if site.species_string in Cathode.standard_anions]

        self.remove_working_ions(remove_sites)

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

        ineq_dimers = []

        if method == "symmops":

            symmops = SpacegroupAnalyzer(self).get_space_group_operations()

            # If no site is provided, consider all inequivalent cation sites
            if site_index is None:

                ineq_cations = self.find_noneq_cations()

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

            # Else only consider the site provided
            else:
                site_dimers = self.find_oxygen_dimers(site_index)

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

        noneq_dimer_lists = [[dimer, ] for dimer in self.find_noneq_dimers()]

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
                if np.linalg.norm(
                        oxy_1.coords - (shared_neighbor_4.coords - oxy_1.coords)
                        - site.coords) < REPRESENTATION_DIST_TOL:
                    representation[5] = site.species_and_occu

                if np.linalg.norm(
                        oxy_1.coords - (shared_neighbor_3.coords - oxy_1.coords)
                        - site.coords) < REPRESENTATION_DIST_TOL:
                    representation[6] = site.species_and_occu

                if np.linalg.norm(
                        oxy_2.coords - (shared_neighbor_4.coords - oxy_2.coords)
                        - site.coords) < REPRESENTATION_DIST_TOL:
                    representation[7] = site.species_and_occu

                if np.linalg.norm(
                        oxy_2.coords - (shared_neighbor_3.coords - oxy_2.coords)
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


class DimerNEBAnalysis(NEBAnalysis):
    """
    Subclass of the NEBAnalysis class in order to change the plotting of the
    barriers, as well as allowing for saving the NEB analysis to a json file.
    """

    def __init__(self, r, energies, forces, structures, spline_options=None,
                 dimer_indices=None):
        super().__init__(
            r, energies, forces, structures, spline_options
        )
        self._dimer_indices = tuple(dimer_indices)

    @property
    def dimer_indices(self):
        return self._dimer_indices

    @dimer_indices.setter
    def dimer_indices(self, indices):
        self._dimer_indices = indices

    @property
    def dimer_distances(self):
        return np.array([s.distance_matrix[self.dimer_indices]
                         for s in self.structures])

    def as_dict(self):
        """
        Dict representation of NEBAnalysis.

        Returns:
            JSON serializable dict representation.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                'r': jsanitize(self.r),
                'energies': jsanitize(self.energies),
                'forces': jsanitize(self.forces),
                'structures': [s.as_dict() for s in self.structures],
                "dimer_indices": self.dimer_indices}

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
    def from_dir(cls, root_dir, relaxation_dirs=None, **kwargs):
        """

        Args:
            root_dir:
            relaxation_dirs:
            **kwargs:

        Returns:

        """
        if relaxation_dirs is not None:
            raise NotImplementedError

        indices = tuple(
            [int(el) for el in
             os.path.abspath(root_dir).split('/')[-1].split('_')
             if all([is_number(c) for c in el])]
        )

        neb = super().from_dir(root_dir, relaxation_dirs, **kwargs)

        # Because the dimer indices are based on the internal indices of the
        # Cathode object, we need to load the cathode json files to
        # determine the distance between the dimers properly.
        image_dirs = [file for file in os.listdir(root_dir)
                      if len(file) == 2 and os.path.isdir(file)]

        structures = [
            Cathode.from_file(os.path.join(image_dir, "final_cathode.json"))
            for image_dir in image_dirs
        ]

        # Sort the data according to the directory numbers
        structure_data = sorted(
            zip(image_dirs, structures),
            key=lambda z: int(z[0])
        )

        dimer_neb = DimerNEBAnalysis(
            r=neb.r,
            energies=neb.energies,
            forces=neb.forces,
            structures=[el[1] for el in structure_data],
            spline_options=neb.spline_options,
            dimer_indices=indices,
        )

        return dimer_neb

    @classmethod
    def from_str(cls, input_string, fmt="json"):
        """
        Initialize a DimerNEBAnalysis from a string.

        Currently only supports 'json' formats.

        Args:
            input_string (str): String from which the object is initialized.
            fmt (str): Format of the string representation.

        Returns:
            (*pybat.core.DimerNEBAnalysis*)
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

    @classmethod
    def from_dict(cls, d):

        return cls(r=d['r'], energies=d["energies"], forces=d["forces"],
                   structures=[LiRichCathode.from_dict(structure) for
                               structure in d["structures"]],
                   dimer_indices=d["dimer_indices"])

    def get_plot(self, normalize_rnx_coodinate=True, label_barrier=True):
        """
        Returns the NEB plot. Uses Henkelman's approach of spline fitting
        each section of the reaction path based on tangent force and energies.

        Args:
            label_barrier (bool): Whether to label the maximum barrier.

        Returns:
            matplotlib.pyplot object.
        """
        plt = pretty_plot(12, 8)

        spline_x = np.arange(0, np.max(self.r), 0.01)
        spline_y = self.spline(spline_x) * 1000

        relative_energies = self.energies - self.energies[0]

        plt.plot(spline_x, spline_y, 'k--',
                 self.r, relative_energies * 1000, 'ro',
                 linewidth=2,
                 markersize=10)

        plt.xlabel("O-O Distance ($\mathrm{\AA}$)")
        plt.xticks(self.r[::2], [str(round(d, 2)) for d in
                                 self.dimer_distances[::2]])
        plt.ylabel("Energy (meV)")
        plt.ylim((np.min(spline_y) - 10, np.max(spline_y) * 1.02 + 20))

        if label_barrier:
            data = zip(spline_x, spline_y)
            barrier = max(data, key=lambda d: d[1])
            plt.plot([0, barrier[0]], [barrier[1], barrier[1]], 'k--')
            plt.annotate('%.0f meV' % (np.max(spline_y) - np.min(spline_y)),
                         xy=(barrier[0] / 2, barrier[1] * 1.02),
                         xytext=(barrier[0] / 2, barrier[1] * 1.02),
                         horizontalalignment='center')

        plt.tight_layout()
        return plt


# SO Plagiarism
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


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
