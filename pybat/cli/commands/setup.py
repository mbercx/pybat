# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import numpy as np
import os
import shutil
import pdb

from pybat.core import Cathode, LiRichCathode
from pybat.sets import BulkSCFSet, BulkRelaxSet, PybatNEBSet
from monty.serialization import loadfn
from pymatgen.core import Structure
from pymatgen.analysis.path_finder import ChgcarPotential, NEBPathfinder
from pymatgen.io.vasp.outputs import Chgcar, Outcar
from pymatgen.io.vasp.sets import MPStaticSet

"""
Setup scripts for the calculations.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Jul 2018"

MODULE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "../../set_configs")

DFT_FUNCTIONAL = "PBE_54"


def _load_yaml_config(filename):
    config = loadfn(os.path.join(MODULE_DIR, "%s.yaml" % filename))
    return config


def scf(structure_file, calculation_dir="", write_chgcar=False,
        dftu_values=None, hse_calculation=False):
    """
    Set up a standard scf calculation. Always uses the tetrahedron method to
    calculate accurate total energies.

    Args:
        structure_file (str): Path to the Cathode structure file, either
            relative or absolute.
        write_chgcar (bool): Write out the charge
        calculation_dir (str): Path to the directory in which to set up the
            VASP calculation.
        hse_calculation (bool): Flag that indicates that a hybrid HSE06
            functional should be used for the geometry optimization.
    """

    # Import the structure as a Cathode instance from the structure file
    structure_file = os.path.abspath(structure_file)
    structure = Cathode.from_file(structure_file).as_ordered_structure()

    # Check if a magnetic moment was not provided for the sites. If not, make
    # sure it is zero for the calculations.
    if "magmom" not in structure.site_properties.keys():
        structure.add_site_property("magmom", [0] * len(structure.sites))

    # Set up the calculation
    user_incar_settings = {}

    if hse_calculation:

        hse_config = _load_yaml_config("HSESet")
        user_incar_settings.update(hse_config["INCAR"])

        if calculation_dir == "":
            # Set up the calculation directory
            current_dir = os.path.dirname(".")
            calculation_dir = os.path.join(current_dir, "hse_scf")

    else:

        dftu_config = _load_yaml_config("DFTUSet")
        user_incar_settings.update(dftu_config["INCAR"])

        if not dftu_values is None:
            user_incar_settings.update({"LDAUU": dftu_values})

        if calculation_dir == "":
            # Set up the calculation directory
            current_dir = os.path.dirname(".")
            calculation_dir = os.path.join(current_dir, "dftu_scf")

    # Set charge density to be written if requested
    if write_chgcar:
        user_incar_settings.update({"LCHARG": True, "LAECHG": True})

        # For the HSE06 calculation, also increase the FFT grids, etc...
        if hse_calculation:
            user_incar_settings.update({"PRECFOCK": "Accurate"})

    # Set up the geometry optimization
    scf_calculation = BulkSCFSet(structure=structure,
                                 user_incar_settings=user_incar_settings,
                                 potcar_functional=DFT_FUNCTIONAL)

    # Write the input files to the geometry optimization directory
    scf_calculation.write_input(calculation_dir)
    shutil.copy(structure_file,
                os.path.join(calculation_dir, "initial_cathode.json"))

    return calculation_dir


def relax(structure_file, calculation_dir="", is_metal=False, dftu_values=None,
          hse_calculation=False):
    """
    Set up a standard geometry optimization calculation of a Cathode
    structure. Optimizes both the atomic positions as well as the unit cell.

    Args:
        structure_file (str): Path to the Cathode structure file, either
            relative or absolute.
        calculation_dir (str): Path to the directory in which to set up the
            VASP calculation.
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV.
        dftu_values (dict): Dictionary of LDAUU values, e.g. either
            {"LDAUU":{"O":{"Fe":5}}} to set LDAUU for Fe to 5 in an oxide,
            or {"LDAUU":{"Fe":5}} to set LDAUU to 5 regardless of the input
            structure.
        hse_calculation (bool): Flag that indicates that a hybrid HSE06
            functional should be used for the geometry optimization.
    """

    # Import the structure as a Cathode instance from the structure file
    structure_file = os.path.abspath(structure_file)
    structure = Cathode.from_file(structure_file).as_ordered_structure()

    # Check if a magnetic moment was not provided for the sites. If not, make
    # sure it is zero for the calculations.
    if "magmom" not in structure.site_properties.keys():
        structure.add_site_property("magmom", [0] * len(structure.sites))

    # Set up the calculation

    user_incar_settings = {}

    if hse_calculation:

        hse_config = _load_yaml_config("HSESet")
        user_incar_settings.update(hse_config["INCAR"])

        if calculation_dir == "":
            # Set up the calculation directory
            current_dir = os.path.dirname(".")
            calculation_dir = os.path.join(current_dir, "hse_relax")

    else:

        dftu_config = _load_yaml_config("DFTUSet")
        user_incar_settings.update(dftu_config["INCAR"])

        if not dftu_values is None:
            user_incar_settings.update({"LDAUU": dftu_values})

        if calculation_dir == "":
            # Set up the calculation directory
            current_dir = os.path.dirname(".")
            calculation_dir = os.path.join(current_dir, "dftu_relax")

    # For metals, add some Methfessel Paxton smearing
    if is_metal:
        user_incar_settings.update({"ISMEAR": 2, "SIGMA": 0.2})

    # Set up the geometry optimization
    geo_optimization = BulkRelaxSet(structure=structure,
                                    user_incar_settings=user_incar_settings,
                                    potcar_functional=DFT_FUNCTIONAL)

    # Write the input files to the geometry optimization directory
    geo_optimization.write_input(calculation_dir)
    shutil.copy(structure_file,
                os.path.join(calculation_dir, "initial_cathode.json"))

    return calculation_dir


def transition(directory, is_metal=False, is_migration=False,
               dftu_values=None, hse_calculation=False,
               optimize_initial=False):
    """
    This script will set up the geometry optimizations for a transition
    structure, i.e. using ISIF = 2. It is assumed that the initial structure
    is already optimized, unless the user specifically requests its
    optimization.

    If requested, a charge density calculation will be set up for the
    "host structure", i.e. the structure with vacancies at the initial and
    final locations of the migrating ion. This is used later to provide an
    estimated path for the nudged elastic band calculations.

    """
    # Make sure the directory is written as an absolute path
    directory = os.path.abspath(directory)

    # Obtain the initial and final Cathodes
    (initial_cathode, final_cathode) = find_transition_cathodes(
        directory)

    # Check if a magnetic moment was not provided for the sites. If not, make
    # sure it is zero for the calculations.
    if "magmom" not in initial_cathode.site_properties.keys():
        initial_cathode.add_site_property("magmom",
                                          [0] * len(initial_cathode.sites))

    if "magmom" not in final_cathode.site_properties.keys():
        final_cathode.add_site_property("magmom",
                                        [0] * len(initial_cathode.sites))

    # Set up the calculations

    user_incar_settings = {"ISIF": 2}

    # For metals, add some Methfessel Paxton smearing
    if is_metal:
        user_incar_settings.update({"ISMEAR": 2, "SIGMA": 0.2})

    # Load the correct INCAR settings for the chosen functional
    if hse_calculation:

        hse_config = _load_yaml_config("HSESet")
        user_incar_settings.update(hse_config["INCAR"])

    else:

        dftu_config = _load_yaml_config("DFTUSet")
        user_incar_settings.update(dftu_config["INCAR"])

        if not dftu_values is None:
            user_incar_settings.update({"LDAUU": dftu_values})

    # If requested, set up the initial structure optimization calculation
    if optimize_initial:
        initial_optimization = BulkRelaxSet(
            structure=initial_cathode.as_ordered_structure(),
            potcar_functional=DFT_FUNCTIONAL,
            user_incar_settings=user_incar_settings
        )
        initial_optimization.write_input(os.path.join(directory, "initial"))
        initial_cathode.to("json", os.path.join(directory, "initial",
                                                "initial_cathode.json"))
    else:
        os.makedirs(os.path.join(directory, "initial"), exist_ok=True)
        initial_cathode.to("json", os.path.join(directory, "initial",
                                                "final_cathode.json"))

    # Set up the final structure optimization calculation
    final_optimization = BulkRelaxSet(
        structure=final_cathode.as_ordered_structure(),
        potcar_functional=DFT_FUNCTIONAL,
        user_incar_settings=user_incar_settings
    )

    final_optimization.write_input(os.path.join(directory, "final"))

    # Write the initial structure of the final Cathode to the optimization
    # directory
    final_cathode.to("json", os.path.join(directory, "final",
                                          "initial_cathode.json"))

    # If the transition is a migration of an atom in the structure, set up the
    # calculation for the charge density, used to find a better initial pathway
    if is_migration:
        migration_site_index = find_migrating_ion(initial_cathode,
                                                  final_cathode)

        host_structure = initial_cathode.copy()
        host_structure.remove_sites([migration_site_index])
        host_scf = MPStaticSet(host_structure,
                               potcar_functional=DFT_FUNCTIONAL)
        host_scf.write_input(os.path.join(directory, "host"))


def neb(directory, nimages=8, is_metal=False, is_migration=False,
        hse_calculation=False):
    """
    Set up the NEB calculation from the initial and final structures.

    """
    directory = os.path.abspath(directory)

    # Extract the optimized initial and final geometries
    initial_dir = os.path.join(directory, "initial")
    final_dir = os.path.join(directory, "final")

    try:
        # Check to see if the initial final_cathode structure is present
        initial_structure = Cathode.from_file(
            os.path.join(initial_dir, "final_cathode.json")
        ).as_ordered_structure()

    except FileNotFoundError:
        # In case the required json file is not present, check to see if
        # there is VASP output which can be used
        initial_structure = Structure.from_file(os.path.join(initial_dir,
                                                             "CONTCAR"))

        # Add the magnetic configuration to the initial structure
        initial_out = Outcar(os.path.join(initial_dir, "OUTCAR"))
        initial_magmom = [site["tot"] for site in initial_out.magnetization]

        try:
            initial_structure.add_site_property("magmom", initial_magmom)
        except ValueError:
            if len(initial_magmom) == 0:
                print("No magnetic moments found in OUTCAR file. Setting "
                      "magnetic moments to zero.")
                initial_magmom = [0] * len(initial_structure)
                initial_structure.add_site_property("magmom", initial_magmom)
            else:
                raise ValueError("Number of magnetic moments in OUTCAR file "
                                 "do not match the number of sites!")
    except:
        raise FileNotFoundError("Could not find required structure "
                                "information in " + initial_dir + ".")

    final_structure = Structure.from_file(os.path.join(final_dir,
                                                       "CONTCAR"))

    # In case the transition is a migration
    if is_migration:
        # Set up the static potential for the Pathfinder from the host charge
        # density
        host_charge_density = Chgcar.from_file(os.path.join(directory, "host"))
        host_potential = ChgcarPotential(host_charge_density)

        migration_site_index = find_migrating_ion(initial_structure,
                                                  final_structure)

        neb_path = NEBPathfinder(start_struct=initial_structure,
                                 end_struct=final_structure,
                                 relax_sites=migration_site_index,
                                 v=host_potential)

        images = neb_path.images
        neb_path.plot_images("neb.vasp")

    else:
        # Linearly interpolate the initial and final structures
        images = initial_structure.interpolate(end_structure=final_structure,
                                               nimages=nimages + 1,
                                               interpolate_lattices=True)

    # TODO Add functionality for NEB calculations with changing lattices

    user_incar_settings = {}

    # Add the standard Methfessel-Paxton smearing for metals
    if is_metal:
        user_incar_settings.update({"ISMEAR": 2, "SIGMA": 0.2})

    neb_calculation = PybatNEBSet(images, potcar_functional=DFT_FUNCTIONAL,
                                  user_incar_settings=user_incar_settings)

    # Set up the NEB calculation
    neb_calculation.write_input(directory)

    # Make a file to visualize the transition
    neb_calculation.visualize_transition()


def dimers(structure_file, dimer_distance=1.4,
           is_metal=False, hse_calculation=False):
    """
    Set up the geometric optimizations for the non-equivalent dimer formations
    in a Li-Rich cathode material.

    Args:
        structure_file (str): Structure file of the cathode material. Note
            that the structure file should be a json format file that is
            derived from the Cathode class, i.e. it should contain the cation
            configuration of the structure.
        dimer_distance (float): Final distance in angstroms between the oxygen
            atoms in oxygen pair.
        hse_calculation (bool): Flag that indicates that the hybrid functional
            HSE06 should be used to calculate the exchange-correlation energy.

    """
    raise NotImplementedError

    # TODO FIX

    # Load the cathode structure
    cathode = LiRichCathode.from_file(structure_file)

    # Find the non-equivalent dimers
    dimers = cathode.find_noneq_dimers()

    # Set up the calculations
    user_incar_settings = {"ISIF": 2}

    # Add the standard Methfessel-Paxton smearing for metals
    if is_metal:
        user_incar_settings.update({"ISMEAR": 2, "SIGMA": 0.2})

    # Load the correct INCAR settings for the chosen functional
    if hse_calculation:

        hse_config = _load_yaml_config("HSESet")
        user_incar_settings.update(hse_config["INCAR"])

    else:

        dftu_config = _load_yaml_config("DFTUSet")
        user_incar_settings.update(dftu_config["INCAR"])

    # Set up the geometry optimization for the initial structure
    try:
        os.mkdir("initial")
    except FileExistsError:
        pass

    initial_optimization = BulkRelaxSet(
        structure=cathode.as_structure(),
        potcar_functional=DFT_FUNCTIONAL,
        user_incar_settings=user_incar_settings
    )

    initial_optimization.write_input("initial")

    # Set up the geometry optimization calculations for the various dimers
    for dimer in dimers:

        # Set up the dimer directory
        dimer_directory = "".join(str(dimer.indices[0]) + "_"
                                  + str(dimer.indices[1]))
        final_dir = os.path.join(dimer_directory, "final")

        try:
            os.mkdir(dimer_directory)
            os.mkdir(final_dir)
        except FileExistsError:
            pass

        # Write the molecule representing the dimer to the dimer directory
        dimer.visualize_dimer_environment(os.path.join(dimer_directory,
                                                       "dimer.xyz"))

        # Set up the dimer structure, i.e. move the oxygen pair closer together
        dimer_structure = cathode.copy()
        dimer_structure.change_site_distance(site_indices=dimer.indices,
                                             distance=dimer_distance)
        if hse_calculation:
            raise NotImplementedError("HSE06 calculation not implemented yet.")

        # Set up the geometry optimization for the dimer structure
        dimer_optimization = BulkRelaxSet(
            structure=dimer_structure.as_structure(),
            potcar_functional=DFT_FUNCTIONAL,
            user_incar_settings=user_incar_settings
        )

        # Write the calculation files to the 'final' directory
        dimer_optimization.write_input(final_dir)


###########
# UTILITY #
###########


def find_transition_cathodes(directory, initial_contains="init.json",
                             final_contains="final.json"):
    """
    Find the initial and final Cathodes for a transition from the files in a
    directory.

    The function demands .json type files, in order to include the magnetic
    moments, as well as the vacant Sites in the Cathode.

    Args:
        directory:
        initial_contains:
        final_contains:

    Returns:
        Tuple of the initial and final pybat.core.Cathode's

    """
    directory = os.path.abspath(directory)

    initial_structure_file = None
    final_structure_file = None

    for item in os.listdir(directory):

        if initial_contains in item \
                and os.path.isfile(os.path.join(directory, item)):
            initial_structure_file = os.path.join(directory, item)

        if final_contains in item \
                and os.path.isfile(os.path.join(directory, item)):
            final_structure_file = os.path.join(directory, item)

    if initial_structure_file:
        initial_cathode = Cathode.from_file(initial_structure_file)
    else:
        raise FileNotFoundError("No suitably named initial structure file in "
                                "directory.")

    if final_structure_file:
        final_cathode = Cathode.from_file(final_structure_file)
    else:
        raise FileNotFoundError("No suitably named final structure file in "
                                "directory.")

    return initial_cathode, final_cathode


def find_migrating_ion(initial_structure, final_structure):
    """
    Tries to find the index of the ion in the structure that is migrating, by
    considering the distance that the sites have moved. The site whose
    coordinates have changed the most is considered to be the migrating ion.

    Note that the site indices of the initial and final structure have to
    correspond for the algorithm to work.

    Args:
        initial_structure:
        final_structure:

    Returns:
        (int) Index of the migrating site

    """

    # Check if both structures have the same number of sites.
    if len(initial_structure.sites) != len(final_structure.sites):
        raise IOError("Provided structures do not have the same number of "
                      "atoms.")

    max_distance = 0
    migrating_site = None

    # TODO Build in some checks, i.e. make sure that the other ions have not
    # moved significantly, and that the migrating ion has moved sufficiently.

    for site_pair in zip(initial_structure.sites, final_structure.sites):

        # TODO This definition of distance is inadequate, i.e. if the site has
        # crossed into another unit cell, the distance will be very big. This
        # can be fixed by using the nearest image!
        distance = np.linalg.norm(site_pair[0].coords - site_pair[1].coords)

        if distance > max_distance:
            max_distance = distance
            migrating_site = site_pair[0]

    return initial_structure.sites.index(migrating_site)
