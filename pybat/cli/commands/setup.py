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


def scf(structure_file, functional=("pbe", {}), calculation_dir="",
        write_chgcar=False):
    """
    Set up a standard scf calculation. Always uses the tetrahedron method to
    calculate accurate total energies.

    Args:
        structure_file (str): Path to the Cathode structure file, either
            relative or absolute.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        calculation_dir (str): Path to the directory in which to set up the
            VASP calculation.
        write_chgcar (bool): Write out the charge

    Returns:
        str: Path to the directory in which the calculation is set up.

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

    # Set up the functional
    if functional[0] != "pbe":
        functional_config = _load_yaml_config(functional[0] + "Set")
        functional_config["INCAR"].update(functional[1])
        user_incar_settings.update(functional_config["INCAR"])

    # Set up the calculation directory
    if calculation_dir == "":
        calculation_dir = os.path.join(os.getcwd(), functional[0])
        if functional[0] == "pbeu":
            calculation_dir += "_" + "".join(k + str(functional[1]["LDAUU"][k]) for k
                                             in functional[1]["LDAUU"].keys())
        calculation_dir += "_scf"

    # Set charge density to be written if requested
    if write_chgcar:
        user_incar_settings.update({"LCHARG": True, "LAECHG": True})

        # For the HSE06 calculation, also increase the FFT grids, etc...
        if functional[0] == "hse":
            user_incar_settings.update({"PRECFOCK": "Accurate"})

    # Set up the BulkSCFSet
    scf_calculation = BulkSCFSet(structure=structure,
                                 user_incar_settings=user_incar_settings,
                                 potcar_functional=DFT_FUNCTIONAL)

    # Write the input files to the SCF calculation directory
    scf_calculation.write_input(calculation_dir)
    shutil.copy(structure_file,
                os.path.join(calculation_dir, "initial_cathode.json"))

    return calculation_dir


def relax(structure_file, functional=("pbe", {}), calculation_dir="",
          is_metal=False):
    """
    Set up a standard geometry optimization calculation of a Cathode
    structure. Optimizes both the atomic positions as well as the unit cell.

    Args:
        structure_file (str): Path to the Cathode structure file, either
            relative or absolute.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        calculation_dir (str): Path to the directory in which to set up the
            VASP calculation.
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV.

    Returns:
        str: Path to the directory in which the calculation is set up.

    """

    # Import the structure as a Cathode instance from the structure file
    structure_file = os.path.abspath(structure_file)
    structure = Cathode.from_file(structure_file).as_ordered_structure()

    # Check if a magnetic moment was not provided for the sites. If not, make
    # sure it is zero for the calculations._
    if "magmom" not in structure.site_properties.keys():
        structure.add_site_property("magmom", [0] * len(structure.sites))

    # Set up the calculation
    user_incar_settings = {}

    # Set up the functional
    if functional[0] != "pbe":
        functional_config = _load_yaml_config(functional[0] + "Set")
        functional_config["INCAR"].update(functional[1])
        user_incar_settings.update(functional_config["INCAR"])

    # Set up the calculation directory
    if calculation_dir == "":
        calculation_dir = os.path.join(os.getcwd(), functional[0])
        if functional[0] == "pbeu":
            calculation_dir += "_" + "".join(k + str(functional[1]["LDAUU"][k]) for k
                                             in functional[1]["LDAUU"].keys())
        calculation_dir += "_relax"

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


def transition(directory, functional=("pbe", {}), is_metal=False,
               is_migration=False, optimize_initial=False):
    """
    Set up the geometry optimizations for a transition for a structure, i.e. using
    ISIF = 2. By default, it is assumed that the initial structure is already
    optimized, unless the user specifically requests its optimization.

    If requested, a charge density calculation will be set up for the
    "host structure", i.e. the structure with vacancies at the initial and
    final locations of the migrating ion. This can be used later to provide an
    estimated path for the nudged elastic band calculations. (WIP)

    Args:
        directory (str): Directory in which the transition calculations should be
            set up.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV.
        is_migration (bool): Flag that indicates that the transition is a migration
            of an atom in the structure.
        optimize_initial (bool): Flag that indicates that the initial structure
            should also be optimized.

    Returns:
        None

    """
    # Make sure the directory is written as an absolute path
    directory = os.path.abspath(directory)

    # Obtain the initial and final Cathodes
    (initial_cathode, final_cathode) = find_transition_cathodes(
        directory)

    # Check if a magnetic moment was not provided for the sites
    if "magmom" not in initial_cathode.site_properties.keys():
        # If not, set it to zero for all sites
        initial_cathode.add_site_property("magmom",
                                          [0] * len(initial_cathode.sites))

    if "magmom" not in final_cathode.site_properties.keys():
        # If not, set it to zero for all sites
        final_cathode.add_site_property("magmom",
                                        [0] * len(initial_cathode.sites))

    # Set up the calculations
    user_incar_settings = {"ISIF": 2}

    # Functional
    if functional[0] != "pbe":
        functional_config = _load_yaml_config(functional[0] + "Set")
        functional_config["INCAR"].update(functional[1])
        user_incar_settings.update(functional_config["INCAR"])

    # For metals, add some Methfessel Paxton smearing
    if is_metal:
        user_incar_settings.update({"ISMEAR": 2, "SIGMA": 0.2})

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


def neb(directory, nimages=7, functional=("pbe", {}), is_metal=False,
        is_migration=False):
    """
    Set up the NEB calculation from the initial and final structures.

    Args:
        directory (str): Directory in which the transition calculations should be
            set up.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        nimages (int): Number of images to use in the NEB calculation.
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV.
        is_migration (bool): Flag that indicates that the transition is a migration
            of an atom in the structure.

    Returns:
        None

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
    except BaseException:
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
    # In case an "middle image" has been provided via which to interpolate
    elif os.path.exists(os.path.join(directory, "middle")):
        print("Found a 'middle' directory in the NEB directory. Interpolating "
              "via middle geometry.")
        # Load the middle image
        middle_structure = Structure.from_file(
            os.path.join(directory, "middle", "CONTCAR")
        )
        # Perform an interpolation via this image
        images_1 = initial_structure.interpolate(end_structure=middle_structure,
                                                 nimages=int((nimages + 1) / 2),
                                                 interpolate_lattices=True)
        images_2 = middle_structure.interpolate(end_structure=final_structure,
                                                nimages=int((nimages) / 2 + 1),
                                                interpolate_lattices=True)

        images = images_1[:-1] + images_2

    else:
        # Linearly interpolate the initial and final structures
        images = initial_structure.interpolate(end_structure=final_structure,
                                               nimages=nimages + 1,
                                               interpolate_lattices=True)

    # TODO Add functionality for NEB calculations with changing lattices

    user_incar_settings = {}

    # Set up the functional
    if functional[0] != "pbe":
        functional_config = _load_yaml_config(functional[0] + "Set")
        functional_config["INCAR"].update(functional[1])
        user_incar_settings.update(functional_config["INCAR"])

    # Add the standard Methfessel-Paxton smearing for metals
    if is_metal:
        user_incar_settings.update({"ISMEAR": 2, "SIGMA": 0.2})

    neb_calculation = PybatNEBSet(images, potcar_functional=DFT_FUNCTIONAL,
                                  user_incar_settings=user_incar_settings)

    # Set up the NEB calculation
    neb_calculation.write_input(directory)

    # Make a file to visualize the transition
    neb_calculation.visualize_transition(os.path.join(directory, "transition.cif"))


def dos(structure_file, chgcar_file, functional, kpoint_density):
    """
    Set up Density of States calculation.

    Args:
        structure_file (str): Structure for which the DOS...
        chgcar_file (str):
        functional (tuple):
        kpoint_density (float):

    Returns:
        calculation_dir (str): Directory in which the calculation is set up.

    """
    raise NotImplementedError


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
        directory (str): Path to the directory in which the initial and final
            cathode structure files should be present.
        initial_contains (str): String that is present in the initial Cathode structure
            file.
        final_contains (str): String that is present in the final Cathode structure
            file.

    Returns:
        tuple: Tuple of the initial and final pybat.core.Cathode's

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
        initial_structure (Cathode): Initial cathode structure.
        final_structure (Cathode): Final cathode structure.

    Returns:
        int: Index of the migrating site.

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
