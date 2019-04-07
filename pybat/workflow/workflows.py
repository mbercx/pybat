# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os

import numpy as np
from fireworks import Firework, LaunchPad, PyTask, Workflow
from pymongo.errors import ServerSelectionTimeoutError
from ruamel.yaml import YAML

from pybat.cli.commands.define import define_dimer, define_migration
from pybat.cli.commands.setup import transition
from pybat.core import Cathode, LiRichCathode, Dimer
from pybat.workflow.firetasks import VaspTask, CustodianTask, ConfigurationTask, \
    EnergyConfTask
from pybat.workflow.fireworks import ScfFirework, RelaxFirework, NebFirework

"""
Package that contains all the Workflows of the pybat package.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2019, Marnik Bercx, University of Antwerp"
__version__ = "alpha"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Mar 2019"

# Load the workflow configuration
CONFIG_FILE = os.path.join(os.path.expanduser("~"), ".pybat_wf_config.yaml")

if os.path.exists(CONFIG_FILE):
    with open(CONFIG_FILE, 'r') as configfile:
        yaml = YAML()
        yaml.default_flow_style = False
        CONFIG = yaml.load(configfile.read())

        try:
            LAUNCHPAD = LaunchPad(
                host=CONFIG["SERVER"].get("host", ""),
                port=int(CONFIG["SERVER"].get("port", 0)),
                name=CONFIG["SERVER"].get("name", ""),
                username=CONFIG["SERVER"].get("username", ""),
                password=CONFIG["SERVER"].get("password", ""),
                ssl=CONFIG["SERVER"].get("ssl", False),
                authsource=CONFIG["SERVER"].get("authsource", None)
            )
        except ServerSelectionTimeoutError:
            raise TimeoutError("Could not connect to server. Please make "
                               "sure the details of the server are correctly "
                               "set up.")

else:
    raise FileNotFoundError("No configuration file found in user's home "
                            "directory. Please use pybat config  "
                            "in order to set up the configuration for "
                            "the workflows.")


# TODO Extend configuration and make the whole configuration setup more user friendly
# Currently the user is not guided to the workflow setup when attempting to use
# pybat workflows, this should change and be tested. Moreover, careful additions should
# be made to make sure all user-specific configuration elements are easily configured
# and implemented in the code.

# TODO Fix the CustodianTask

# TODO Add UnitTests!
# It's really getting time to do this. Think about what unit tests you need and make a
# test suite.

def scf_workflow(structure_file, functional=("pbe", {}), directory="",
                 write_chgcar=False, in_custodian=False, number_nodes=None):
    """
    Set up a self consistent field calculation (SCF) workflow and add it to the
    launchpad of the mongoDB server defined in the config file.

    Args:
        structure_file (str): Path to the geometry file of the structure.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        directory (str): Directory in which the SCF calculation should be performed.
        write_chgcar (bool): Flag that indicates whether the CHGCAR file should
            be written.
        in_custodian (bool): Flag that indicates whether the calculation should be
            run inside a Custodian.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker.

    Returns:
        None

    """

    # Set up the calculation directory
    if directory == "":
        directory = os.path.join(os.getcwd(), functional[0])
        if functional[0] == "pbeu":
            directory += "_" + "".join(k + str(functional[1]["LDAUU"][k]) for k
                                       in functional[1]["LDAUU"].keys())
        directory += "_scf"

    # Set up the SCF Firework
    scf_firework = ScfFirework(
        structure_file=structure_file, functional=functional,
        directory=directory, write_chgcar=write_chgcar,
        in_custodian=in_custodian, number_nodes=number_nodes
    )

    # Set up a clear name for the workflow
    cathode = LiRichCathode.from_file(structure_file)
    workflow_name = str(cathode.composition.reduced_formula).replace(" ", "")
    workflow_name += str(functional)

    # Create the workflow
    workflow = Workflow(fireworks=[scf_firework, ],
                        name=workflow_name)

    LAUNCHPAD.add_wf(workflow)


def relax_workflow(structure, functional=("pbe", {}), directory="",
                   is_metal=False, in_custodian=False, number_nodes=None):
    """
    Set up a geometry optimization workflow and add it to the launchpad of the
    mongoDB server defined in the config file.

    Args:
        structure (pymatgen.Structure): Structure for which to set up the geometry
            optimization workflow.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        directory (str): Directory in which the SCF calculation should be performed.
        is_metal (bool): Flag that indicates whether the material for which the
            geometry optimization should be performed is metallic. Determines the
            smearing method used.
        in_custodian (bool): Flag that indicates wheter the calculation should be
            run inside a Custodian.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker.

    Returns:
        None

    """

    # Set up the calculation directory
    if directory == "":
        directory = os.path.join(os.getcwd(), functional[0])
        if functional[0] == "pbeu":
            directory += "_" + "".join(k + str(functional[1]["LDAUU"][k]) for k
                                       in functional[1]["LDAUU"].keys())
        directory += "_relax"

    # Set up the geometry optimization Firework
    relax_firework = RelaxFirework(structure=structure,
                                   functional=functional,
                                   directory=directory,
                                   is_metal=is_metal,
                                   in_custodian=in_custodian,
                                   number_nodes=number_nodes)

    # Set up a clear name for the workflow
    workflow_name = str(structure.composition.reduced_formula).replace(" ", "")
    workflow_name += str(functional)

    # Create the workflow
    workflow = Workflow(fireworks=[relax_firework, ],
                        name=workflow_name)

    LAUNCHPAD.add_wf(workflow)


def dimer_workflow(structure_file, dimer_indices=(0, 0), distance=0,
                   functional=("pbe", {}), is_metal=False, in_custodian=False,
                   number_nodes=None):
    """
    Set up a workflow that calculates the thermodynamics for a dimer
    formation in the current directory.

    Can later be expanded to also include kinetic barrier calculation.

    Args:
        structure_file (str): Structure file of the cathode material. Note
            that the structure file should be a json format file that is
            derived from the Cathode class, i.e. it should contain the cation
            configuration of the structure.
        dimer_indices (tuple): Indices of the oxygen sites which are to form a
            dimer. If no indices are provided, the user will be prompted.
        distance (float): Final distance between the oxygen atoms. If no
            distance is provided, the user will be prompted.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV. Defaults to False.
        in_custodian (bool): Flag that indicates that the calculations
            should be run within a Custodian. Defaults to False.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker.

    """
    # TODO Change naming scheme

    # Let the user define a dimer, unless one is provided
    dimer_dir = define_dimer(structure_file=structure_file,
                             dimer_indices=dimer_indices,
                             distance=distance,
                             write_cif=True)

    # Set up the FireTask that sets up the transition calculation
    setup_transition = PyTask(
        func="pybat.cli.commands.setup.transition",
        kwargs={"directory": dimer_dir,
                "functional": functional,
                "is_metal": is_metal,
                "is_migration": False}
    )

    # Create the PyTask that runs the calculation
    if in_custodian:
        vasprun = CustodianTask(directory=os.path.join(dimer_dir, "final"))
    else:
        vasprun = VaspTask(directory=os.path.join(dimer_dir, "final"))

    # Extract the final cathode from the geometry optimization
    get_cathode = PyTask(
        func="pybat.cli.commands.get.get_cathode",
        kwargs={"directory": os.path.join(dimer_dir, "final"),
                "write_cif": True}
    )

    # Add number of nodes to spec, or "none"
    firework_spec = {"_launch_dir": os.getcwd()}
    if number_nodes is None:
        firework_spec.update({"_category": "none"})
    else:
        firework_spec.update({"_category": str(number_nodes) + "nodes"})

    transition_firework = Firework(tasks=[setup_transition, vasprun, get_cathode],
                                   name="Dimer Geometry optimization",
                                   spec=firework_spec)

    # Set up the SCF calculation directory
    scf_dir = os.path.join(dimer_dir, "scf_final")

    final_cathode = os.path.join(dimer_dir, "final", "final_cathode.json")

    # Set up the SCF calculation
    scf_firework = ScfFirework(
        structure_file=final_cathode, functional=functional,
        directory=scf_dir, write_chgcar=False, in_custodian=in_custodian,
        number_nodes=number_nodes
    )

    workflow = Workflow(fireworks=[transition_firework, scf_firework],
                        name=structure_file + dimer_dir.split("/")[-1],
                        links_dict={transition_firework: [scf_firework]})

    LAUNCHPAD.add_wf(workflow)


def migration_workflow(structure_file, migration_indices=(0, 0),
                       functional=("pbe", {}), is_metal=False,
                       in_custodian=False, number_nodes=None):
    """
    Set up a workflow that calculates the thermodynamics for a migration in
    the current directory.

    Can later be expanded to also include kinetic barrier calculation.

    Args:
        structure_file (str): Structure file of the cathode material. Note
            that the structure file should be a json format file that is
            derived from the Cathode class, i.e. it should contain the cation
            configuration of the structure.
        migration_indices (tuple): Tuple of the indices which designate the
            migrating site and the vacant site to which the cation will
            migrate. If no indices are provided, the user will be prompted.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV. Defaults to False.
        in_custodian (bool): Flag that indicates that the calculations
            should be run within a Custodian. Defaults to False.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker.

    """

    # TODO Add setup steps to the workflow
    # In case adjustments need to made to the setup of certain calculations,
    #  after which the calculation needs to be rerun, not adding the setup
    # steps to the workflow means that these will have to be rerun manually,
    #  instead of simply relying on the fireworks commands.

    # Let the user define a migration
    migration_dir = define_migration(structure_file=structure_file,
                                     migration_indices=migration_indices,
                                     write_cif=True)

    # Set up the transition calculation
    transition(directory=migration_dir,
               functional=functional,
               is_metal=is_metal,
               is_migration=False)

    # Create the PyTask that runs the calculation
    if in_custodian:
        vasprun = CustodianTask(directory=os.path.join(migration_dir, "final"))
    else:
        vasprun = VaspTask(directory=os.path.join(migration_dir, "final"))

    # Add number of nodes to spec, or "none"
    firework_spec = {"_launch_dir": os.getcwd()}
    if number_nodes is None:
        firework_spec.update({"_category": "none"})
    else:
        firework_spec.update({"_category": str(number_nodes) + "nodes"})

    transition_firework = Firework(tasks=[vasprun],
                                   name="Migration Geometry optimization",
                                   spec=firework_spec)

    workflow = Workflow(fireworks=[transition_firework],
                        name=structure_file + migration_dir.split("/")[-1])

    LAUNCHPAD.add_wf(workflow)


def neb_workflow(directory, nimages=7, functional=("pbe", {}), is_metal=False,
                 is_migration=False, in_custodian=False,
                 number_nodes=None):
    """
    Set up a workflow that calculates the kinetic barrier between two geometries.

    # TODO
    TEMPORARY? Should NEB be integrated in other workflows? If so, should we still
    have a separate NEB workflow?

    Args:
        directory (str): Directory in which the NEB calculation should be performed.
        nimages (int): Number of images to use for the NEB calculation.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV. Defaults to False.
        is_migration (bool): Flag that indicates that the transition is a migration
            of an atom in the structure.
        in_custodian (bool): Flag that indicates that the calculations
            should be run within a Custodian. Defaults to False.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker. Defaults to the number of images.

    """
    # If no number of nodes is specified, take the number of images
    if number_nodes is None:
        number_nodes = nimages

    # Create the Firework that sets up and runs the NEB
    neb_firework = NebFirework(
        directory=directory,
        nimages=nimages,
        functional=functional,
        is_metal=is_metal,
        is_migration=is_migration,
        in_custodian=in_custodian,
        number_nodes=number_nodes
    )

    # Add number of nodes to spec, or "none"
    firework_spec = {"_launch_dir": os.getcwd()}
    if number_nodes is None:
        firework_spec.update({"_category": "none"})
    else:
        firework_spec.update({"_category": str(number_nodes) + "nodes"})

    cathode = Cathode.from_file(
        os.path.join(directory, "final", "initial_cathode.json")
    )
    dir_name = os.path.abspath(directory).split("/")[-1]
    workflow_name = str(cathode.composition).replace(" ", "") + " " + dir_name

    workflow = Workflow(fireworks=[neb_firework, ],
                        name=workflow_name)

    LAUNCHPAD.add_wf(workflow)


def configuration_workflow(structure_file, substitution_sites=None, element_list=None,
                           sizes=None, concentration_restrictions=None,
                           max_configurations=None, functional=("pbe", {}),
                           directory=None, in_custodian=False, number_nodes=None):
    """
    Set up a workflow for a set of atomic configurations, which includes a geometric
    optimization as well as a SCF calculation based on the final geometry.

    Args:
        structure_file (str): Structure file of the cathode material. Note
            that the structure file should be a json format file that is
            derived from the Cathode class, i.e. it should contain the cation
            configuration of the structure.
        substitution_sites (list): List of site indices or pymatgen.Sites to be
            substituted.
        element_list (list): List of string representations of the cation elements
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
        max_configurations (int): Maximum number of new configurations to generate.
            Note that the function detects all the cathode.json files present
            in the directory tree and ignores the corresponding configurations.
            max_configurations is the maximum number of new configurations that need
            to be generated, i.e. on top of the configurations already present in the
            directory tree in the form of cathode.json files.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        directory (str): Path to the directory in which the configurations and
            calculations should be set up.
        in_custodian (bool): Flag that indicates that the calculations
            should be run within a Custodian. Defaults to False.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker.

    Returns:
        None

    """

    # Load the cathode from the structure file
    cathode = Cathode.from_file(structure_file)

    # Check for the required input, and request if necessary
    if not substitution_sites or not element_list or not sizes:
        print(cathode)
        print()
    if not substitution_sites:
        substitution_sites = [int(i) for i in input(
            "Please provide the substitution site indices, separated by a space: "
        ).split(" ")]
    if not element_list:
        element_list = [i for i in input(
            "Please provide the substitution elements, separated by a space: "
        ).split(" ")]
    if not sizes:
        sizes = [int(i) for i in input(
            "Please provide the possible unit cell sizes, separated by a space: "
        ).split(" ")]

    # Set up the directory
    if directory == "":
        directory = os.getcwd()
    directory = os.path.abspath(directory)

    configuration_task = ConfigurationTask(
        structure=cathode,
        directory=directory,
        substitution_sites=list(substitution_sites),
        element_list=element_list,
        sizes=list(sizes),
        concentration_restrictions=concentration_restrictions,
        max_configurations=max_configurations
    )

    energy_task = EnergyConfTask(
        functional=functional,
        in_custodian=in_custodian,
        number_nodes=number_nodes
    )

    # Set up a (sort of) clear name for the workflow
    workflow_name = str(cathode.composition.reduced_formula).replace(" ", "")
    workflow_name += " " + str(element_list)
    workflow_name += " " + str(functional)

    configuration_fw = Firework(tasks=[configuration_task, energy_task],
                                name="Configuration Setup",
                                spec={"_category": "none"})

    # Create the workflow
    workflow = Workflow(
        fireworks=[configuration_fw],
        name=workflow_name
    )

    LAUNCHPAD.add_wf(workflow)


def noneq_dimers_workflow(structure_file, distance, functional=("pbe", {}),
                          is_metal=False, in_custodian=False, number_nodes=None):
    """
    Run dimer calculations for all the nonequivalent dimers in a structure.

    Args:
        structure_file (str): Structure file of the cathode material. Note
            that the structure file should be a json format file that is
            derived from the Cathode class, i.e. it should contain the cation
            configuration of the structure.
        distance (float): Final distance between the oxygen atoms. If no
            distance is provided, the user will be prompted.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV. Defaults to False.
        in_custodian (bool): Flag that indicates that the calculations
            should be run within a Custodian. Defaults to False.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker.

    Returns:
        None

    """

    lirich = LiRichCathode.from_file(structure_file)
    dimer_lists = lirich.list_noneq_dimers()

    for dimer_list in dimer_lists:

        # Find the dimer closest to the center of the lattice. Just for
        # visualization purposes.
        central_dimer = [(), 1e10]

        for dimer in dimer_list:

            dimer_center = Dimer(lirich, dimer).center
            lattice_center = np.sum(lirich.lattice.matrix, 0) / 3

            dist_to_center = np.linalg.norm(dimer_center - lattice_center)

            if dist_to_center < central_dimer[1]:
                central_dimer = [dimer, dist_to_center]

        dimer_workflow(structure_file=structure_file,
                       dimer_indices=central_dimer[0],
                       distance=distance,
                       functional=functional,
                       is_metal=is_metal,
                       in_custodian=in_custodian,
                       number_nodes=number_nodes)


def site_dimers_workflow(structure_file, site_index, distance,
                         functional=("pbe", {}), is_metal=False,
                         in_custodian=False, number_nodes=None):
    """
    Run dimer calculations for all the dimers around a site.

    Args:
        structure_file (str): Structure file of the cathode material. Note
            that the structure file should be a json format file that is
            derived from the Cathode class, i.e. it should contain the cation
            configuration of the structure.
        site_index (int): Index of the site around which the dimers should
            be investigated. Corresponds to the internal Python index.
        distance (float): Final distance between the oxygen atoms. If no
            distance is provided, the user will be prompted.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV. Defaults to False.
        in_custodian (bool): Flag that indicates that the calculations
            should be run within a Custodian. Defaults to False.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker.

    Returns:
        None

    """

    lirich = LiRichCathode.from_file(structure_file)
    dimer_list = lirich.find_noneq_dimers(int(site_index))

    for dimer in dimer_list:
        dimer_workflow(structure_file=structure_file,
                       dimer_indices=dimer,
                       distance=distance,
                       functional=functional,
                       is_metal=is_metal,
                       in_custodian=in_custodian,
                       number_nodes=number_nodes)


# region * Utility scripts

def find_all(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            result.append(os.path.join(root, name))
    return result


def find_all_cathode_hashes(path):
    return [Cathode.from_file(file).__hash__() for file in find_all("cathode.json", path)]


def find_hash_dict(path):
    path = os.path.abspath(path)
    return {Cathode.from_file(file).__hash__(): file.replace(path, "").replace(
        "cathode.json", "")
        for file in find_all("cathode.json", path)}


def generate_conf_dir(directory, element_list, configuration, number):
    if "Vac" in element_list:
        # Set up Li configuration directory
        conf_dir = os.path.join(
            directory, "tm_conf_1",
            str(round(configuration.concentration, 3)),
            "workion_conf" + str(number), "prim"
        )
    else:
        # Set up TM configuration directory
        try:
            conf_dir = os.path.join(
                directory, "tm_conf_" + str(number),
                str(round(configuration.concentration, 3)), "workion_conf1",
                "prim"
            )
        except ZeroDivisionError:
            conf_dir = os.path.join(
                directory, "tm_conf_" + str(number), "prim"
            )

    return conf_dir

# endregion

# region * Token workflows for testing


# endregion
