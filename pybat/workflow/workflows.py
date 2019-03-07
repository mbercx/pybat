# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import ast
import numpy as np

from pybat.workflow.firetasks import VaspTask, CustodianTask
from pybat.workflow.fireworks import ScfFirework, RelaxFirework, NebFirework

from pybat.core import Cathode, LiRichCathode, Dimer
from pybat.cli.commands.define import define_dimer, define_migration
from pybat.cli.commands.setup import transition

from ruamel.yaml import YAML
from pymongo.errors import ServerSelectionTimeoutError
from fireworks import Firework, LaunchPad, PyTask, Workflow, FWAction

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


def relax_workflow(structure_file, functional=("pbe", {}), directory="",
                   is_metal=False, in_custodian=False, number_nodes=None):
    """
    Set up a geometry optimization workflow and add it to the launchpad of the
    mongoDB server defined in the config file.

    Args:
        structure_file (str): Path to the geometry file of the structure.
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
    relax_firework = RelaxFirework(structure_file=structure_file,
                                   functional=functional,
                                   directory=directory,
                                   is_metal=is_metal,
                                   in_custodian=in_custodian,
                                   number_nodes=number_nodes)

    # Set up a clear name for the workflow
    cathode = LiRichCathode.from_file(structure_file)
    workflow_name = str(cathode.composition.reduced_formula).replace(" ", "")
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


def configuration_workflow(structure_file, substitution_sites=None, cation_list=None,
                           sizes=None, concentration_restrictions=None,
                           max_configurations=None, functional=("pbe", {}),
                           directory="", in_custodian=False, number_nodes=None):
    # Load the cathode from the structure file
    cat = Cathode.from_file(structure_file)

    # Check for the required input, and request if necessary
    print(cat)
    print()
    if not substitution_sites:
        substitution_sites = [int(i) for i in input(
            "Please provide the substitution site indices, separated by a space: "
        ).split(" ")]
    if not cation_list:
        cation_list = [i for i in input(
            "Please provide the substitution elements, separated by a space: "
        ).split(" ")]
    if not sizes:
        sizes = [int(i) for i in input(
            "Please provide the possible unit cell sizes, separated by a space: "
        ).split(" ")]
    if not concentration_restrictions:
        concentration_restrictions = ast.literal_eval(input(
            "Please provide the concentration restrictions, written as you would "
            "define a dictionary, or None: "))
    if not max_configurations:
        max_configurations = int(input(
            "Please provide the maximum configurations, as an integer: "))
        if max_configurations == 0:
            max_configurations = None

    configurations = cat.get_cation_configurations(
        substitution_sites=substitution_sites,
        cation_list=cation_list,
        sizes=sizes,
        concentration_restrictions=concentration_restrictions,
        max_configurations=max_configurations
    )
    print("Found " + str(len(configurations)) + " configurations.")

    if directory == "":
        directory = os.getcwd()

    functional_dir = functional[0]
    if functional[0] == "pbeu":
        functional_dir += "_" + "".join(k + str(functional[1]["LDAUU"][k]) for k
                                        in functional[1]["LDAUU"].keys())

    firework_list = []
    # TODO add functionality to create new configurations directories if present
    # These scripts do not consider the fact that there already may be configuration
    # directories present. This needs to be changed.

    # Because of the directory structure, we need to differentiate between TM
    # configurations and Li/Vac configurations #TODO
    if "Vac" in cation_list:
        # Set up Li configuration study
        for conf_number, configuration in enumerate(configurations):
            conf_dir = os.path.join(
                os.path.abspath(directory), "tm_conf_1",
                str(round(configuration.concentration, 3)),
                "workion_conf" + str(conf_number), "prim"
            )
            if not os.path.exists(conf_dir):
                os.makedirs(conf_dir)
            configuration.to("json", os.path.join(conf_dir, "cathode.json"))
            relax_dir = os.path.join(conf_dir, functional_dir + "_relax")
            scf_dir = os.path.join(conf_dir, functional_dir + "_scf")

            scf_firework = ScfFirework(
                structure_file=os.path.join(relax_dir, "final_cathode.json"),
                functional=functional,
                directory=scf_dir,
                write_chgcar=False,
                in_custodian=in_custodian,
                number_nodes=number_nodes
            )
            fw_action = FWAction(additions=scf_firework)

            firework_list.append(RelaxFirework(
                structure_file=os.path.join(conf_dir, "cathode.json"),
                functional=functional,
                directory=relax_dir,
                in_custodian=in_custodian,
                number_nodes=number_nodes,
                fw_action=fw_action
            ))
    else:
        # Set up TM configuration study
        for conf_number, configuration in enumerate(configurations):
            try:
                conf_dir = os.path.join(
                    os.path.abspath(directory), "tm_conf_" + str(conf_number),
                    str(round(configuration.concentration, 3)), "workion_conf1", "prim"
                )
            except ZeroDivisionError:
                conf_dir = os.path.join(
                    os.path.abspath(directory), "tm_conf_" + str(conf_number), "prim"
                )
            if not os.path.exists(conf_dir):
                os.makedirs(conf_dir)
            configuration.to("json", os.path.join(conf_dir, "cathode.json"))
            relax_dir = os.path.join(conf_dir, functional_dir + "_relax")
            scf_dir = os.path.join(conf_dir, functional_dir + "_scf")

            scf_firework = ScfFirework(
                structure_file=os.path.join(relax_dir, "final_cathode.json"),
                functional=functional,
                directory=scf_dir,
                write_chgcar=False,
                in_custodian=in_custodian,
                number_nodes=number_nodes
            )
            fw_action = FWAction(additions=scf_firework)

            firework_list.append(RelaxFirework(
                structure_file=os.path.join(conf_dir, "cathode.json"),
                functional=functional,
                directory=relax_dir,
                in_custodian=in_custodian,
                number_nodes=number_nodes,
                fw_action=fw_action
            ))

    # Set up a clear name for the workflow
    workflow_name = str(cat.composition.reduced_formula).replace(" ", "")
    workflow_name += " " + str(cation_list)
    workflow_name += " " + str(functional)

    # Create the workflow
    workflow = Workflow(fireworks=firework_list,
                        name=workflow_name)

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

# region * Token workflows for testing


# endregion
