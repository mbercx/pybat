# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os

import numpy as np
from fireworks import Firework, PyTask, Workflow

from pybat.cli.commands.define import define_migration
from pybat.cli.commands.setup import transition
from pybat.core import Cathode, Dimer
from pybat.workflow.firetasks import VaspTask, CustodianTask, ConfigurationTask, \
    EnergyConfTask
from pybat.workflow.fireworks import PybatStaticFW, PybatOptimizeFW, NebFirework

"""
Package that contains all the Workflows of the pybat package.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2019, Marnik Bercx, University of Antwerp"
__version__ = "alpha"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Mar 2019"


# TODO Fix the CustodianTask

# TODO Add UnitTests!
# It's really getting time to do this. Think about what unit tests you need and make a
# test suite.

def get_wf_static(structure, directory, functional=("pbe", {}),
                  write_chgcar=False, in_custodian=False, number_nodes=None):
    """
    Set up a workflow for a standard static calculation.

    Args:
        structure (pymatgen.Structure): Structure for which to set up the static
            calculation workflow.
        directory (str): Directory in which the static calculation should be performed.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        write_chgcar (bool): Flag that indicates whether the CHGCAR file should
            be written.
        in_custodian (bool): Flag that indicates whether the calculation should be
            run inside a Custodian.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker.

    Returns:
        fireworks.Workflow

    """

    # Set up the static Firework
    static_fw = PybatStaticFW(
        structure=structure, functional=functional,
        directory=directory, write_chgcar=write_chgcar,
        in_custodian=in_custodian, number_nodes=number_nodes
    )

    # Set up a clear name for the workflow
    workflow_name = str(structure.composition.reduced_formula).replace(" ", "")
    workflow_name += " " + str(functional)

    return Workflow(fireworks=[static_fw, ],
                    name=workflow_name)


def get_wf_optimize(structure, directory, functional=("pbe", {}),
                    is_metal=False, in_custodian=False, number_nodes=None,
                    fw_action=None):
    """
    Set up a geometry optimization workflow.

    Args:
        structure (pymatgen.Structure): Structure for which to set up the geometry
            optimization workflow.
        directory (str): Directory in which the geometry optimization should be performed.
        functional (tuple): Tuple with the functional details. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        is_metal (bool): Flag that indicates whether the material for which the
            geometry optimization should be performed is metallic. Determines the
            smearing method used.
        in_custodian (bool): Flag that indicates wheter the calculation should be
            run inside a Custodian.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker.
        fw_action (fireworks.FWAction): FWAction to return after the final
                PulayTask is completed.

    Returns:
        None

    """

    # Set up the geometry optimization Firework
    optimize_fw = PybatOptimizeFW(structure=structure,
                                  functional=functional,
                                  directory=directory,
                                  is_metal=is_metal,
                                  in_custodian=in_custodian,
                                  number_nodes=number_nodes,
                                  fw_action=fw_action)

    # Set up a clear name for the workflow
    workflow_name = str(structure.composition.reduced_formula).replace(" ", "")
    workflow_name += " " + str(functional)

    # Create the workflow
    return Workflow(fireworks=[optimize_fw, ],
                    name=workflow_name)


def get_wf_configurations(structure, directory, substitution_sites=None,
                          element_list=None, sizes=None, concentration_restrictions=None,
                          max_configurations=None, configuration_list=None,
                          include_existing=True, functional=("pbe", {}),
                          in_custodian=False, number_nodes=None):
    """
    Set up a workflow for a set of atomic configurations, which includes a geometric
    optimization as well as a static calculation based on the final geometry.

    Args:
        structure (pymatgen.Structure): Structure for which to set up the configuration
            Workflow.
        directory (str): Path to the directory in which the configurations and
            calculations should be set up.
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
        max_configurations (int): Maximum number of configurations for which the total
            energy should be calculated. Note that that in case include_existing is
            set to False, the existing configurations in the directory tree are not
            included in this number.
        functional (tuple): Tuple with the functional details. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        include_existing (bool): Include the existing configurations in the directory
            tree for the calculations.
        in_custodian (bool): Flag that indicates that the calculations
            should be run within a Custodian. Defaults to False.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker.

    Returns:
        None

    """
    substitution_sites = substitution_sites or []
    sizes = sizes or []

    configuration_task = ConfigurationTask(
        structure=structure,
        directory=directory,
        substitution_sites=list(substitution_sites),
        element_list=element_list,
        sizes=list(sizes),
        concentration_restrictions=concentration_restrictions,
        max_configurations=max_configurations,
        include_existing=include_existing,
        configuration_list=configuration_list
    )

    energy_task = EnergyConfTask(
        functional=functional,
        in_custodian=in_custodian,
        number_nodes=number_nodes
    )

    # Set up a (sort of) clear name for the workflow
    workflow_name = str(structure.composition.reduced_formula).replace(" ", "")
    workflow_name += " " + str(element_list)
    workflow_name += " " + str(functional)

    configuration_fw = Firework(tasks=[configuration_task, energy_task],
                                name="Configuration Setup",
                                spec={"_category": "none"})

    # Create the workflow
    return Workflow(
        fireworks=[configuration_fw],
        name=workflow_name
    )


def get_wf_migration(structure, migration_indices=(0, 0),
                     functional=("pbe", {}), is_metal=False,
                     in_custodian=False, number_nodes=None):
    """
    Set up a workflow that calculates the reaction energy for a migration in
    the current directory.

    Can later be expanded to also include kinetic barrier calculation.

    Args:
        structure (pybat.core.Cathode): Cathode for which to calculate the reaction
            energy of the migration.
        migration_indices (tuple): Tuple of the indices which designate the
            migrating site and the vacant site to which the cation will
            migrate. If no indices are provided, the user will be prompted.
        functional (tuple): Tuple with the functional details. The first element
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

    # TODO Can currently not be executed from jupyter notebook

    # Let the user define a migration
    migration_dir = define_migration(structure=structure,
                                     site=migration_indices[0],
                                     final_site=migration_indices[1],
                                     write_cif=True)

    # Set up the transition calculation
    transition(directory=migration_dir,
               functional=functional,
               is_metal=is_metal)

    # Create the PyTask that runs the calculation
    if in_custodian:
        vasprun = CustodianTask(directory=os.path.join(migration_dir, "final"))
    else:
        vasprun = VaspTask(directory=os.path.join(migration_dir, "final"))

    # Extract the final cathode from the geometry optimization
    get_cathode = PyTask(
        func="pybat.cli.commands.get.get_cathode",
        kwargs={"directory": os.path.join(migration_dir, "final"),
                "write_cif": True}
    )

    # Add number of nodes to spec, or "none"
    firework_spec = {}
    if number_nodes is None:
        firework_spec.update({"_category": "none"})
    else:
        firework_spec.update({"_category": str(number_nodes) + "nodes"})

    transition_firework = Firework(tasks=[vasprun, get_cathode],
                                   name="Migration Geometry optimization",
                                   spec=firework_spec)

    # Set up the static calculation directory
    static_dir = os.path.join(migration_dir, "static_final")

    final_cathode = os.path.join(migration_dir, "final", "final_cathode.json")

    # Set up the static calculation
    static_fw = PybatStaticFW(
        structure=final_cathode, functional=functional,
        directory=static_dir, write_chgcar=False, in_custodian=in_custodian,
        number_nodes=number_nodes
    )

    struc_name = str(structure.composition.reduced_composition).replace(" ", "")

    return Workflow(
        fireworks=[transition_firework, static_fw],
        name=struc_name + " " + migration_dir.split("/")[-1],
        links_dict={transition_firework: [static_fw]}
    )


def get_wf_neb(directory, nimages=7, functional=("pbe", {}), is_metal=False,
               in_custodian=False, number_nodes=None):
    """
    Set up a workflow that calculates the kinetic barrier between two geometries.

    # TODO
    TEMPORARY? Should NEB be integrated in other workflows? If so, should we still
    have a separate NEB workflow?

    Args:
        directory (str): Directory in which the NEB calculation should be performed.
        nimages (int): Number of images to use for the NEB calculation.
        functional (tuple): Tuple with the functional details. The first element
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
            it is picked up by the right Fireworker. Defaults to the number of images.

    """
    directory = os.path.abspath(directory)

    # If no number of nodes is specified, take the number of images
    if number_nodes is None:
        number_nodes = nimages

    # Create the Firework that sets up and runs the NEB
    neb_firework = NebFirework(
        directory=directory,
        nimages=nimages,
        functional=functional,
        is_metal=is_metal,
        in_custodian=in_custodian,
        number_nodes=number_nodes
    )

    # Add number of nodes to spec, or "none"
    firework_spec = {}
    if number_nodes is None:
        firework_spec.update({"_category": "none"})
    else:
        firework_spec.update({"_category": str(number_nodes) + "nodes"})

    cathode = Cathode.from_file(
        os.path.join(directory, "final", "initial_cathode.json")
    )
    dir_name = os.path.abspath(directory).split("/")[-1]
    workflow_name = str(cathode.composition).replace(" ", "") + " " + dir_name

    return Workflow(
        fireworks=[neb_firework, ],
        name=workflow_name
    )


def get_wf_dimer(structure, directory, dimer_indices, distance,
                 functional=("pbe", {}), is_metal=False, in_custodian=False,
                 number_nodes=None):
    """
    Set up a workflow that calculates the thermodynamics for a dimer
    formation in the current directory.

    Can later be expanded to also include kinetic barrier calculation.

    Args:
        structure (pybat.core.LiRichCathode): LiRichCathode for which to perform a
            dimer workflow. Should be a LiRichCathode, as only for this class the dimer
            scripts are defined.
        directory (str): Path to the directory in which the dimer calculation should
            be run.
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
    setup_dimer = PyTask(
        func="pybat.cli.commands.define.dimer",
        kwargs={"structure": structure,
                "directory": directory,
                "dimer_indices": dimer_indices,
                "distance": distance,
                "write_cif": True}
    )

    # Set up the FireTask that sets up the transition calculation
    setup_transition = PyTask(
        func="pybat.cli.commands.setup.transition",
        kwargs={"directory": directory,
                "functional": functional,
                "is_metal": is_metal}
    )

    # Create the PyTask that runs the calculation
    if in_custodian:
        vasprun = CustodianTask(directory=os.path.join(directory, "final"))
    else:
        vasprun = VaspTask(directory=os.path.join(directory, "final"))

    # Extract the final cathode from the geometry optimization
    get_cathode = PyTask(
        func="pybat.cli.commands.get.get_cathode",
        kwargs={"directory": os.path.join(directory, "final"),
                "write_cif": True}
    )

    # Add number of nodes to spec, or "none"
    firework_spec = {}
    if number_nodes is None:
        firework_spec.update({"_category": "none"})
    else:
        firework_spec.update({"_category": str(number_nodes) + "nodes"})

    dimer_firework = Firework(tasks=[setup_dimer, setup_transition, vasprun, get_cathode],
                              name="Dimer Geometry optimization",
                              spec=firework_spec)

    # Set up the static calculation directory
    static_dir = os.path.join(directory, "static_final")
    final_cathode = os.path.join(directory, "final", "final_cathode.json")

    # Set up the static calculation
    static_fw = PybatStaticFW(
        structure=final_cathode, functional=functional,
        directory=static_dir, write_chgcar=False, in_custodian=in_custodian,
        number_nodes=number_nodes
    )

    wf_name = str(structure.composition.reduced_composition).replace(" ", "")
    wf_name += "\n dimer" + str(dimer_indices)

    return Workflow(
        fireworks=[dimer_firework, static_fw],
        name=wf_name,
        links_dict={dimer_firework: [static_fw]}
    )


def get_wfs_noneq_dimers(structure, distance, functional=("pbe", {}),
                         is_metal=False, in_custodian=False, number_nodes=None):
    """
    Run dimer calculations for all the nonequivalent dimers in a structure.

    Args:
        structure (pybat.core.LiRichCathode): LiRichCathode for which to perform a
            non-equivalent dimer workflow. Should be a LiRichCathode, as only for this
            class the dimer scripts are defined.
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

    dimer_lists = structure.list_noneq_dimers()
    workflows = []

    for dimer_list in dimer_lists:

        # Find the dimer closest to the center of the lattice. Just for
        # visualization purposes.
        central_dimer = [(), 1e10]

        for dimer in dimer_list:

            dimer_center = Dimer(structure, dimer).center
            lattice_center = np.sum(structure.lattice.matrix, 0) / 3

            dist_to_center = np.linalg.norm(dimer_center - lattice_center)

            if dist_to_center < central_dimer[1]:
                central_dimer = [dimer, dist_to_center]

        workflows.append(
            get_wf_dimer(structure=structure,
                         dimer_indices=central_dimer[0],
                         distance=distance,
                         functional=functional,
                         is_metal=is_metal,
                         in_custodian=in_custodian,
                         number_nodes=number_nodes)
        )

    return workflows


# TODO Add to both multiple dimer workflows
# # Create the dimer directory
# dimer_dir = os.path.join(
#     os.getcwd(), "dimer_" + "_".join([str(el) for el in dimer_indices])
# )
# try:
#     os.mkdir(dimer_dir)
# except FileExistsError:
#     warnings.warn("Warning: " + dimer_dir + " already exists, "
#                                             "overwriting...")

def get_wfs_site_dimers(structure, site_index, distance,
                        functional=("pbe", {}), is_metal=False,
                        in_custodian=False, number_nodes=None):
    """
    Run dimer calculations for all the dimers around a site.

    Args:
        structure (pybat.core.LiRichCathode): Structure of the cathode material. Should
            be a LiRichCathode, as only for this class the dimer scripts are defined.
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

    dimer_list = structure.find_noneq_dimers(int(site_index))
    workflows = []

    for dimer in dimer_list:
        workflows.append(
            get_wf_dimer(structure=structure,
                         dimer_indices=dimer,
                         distance=distance,
                         functional=functional,
                         is_metal=is_metal,
                         in_custodian=in_custodian,
                         number_nodes=number_nodes)
        )

    return workflows


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
