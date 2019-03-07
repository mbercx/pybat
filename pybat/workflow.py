# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import subprocess
import ast
import pdb

import numpy as np

from pymatgen import Structure

from pybat.core import Cathode, LiRichCathode, Dimer
from pybat.cli.commands.define import define_dimer, define_migration
from pybat.cli.commands.setup import transition

from ruamel.yaml import YAML
from pymongo.errors import ServerSelectionTimeoutError
from custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler, \
    UnconvergedErrorHandler
from custodian.vasp.jobs import VaspJob
from fireworks import FiretaskBase, Firework, LaunchPad, PyTask, Workflow, FWAction, \
    ScriptTask

"""
Workflow setup for the pybat package.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Jul 2018"

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

# TODO Add good description
PULAY_TOLERANCE = 1e-2


# TODO Extend configuration and make the whole configuration setup more user friendly
# Currently the user is not guided to the workflow setup when attempting to use
# pybat workflows, this should change and be tested. Moreover, careful additions should
# be made to make sure all user-specific configuration elements are easily configured
# and implemented in the code.

# TODO Create methods that return FireWorks, so the workflow methods can be modularized
# At this point, it's becoming clear that the workflows are getting more and more
# extensive, and are simply becoming less easy to grasp. It might be useful to create
# methods that set up the FireWorks (e.g. for a relaxation, SCF calculations), and
# then call upon these methods in the workflow methods.

# TODO Fix the CustodianTask

# TODO Add UnitTests!
# It's really getting time to do this. Think about what unit tests you need and make a
# test suite.

# region * Region 1 - Firetasks

class VaspTask(FiretaskBase):
    """
    Firetask that represents a VASP calculation run.

    Required parameters:
        directory (str): Directory in which the VASP calculation should be run.

    """
    required_params = ["directory"]
    _fw_name = "{{pybat.workflow.VaspTask}}"

    def run_task(self, fw_spec):
        os.chdir(self["directory"])
        subprocess.run(fw_spec["_fw_env"]["vasp_command"], shell=True)


class CustodianTask(FiretaskBase):
    """
    Firetask that represents a calculation run inside a Custodian.

    Required parameters:
        directory (str): Directory in which the VASP calculation should be run.

    """
    required_params = ["directory"]
    _fw_name = "{{pybat.workflow.CustodianTask}}"

    def run_task(self, fw_spec):
        directory = os.path.abspath(self["directory"])
        os.chdir(directory)

        output = os.path.join(directory, "out")
        # TODO Make the output file more general
        vasp_cmd = fw_spec["_fw_env"]["vasp_command"]

        handlers = [VaspErrorHandler(output_filename=output),
                    UnconvergedErrorHandler(output_filename=output)]

        jobs = [VaspJob(vasp_cmd=vasp_cmd,
                        output_file=output,
                        stderr_file=output,
                        auto_npar=False)]

        c = Custodian(handlers, jobs, max_errors=10)
        c.run()


class PulayTask(FiretaskBase):
    """
    Check if the lattice vectors of a structure have changed significantly during
    the geometry optimization, which could indicate that there where Pulay stresses
    present. If so, start a new geometry optimization with the final structure.

    Required parameters:
        directory (str): Directory in which the geometry optimization calculation
            was run.

    Optional parameters:
        in_custodian (bool): Flag that indicates whether the calculation should be
            run inside a Custodian.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker.
        tolerance (float): Tolerance that indicates the maximum change in norm for the
            matrix defined by the cartesian coordinates of the lattice vectors.
            If the norm changes more than the tolerance, another geometry optimization
            is performed starting from the final geometry.

    """
    required_params = ["directory"]
    option_params = ["in_custodian", "number_nodes", "tolerance", "fw_action"]
    _fw_name = "{{pybat.workflow.PulayTask}}"

    def run_task(self, fw_spec):
        """

        Args:
            fw_spec:

        Returns:
            FWAction

        """
        # Extract the parameters into variables; this makes for cleaner code IMO
        directory = self["directory"]
        in_custodian = self.get("in_custodian", False)
        number_nodes = self.get("number_nodes", None)
        tolerance = self.get("tolerance", PULAY_TOLERANCE)
        fw_action = self.get('fw_action', {})

        pdb.set_trace()

        # Check if the lattice vectors have changed significantly
        initial_structure = Structure.from_file(
            os.path.join(directory, "POSCAR")
        )
        final_structure = Structure.from_file(
            os.path.join(directory, "CONTCAR")
        )

        sum_differences = np.linalg.norm(
            initial_structure.lattice.matrix - final_structure.lattice.matrix
        )

        # If the difference is small, return an empty FWAction
        if sum_differences < tolerance:
            return FWAction.from_dict(fw_action)

        # Else, set up another geometry optimization
        else:
            print("Lattice vectors have changed significantly during geometry "
                  "optimization. Performing another full geometry optimization to "
                  "make sure there were no Pulay stresses present.\n\n")

            # Create the ScriptTask that copies the CONTCAR to the POSCAR
            copy_contcar = ScriptTask.from_str(
                "cp " + os.path.join(directory, "CONTCAR") +
                " " + os.path.join(directory, "POSCAR")
            )

            # Create the PyTask that runs the calculation
            if in_custodian:
                vasprun = CustodianTask(directory=directory)
            else:
                vasprun = VaspTask(directory=directory)

            # Create the PyTask that check the Pulay stresses again
            pulay_task = PulayTask(
                directory=directory, in_custodian=in_custodian,
                number_nodes=number_nodes, tolerance=tolerance,
                fw_action=fw_action
            )

            # Add number of nodes to spec, or "none"
            firework_spec = {"_launch_dir": os.getcwd()}
            if number_nodes is None:
                firework_spec.update({"_category": "none"})
            else:
                firework_spec.update({"_category": str(number_nodes) + "nodes"})

            # Combine the two FireTasks into one FireWork
            relax_firework = Firework(tasks=[copy_contcar, vasprun, pulay_task],
                                      name="Pulay Step",
                                      spec=firework_spec)

            return FWAction(additions=relax_firework)


# endregion

# region * Region 2 - Fireworks


def create_scf_fw(structure_file, functional, directory, write_chgcar=False,
                  in_custodian=False,
                  number_nodes=None):
    """
    Create a FireWork for performing an SCF calculation.

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
        Firework: A firework that represents an SCF calculation.

    """
    # Create the PyTask that sets up the calculation
    setup_scf = PyTask(
        func="pybat.cli.commands.setup.scf",
        kwargs={"structure_file": structure_file,
                "functional": functional,
                "calculation_dir": directory,
                "write_chgcar": write_chgcar}
    )

    # Create the PyTask that runs the calculation
    if in_custodian:
        vasprun = CustodianTask(directory=directory)
    else:
        vasprun = VaspTask(directory=directory)

    # Add number of nodes to spec, or "none"
    firework_spec = {"_launch_dir": os.getcwd()}
    if number_nodes is None:
        firework_spec.update({"_category": "none"})
    else:
        firework_spec.update({"_category": str(number_nodes) + "nodes"})

    # Combine the two FireTasks into one FireWork
    scf_firework = Firework(tasks=[setup_scf, vasprun],
                            name="SCF calculation",
                            spec=firework_spec)

    return scf_firework


class RelaxFirework(Firework):

    def __init__(self, structure_file, functional, directory, is_metal=False,
                 in_custodian=False, number_nodes=None, fw_action=None):

        # Create the PyTask that sets up the calculation
        setup_relax = PyTask(
            func="pybat.cli.commands.setup.relax",
            kwargs={"structure_file": structure_file,
                    "functional": functional,
                    "calculation_dir": directory,
                    "is_metal": is_metal}
        )

        # Create the PyTask that runs the calculation
        if in_custodian:
            vasprun = CustodianTask(directory=directory)
        else:
            vasprun = VaspTask(directory=directory)

        # Extract the final cathode from the geometry optimization
        get_cathode = PyTask(
            func="pybat.cli.commands.get.get_cathode",
            kwargs={"directory": os.path.join(directory),
                    "write_cif": True}
        )

        # Create the PyTask that check the Pulay stresses
        pulay_task = PulayTask(directory=directory,
                               in_custodian=in_custodian,
                               number_nodes=number_nodes,
                               tol=PULAY_TOLERANCE,
                               fw_action=fw_action)

        # Only add number of nodes to spec if specified
        firework_spec = {"_launch_dir": os.getcwd()}
        if number_nodes is None:
            firework_spec.update({"_category": "none"})
        else:
            firework_spec.update({"_category": str(number_nodes) + "nodes"})

        # Combine the FireTasks into one FireWork
        super(RelaxFirework, self).__init__(
            tasks=[setup_relax, vasprun, get_cathode, pulay_task],
            name="Geometry optimization", spec=firework_spec
        )


def create_relax_fw(structure_file, functional, directory, is_metal=False,
                    in_custodian=False, number_nodes=None, fw_action=None):
    """
    Create a FireWork for performing an Geometry Optimization.

    Args:
        structure_file (str): Path to the geometry file of the structure.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        directory (str): Directory in which the SCF calculation should be performed.
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV.
        in_custodian (bool): Flag that indicates whether the calculation should be
            run inside a Custodian.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker.

    Returns:
        Firework: A firework that represents an SCF calculation.

    """
    # Create the PyTask that sets up the calculation
    setup_relax = PyTask(
        func="pybat.cli.commands.setup.relax",
        kwargs={"structure_file": structure_file,
                "functional": functional,
                "calculation_dir": directory,
                "is_metal": is_metal}
    )

    # Create the PyTask that runs the calculation
    if in_custodian:
        vasprun = CustodianTask(directory=directory)
    else:
        vasprun = VaspTask(directory=directory)

    # Extract the final cathode from the geometry optimization
    get_cathode = PyTask(
        func="pybat.cli.commands.get.get_cathode",
        kwargs={"directory": os.path.join(directory),
                "write_cif": True}
    )

    # Create the PyTask that check the Pulay stresses
    pulay_task = PulayTask(directory=directory,
                           in_custodian=in_custodian,
                           number_nodes=number_nodes,
                           tol=PULAY_TOLERANCE,
                           fw_action=fw_action)

    # Only add number of nodes to spec if specified
    firework_spec = {"_launch_dir": os.getcwd()}
    if number_nodes is None:
        firework_spec.update({"_category": "none"})
    else:
        firework_spec.update({"_category": str(number_nodes) + "nodes"})

    # Combine the FireTasks into one FireWork
    relax_firework = Firework(tasks=[setup_relax, vasprun, get_cathode, pulay_task],
                              name="Geometry optimization",
                              spec=firework_spec)

    return relax_firework


def create_neb_fw(directory, nimages, functional, is_metal=False, is_migration=False,
                  in_custodian=False, number_nodes=None):
    """
    Create a FireWork for performing an NEB calculation.

    Args:
        directory (str): Directory in which the NEB calculation should be performed.
        nimages (int): Number of images to use for the NEB calculation.
        functional (tuple): Tuple with the functional choices. The first element
            contains a string that indicates the functional used ("pbe", "hse", ...),
            whereas the second element contains a dictionary that allows the user
            to specify the various functional tags.
        in_custodian (bool): Flag that indicates whether the calculation should be
            run inside a Custodian.
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV.
        is_migration (bool): Flag that indicates that the transition is a migration
            of an atom in the structure.
        number_nodes (int): Number of nodes that should be used for the calculations.
            Is required to add the proper `_category` to the Firework generated, so
            it is picked up by the right Fireworker.

    Returns:
        Firework: A firework that represents an NEB calculation.

    """
    # Create the PyTask that sets up the calculation
    setup_neb = PyTask(
        func="pybat.cli.commands.setup.neb",
        kwargs={"directory": directory,
                "nimages": nimages,
                "functional": functional,
                "is_metal": is_metal,
                "is_migration": is_migration}
    )

    # Create the PyTask that runs the calculation
    if in_custodian:
        vasprun = CustodianTask(directory=directory)
    else:
        vasprun = VaspTask(directory=directory)

    # Add number of nodes to spec, or "none"
    firework_spec = {"_launch_dir": os.getcwd()}
    if number_nodes == None:
        firework_spec.update({"_category": "none"})
    else:
        firework_spec.update({"_category": str(number_nodes) + "nodes"})

    # Combine the two FireTasks into one FireWork
    neb_firework = Firework(tasks=[setup_neb, vasprun],
                            name="NEB calculation",
                            spec=firework_spec)

    return neb_firework


# endregion

# region * Region 3 - Workflows


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
    scf_firework = create_scf_fw(
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
    relax_firework = create_relax_fw(structure_file=structure_file,
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

    relax_firework = Firework(tasks=[setup_transition, vasprun, get_cathode],
                              name="Dimer Geometry optimization",
                              spec=firework_spec)

    # Set up the SCF calculation directory
    scf_dir = os.path.join(dimer_dir, "scf_final")

    final_cathode = os.path.join(dimer_dir, "final", "final_cathode.json")

    # Set up the SCF calculation
    scf_firework = create_scf_fw(
        structure_file=final_cathode, functional=functional,
        directory=scf_dir, write_chgcar=False, in_custodian=in_custodian,
        number_nodes=number_nodes
    )

    workflow = Workflow(fireworks=[relax_firework, scf_firework],
                        name=structure_file + dimer_dir.split("/")[-1],
                        links_dict={relax_firework: [scf_firework]})

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

    relax_firework = Firework(tasks=[vasprun],
                              name="Migration Geometry optimization",
                              spec=firework_spec)

    workflow = Workflow(fireworks=[relax_firework],
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
    neb_firework = create_neb_fw(
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

    # TODO Improve naming scheme of workflow
    workflow = Workflow(fireworks=[neb_firework, ],
                        name=directory.split("/")[-1])

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

            scf_firework = create_scf_fw(
                structure_file=os.path.join(relax_dir, "final_cathode.json"),
                functional=functional,
                directory=scf_dir,
                write_chgcar=False,
                in_custodian=in_custodian,
                number_nodes=number_nodes
            )
            fw_action = FWAction(additions=scf_firework)

            pdb.set_trace()

            firework_list.append(create_relax_fw(
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

            scf_firework = create_scf_fw(
                structure_file=os.path.join(relax_dir, "final_cathode.json"),
                functional=functional,
                directory=scf_dir,
                write_chgcar=False,
                in_custodian=in_custodian,
                number_nodes=number_nodes
            )
            fw_action = FWAction(additions=scf_firework)

            firework_list.append(create_relax_fw(
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


# endregion

# region * Region 4 - Studies
#
# Studies are a collection of Workflows


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


# endregion

# region * Token tasks for testing

class MiddleTask(FiretaskBase):
    required_params = ["message"]
    option_params = ["fw_action"]
    _fw_name = "{{pybat.workflow.MiddleTask}}"

    def run_task(self, fw_spec):

        if self["message"] == "next":

            print("The message was next!")
            print(self["message"])

            fw = Firework([ScriptTask.from_str("echo next"),
                           MiddleTask(message="other",
                                      fw_action=self.get("fw_action", FWAction()))])

            return FWAction(additions=fw)

        else:

            print("The message was Something Else!")
            print(self["message"])
            print()
            print(self.get("fw_action", FWAction()))
            print()

            return FWAction.from_dict(self.get("fw_action", {}))


class FirstFirework(Firework):

    def __init__(self, final_message):
        super(FirstFirework, self).__init__(
            tasks=[ScriptTask.from_str("echo 'Here we go!'"),
                   MiddleTask(message="next",
                              fw_action=FWAction(additions=Firework([ScriptTask.from_str(
                                  "echo '" + final_message + "'")])))])
