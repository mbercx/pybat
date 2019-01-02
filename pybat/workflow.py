# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import subprocess
import shlex

import numpy as np

from pybat.core import LiRichCathode, Dimer
from pybat.cli.commands.define import define_dimer, define_migration
from pybat.cli.commands.setup import transition

from ruamel.yaml import YAML
from pymongo.errors import ServerSelectionTimeoutError
from custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler, \
    UnconvergedErrorHandler
from custodian.vasp.jobs import VaspJob
from fireworks import Firework, LaunchPad, PyTask, Workflow, FWAction, ScriptTask

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
                host=CONFIG["SERVER"].get("host", default=""),
                port=int(CONFIG["SERVER"].get("port", default=0)),
                name=CONFIG["SERVER"].get("name", default=""),
                username=CONFIG["SERVER"].get("username", default=""),
                password=CONFIG["SERVER"].get("password", default=""),
                ssl=CONFIG["SERVER"].get("ssl", default=False),
                authsource=CONFIG["SERVER"].get("authsource", default=None)

            )
        except ServerSelectionTimeoutError:
            raise TimeoutError("Could not connect to server. Please make "
                               "sure the details of the server are correctly "
                               "set up.")

        VASP_RUN_SCRIPT = CONFIG["WORKFLOW"].get("script_path", default="")
        VASP_RUN_COMMAND = "bash " + VASP_RUN_SCRIPT
else:
    raise FileNotFoundError("No configuration file found in user's home "
                            "directory. Please use pybat workflow setup "
                            "in order to set up the configuration for "
                            "the workflows.")


# TODO Create methods that return FireWorks, so the workflow methods can be modularized
# At this point, it's becoming clear that the workflows are getting more and more
# extensive, and are simply becoming less easy to grasp. It might be useful to create
# methods that set up the FireWorks (e.g. for a relaxation, SCF calculations), and
# then call upon these methods in the workflow methods.

# TODO Turn run_vasp/run_custodian into FireTasks

# TODO Fix the custodian issue
# Currently custodian does not terminate the previous job properly. This may be related
# to the fact that the vasp run command is called in a script, and so custodian can
# only terminate the script, not the actual vasp run.

def run_vasp(directory, number_nodes):
    """
    Method that simply runs VASP in the directory that is specified. Mainly
    used to set up a PyTask that uses outputs from PyTasks in previous
    FireWorks to run VASP in the appropriate directory.

    Args:
        directory (str): Absolute path to the directory in which VASP should be
            run.
        number_nodes (str):  #TODO

    """
    # Workaround for making number of nodes work on breniac #TODO

    os.chdir(directory)
    subprocess.call(VASP_RUN_SCRIPT)


def run_custodian(directory):
    """
    Run VASP under supervision of a custodian in a certain directory.

    Args:
        directory (str): Absolute path to the directory in which VASP should be
            run.
    """

    directory = os.path.abspath(directory)
    os.chdir(directory)

    output = os.path.join(directory, "out")
    vasp_cmd = shlex.split(VASP_RUN_COMMAND)

    handlers = [VaspErrorHandler(output_filename=output),
                UnconvergedErrorHandler(output_filename=output)]

    jobs = [VaspJob(vasp_cmd=vasp_cmd,
                    output_file=output,
                    stderr_file=output,
                    auto_npar=False)]

    c = Custodian(handlers, jobs, max_errors=10)
    c.run()


def create_scf_fw(structure_file, functional, directory, write_chgcar, in_custodian,
                  number_nodes):
    """
    Create a FireWork for performing an SCF calculation

    Returns:

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
        vasprun = PyTask(
            func="pybat.workflow.run_custodian",
            kwargs={"directory": directory}
        )
    else:
        vasprun = PyTask(
            func="pybat.workflow.run_vasp",
            kwargs={"directory": directory,
                    "number_nodes": number_nodes}
        )

    # Combine the two FireTasks into one FireWork
    scf_firework = Firework(tasks=[setup_scf, vasprun],
                            name="SCF calculation",
                            spec={"_launch_dir": os.getcwd(),
                                  "_category": number_nodes})

    return scf_firework


def pulay_check(directory, in_custodian, number_nodes, tol=1e-2):
    """
    Check if the lattice vectors of a structure have changed significantly during
    the geometry optimization, which could indicate that there where Pulay stresses
    present. If so, start a new geometry optimization with the final structure.

    Args:
        directory (str): Directory in which the geometry optimization calculation
            was run.
        in_custodian (bool): Flag that indicates wheter the calculation should be
            run inside a Custodian.
        number_nodes (): #TODO
        tol (float): Tolerance that indicates the maximum change in norm for the
            matrix defined by the cartesian coordinates of the lattice vectors.
            If the norm changes more than the tolerance, another geometry optimization
            is performed starting from the final geometry.

    Returns:
        FWAction

    """
    # Check if the lattice vectors have changed significantly

    initial_cathode = LiRichCathode.from_file(
        os.path.join(directory, "POSCAR")
    )
    final_cathode = LiRichCathode.from_file(
        os.path.join(directory, "CONTCAR")
    )

    sum_differences = np.linalg.norm(
        initial_cathode.lattice.matrix - final_cathode.lattice.matrix
    )

    if sum_differences < tol:
        return FWAction()

    else:
        print("Lattice vectors have changed significantly during geometry "
              "optimization. Performing another full geometry optimization to "
              "make sure there were no Pulay stresses present.")

        # Create the ScriptTask that copies the CONTCAR to the POSCAR
        copy_contcar = ScriptTask.from_str(
            "cp " + os.path.join(directory, "CONTCAR") +
            " " + os.path.join(directory, "POSCAR")
        )

        # Create the PyTask that runs the calculation
        if in_custodian:
            vasprun = PyTask(
                func="pybat.workflow.run_custodian",
                kwargs={"directory": directory}
            )
        else:
            vasprun = PyTask(
                func="pybat.workflow.run_vasp",
                kwargs={"directory": directory,
                        "number_nodes": number_nodes}
            )

        # Create the PyTask that check the Pulay stresses again
        pulay_task = PyTask(
            func="pybat.workflow.check_pulay",
            kwargs={"directory": directory,
                    "in_custodian": in_custodian,
                    "number_nodes": number_nodes}
        )

        # Combine the two FireTasks into one FireWork
        relax_firework = Firework(tasks=[copy_contcar, vasprun, pulay_task],
                                  name="Pulay Step",
                                  spec={"_launch_dir": directory,
                                        "_category": number_nodes})

        return FWAction(additions=relax_firework)


def scf_workflow(structure_file, functional=("pbe", {}), directory="",
                 write_chgcar=False, in_custodian=False):
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
        in_custodian (bool): Flag that indicates wheter the calculation should be
            run inside a Custodian.

    Returns:
        None

    """

    # Set up the directory in which to perform the calculation
    current_dir = os.getcwd()

    if functional[0] == "hse":
        number_nodes = "4nodes"
    else:
        number_nodes = "1node"

    # If no directory was provided, set it up according to the functional
    if directory == "":
        directory = os.path.join(current_dir, functional[0] + "_scf")

    # Combine the two FireTasks into one FireWork
    scf_firework = create_scf_fw(
        structure_file=structure_file, functional=functional,
        directory=directory, write_chgcar=write_chgcar,
        in_custodian=in_custodian, number_nodes=number_nodes
    )

    # Set up a clear name for the workflow
    cathode = LiRichCathode.from_file(structure_file)
    workflow_name = str(cathode.composition.reduced_formula).replace(" ", "")
    workflow_name += print(functional)

    # Create the workflow
    workflow = Workflow(fireworks=[scf_firework, ],
                        name=workflow_name)

    LAUNCHPAD.add_wf(workflow)


def relax_workflow(structure_file, functional=("pbe", {}), directory="",
                   is_metal=False, in_custodian=False):
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

    Returns:
        None

    """

    # Set up the directory in which to perform the calculation
    current_dir = os.getcwd()

    if functional[0] == "hse":
        number_nodes = "4nodes"
    else:
        number_nodes = "1node"

    # If no directory was provided, set it up according to the functional
    if directory == "":
        directory = os.path.join(current_dir, functional[0] + "_relax")
    else:
        directory = os.path.join(current_dir, directory)

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
        vasprun = PyTask(
            func="pybat.workflow.run_custodian",
            kwargs={"directory": directory}
        )
    else:
        vasprun = PyTask(
            func="pybat.workflow.run_vasp",
            kwargs={"directory": directory,
                    "number_nodes": number_nodes}
        )

    # Create the PyTask that check the Pulay stresses
    pulay_task = PyTask(
        func="pybat.workflow.check_pulay",
        kwargs={"directory": directory,
                "in_custodian": in_custodian,
                "number_nodes": number_nodes}
    )

    # Combine the FireTasks into one FireWork
    relax_firework = Firework(tasks=[setup_relax, vasprun, pulay_task],
                              name="Geometry optimization",
                              spec={"_launch_dir": current_dir,
                                    "_category": number_nodes})

    # Set up a clear name for the workflow
    cathode = LiRichCathode.from_file(structure_file)
    workflow_name = str(cathode.composition.reduced_formula).replace(" ", "")
    workflow_name += print(functional)

    # Create the workflow
    workflow = Workflow(fireworks=[relax_firework, ],
                        name=workflow_name)

    LAUNCHPAD.add_wf(workflow)


def dimer_workflow(structure_file, dimer_indices=(0, 0), distance=0,
                   functional=("pbe", {}), is_metal=False, in_custodian=False):
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
    """
    # TODO Change naming scheme

    if functional[0] == "hse":
        number_nodes = "4nodes"
    else:
        number_nodes = "1node"

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

    # Set up the FireTask for the custodian run, if requested (lala)
    if in_custodian:
        run_relax = PyTask(
            func="pybat.workflow.run_custodian",
            kwargs={"directory": os.path.join(dimer_dir, "final")}
        )
    else:
        run_relax = PyTask(
            func="pybat.workflow.run_vasp",
            kwargs={"directory": os.path.join(dimer_dir, "final"),
                    "number_nodes": number_nodes}
        )

    # Extract the final cathode from the geometry optimization
    get_cathode = PyTask(
        func="pybat.cli.commands.get.get_cathode",
        kwargs={"directory": os.path.join(dimer_dir, "final"),
                "write_cif": True, }
    )

    relax_firework = Firework(tasks=[setup_transition, run_relax, get_cathode],
                              name="Dimer Geometry optimization",
                              spec={"_launch_dir": dimer_dir,
                                    "_category": "2nodes"})

    # Set up the SCF calculation following the relaxation, in order to get an accurate
    # total energy
    # Create the PyTask that sets up the calculation
    if functional[0] == "pbeu":
        scf_dir = os.path.join(
            dimer_dir,
            "pbeu_" +
            "".join([k + str(functional[1][k]) for k in functional[1].keys()]) +
            "_scf"
        )
    else:
        scf_dir = os.path.join(dimer_dir, functional[0] + "_scf")

    final_cathode = os.path.join(dimer_dir, "final", "final_cathode.json")

    scf_firework = create_scf_fw(
        structure_file=final_cathode, functional=functional,
        directory=scf_dir, write_chgcar=False, in_custodian=in_custodian,
        number_nodes=number_nodes)

    workflow = Workflow(fireworks=[relax_firework, scf_firework],
                        name=structure_file + dimer_dir.split("/")[-1],
                        links_dict={relax_firework: [scf_firework]})

    LAUNCHPAD.add_wf(workflow)


def migration_workflow(structure_file, migration_indices=(0, 0),
                       functional=("pbe", {}), is_metal=False,
                       in_custodian=False):
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
    """

    if functional[0] == "hse":
        number_nodes = "4nodes"
    else:
        number_nodes = "1node"

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

    # Set up the FireTask for the custodian run, if requested
    if in_custodian:
        run_relax = PyTask(
            func="pybat.workflow.run_custodian",
            kwargs={"directory": os.path.join(migration_dir, "final")}
        )
    else:
        run_relax = PyTask(
            func="pybat.workflow.run_vasp",
            kwargs={"directory": os.path.join(migration_dir, "final"),
                    "number_nodes": number_nodes}
        )

    relax_firework = Firework(tasks=[run_relax],
                              name="Migration Geometry optimization",
                              spec={"_launch_dir": migration_dir,
                                    "_category": "2nodes"})

    workflow = Workflow(fireworks=[relax_firework],
                        name=structure_file + migration_dir.split("/")[-1])

    LAUNCHPAD.add_wf(workflow)


def noneq_dimers_workflow(structure_file, distance, functional=("pbe", {}),
                          is_metal=False, in_custodian=False):
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
                       in_custodian=in_custodian)


def site_dimers_workflow(structure_file, site_index, distance,
                         functional=("pbe", {}), is_metal=False,
                         in_custodian=False):
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
    """

    lirich = LiRichCathode.from_file(structure_file)
    dimer_list = lirich.find_oxygen_dimers(int(site_index))

    for dimer in dimer_list:
        dimer_workflow(structure_file=structure_file,
                       dimer_indices=dimer,
                       distance=distance,
                       functional=functional,
                       is_metal=is_metal,
                       in_custodian=in_custodian)
