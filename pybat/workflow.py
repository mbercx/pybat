# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import subprocess
import shlex

from ruamel.yaml import YAML
from pymongo.errors import ServerSelectionTimeoutError
from pybat.core import LiRichCathode
from pybat.cli.commands.define import define_dimer, define_migration
from pybat.cli.commands.setup import transition
from custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler, \
    UnconvergedErrorHandler
from custodian.vasp.jobs import VaspJob
from fireworks import Firework, LaunchPad, PyTask, Workflow

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
                password=CONFIG["SERVER"].get("password", default="")
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


def run_vasp(directory):
    """
    Method that simply runs VASP in the directory that is specified. Mainly
    used to set up a PyTask that uses outputs from PyTasks in previous
    FireWorks to run VASP in the appropriate directory.

    Args:
        directory (str): Absolute path to the directory in which VASP should be
            run.
    """

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

def relax_workflow(structure_file, directory="", is_metal=False,
                   hse_calculation=False, in_custodian=False):
    """

    Args:
        structure_file:
        is_metal:
        hse_calculation:
        in_custodian:

    Returns:

    """

    # Set up the directory in which to perform the calculation
    current_dir = os.getcwd()

    # If no directory was provided
    if directory == "":

        if hse_calculation:
            directory = os.path.join(current_dir, "hse_relax")
        else:
            directory = os.path.join(current_dir, "dftu_relax")

    # Create the PyTask that sets up the calculation
    setup_relax = PyTask(
        func="pybat.cli.commands.setup.setup_relax",
        kwargs={"structure_file": structure_file,
                "calculation_dir": directory,
                "is_metal":is_metal,
                "hse_calculation": hse_calculation}
    )

    # Create the PyTask that runs the calculation
    if in_custodian:
        run_vasp = PyTask(
            func="pybat.workflow.run_custodian",
            kwargs={"directory": directory}
        )
    else:
        run_vasp = PyTask(
            func="pybat.workflow.run_vasp",
            kwargs={"directory": directory}
        )

    # Combine the two FireTasks into one FireWork
    relax_firework = Firework(tasks=[setup_relax, run_vasp],
                              name="Geometry optimization",
                              spec={"_launch_dir": current_dir})

    # Create the workflow
    workflow = Workflow(fireworks=[relax_firework,],
                        name=structure_file + " Geometry optimization")

    LAUNCHPAD.add_wf(workflow)


def dimer_workflow(structure_file, dimer_indices=(0, 0), distance=0,
                   is_metal=False, hse_calculation=False, in_custodian=False):
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
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV. Defaults to False.
        hse_calculation (bool): Flag that indicates that the hybrid functional
            HSE06 should be used to calculate the exchange-correlation
            energy. Defaults to False.
        in_custodian (bool): Flag that indicates that the calculations
            should be run within a Custodian. Defaults to False.
    """

    # Let the user define a dimer
    dimer_dir = define_dimer(structure_file=structure_file,
                             dimer_indices=dimer_indices,
                             distance=distance,
                             write_cif=True)

    # Set up the transition calculation
    transition(directory=dimer_dir,
               is_metal=is_metal,
               is_migration=False,
               hse_calculation=hse_calculation)

    # Set up the FireTask for the custodian run, if requested
    if in_custodian:
        run_relax = PyTask(
            func="pybat.workflow.run_custodian",
            kwargs={"directory": os.path.join(dimer_dir, "final")}
        )
    else:
        run_relax = PyTask(
            func="pybat.workflow.run_vasp",
            kwargs={"directory": os.path.join(dimer_dir, "final")}
        )

    relax_firework = Firework(tasks=[run_relax],
                              name="Dimer Geometry optimization",
                              spec={"_launch_dir": dimer_dir,
                                    "_category": "2nodes"})

    workflow = Workflow(fireworks=[relax_firework],
                        name=structure_file + dimer_dir.split("/")[-1])

    LAUNCHPAD.add_wf(workflow)


def migration_workflow(structure_file, migration_indices=(0, 0),
                       is_metal=False, hse_calculation=False,
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
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV. Defaults to False.
        hse_calculation (bool): Flag that indicates that the hybrid functional
            HSE06 should be used to calculate the exchange-correlation
            energy. Defaults to False.
        in_custodian (bool): Flag that indicates that the calculations
            should be run within a Custodian. Defaults to False.
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
               is_metal=is_metal,
               is_migration=False,
               hse_calculation=hse_calculation)

    # Set up the FireTask for the custodian run, if requested
    if in_custodian:
        run_relax = PyTask(
            func="pybat.workflow.run_custodian",
            kwargs={"directory": os.path.join(migration_dir, "final")}
        )
    else:
        run_relax = PyTask(
            func="pybat.workflow.run_vasp",
            kwargs={"directory": os.path.join(migration_dir, "final")}
        )

    relax_firework = Firework(tasks=[run_relax],
                              name="Migration Geometry optimization",
                              spec={"_launch_dir": migration_dir,
                                    "_category": "2nodes"})

    workflow = Workflow(fireworks=[relax_firework],
                        name=structure_file + migration_dir.split("/")[-1])

    LAUNCHPAD.add_wf(workflow)


def site_dimers_workflow(structure_file, site_index, distance, is_metal=False,
                         hse_calculation=False, in_custodian=False):
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
        is_metal (bool): Flag that indicates the material being studied is a
            metal, which changes the smearing from Gaussian to second order
            Methfessel-Paxton of 0.2 eV. Defaults to False.
        hse_calculation (bool): Flag that indicates that the hybrid functional
            HSE06 should be used to calculate the exchange-correlation
            energy. Defaults to False.
        in_custodian (bool): Flag that indicates that the calculations
            should be run within a Custodian. Defaults to False.
    """

    lirich = LiRichCathode.from_file(structure_file)
    dimer_list = lirich.find_oxygen_dimers(int(site_index))

    for dimer in dimer_list:
        dimer_workflow(structure_file=structure_file,
                      dimer_indices=dimer,
                      distance=distance,
                      is_metal=is_metal,
                      hse_calculation=hse_calculation,
                      in_custodian=in_custodian)
