# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import subprocess
import shlex

from custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler, \
    UnconvergedErrorHandler
from custodian.vasp.jobs import VaspJob

from fireworks import Firework, LaunchPad, PyTask, FWorker, \
    Workflow
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

"""
Workflow setup for the pybat package.

"""

# Path to the
VASP_RUN_SCRIPT = "/user/antwerpen/202/vsc20248/local/scripts/job_workflow.sh"
VASP_RUN_COMMAND = "mpirun -genv LD_BIND_NOW=1 vasp_std"

# Set up the Launchpad for the workflows
LAUNCHPAD = LaunchPad(host="ds135179.mlab.com", port=35179, name="quotas",
                      username="mbercx", password="quotastests")

def run_vasp(directory):
    """
    Method that simply runs VASP in the directory that is specified. Mainly
    used to set up a PyTask that uses outputs from PyTasks in previous
    FireWorks to run VASP in the appropriate directory.

    Args:
        directory: Absolute path to the directory in which VASP should be run.

    Returns:
        None

    """

    os.chdir(directory)
    subprocess.call(VASP_RUN_SCRIPT)


def run_custodian(directory):
    """
    Run VASP under supervision of a custodian in a certain directory.

    Args:
        directory:

    Returns:

    """

    directory = os.path.abspath(directory)
    os.chdir(directory)

    output = os.path.join(directory, "out")
    vasp_cmd = shlex.split(VASP_RUN_COMMAND)

    handlers = [VaspErrorHandler(output_filename=output),
                UnconvergedErrorHandler(output_filename=output)]

    jobs = [VaspJob(vasp_cmd=vasp_cmd,
                    output_file=output,
                    stderr_file=output)]

    c = Custodian(handlers, jobs, max_errors=10)
    c.run()


def transition_workflow(structure_file, dimer_indices=(0, 0), distance=0):
    """
    Set up a workflow that calculates the geometry of a transition, within a
    Custodian.

    Returns:

    ""

    # Define the transition, based on the original structure
    from pybat.cli.commands.define import define_dimer

    define_dimer(structure_file=structure_file,
                 dimer_indices=dimer_indices,
                 distance=distance,
                 write_cif=True)

    # Set up the calculation



    # Run the calculation



