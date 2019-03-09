# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import subprocess
import numpy as np

from pymatgen import Structure
from custodian import Custodian
from custodian.vasp.jobs import VaspJob
from custodian.vasp.handlers import VaspErrorHandler, UnconvergedErrorHandler
from fireworks import Firework, FWAction, ScriptTask

from fireworks import FiretaskBase

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2019, Marnik Bercx, University of Antwerp"
__version__ = "pre-alpha"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Mar 2019"


class VaspTask(FiretaskBase):
    """
    Firetask that represents a VASP calculation run.

    Required parameters:
        directory (str): Directory in which the VASP calculation should be run.

    """
    required_params = ["directory"]
    _fw_name = "{{pybat.workflow.firetasks.VaspTask}}"

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
    _fw_name = "{{pybat.workflow.firetasks.CustodianTask}}"

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
    _fw_name = "{{pybat.workflow.firetasks.PulayTask}}"

    # Standard tolerance for deciding to perform another geometry optimization.
    # Basically, PulayTask calculates the 2-norm of the absolute matrix taken from the
    # difference between the initial and final matrices of the lattice vectors of the
    # structure.
    pulay_tolerance = 1e-2

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
        tolerance = self.get("tolerance", PulayTask.pulay_tolerance)
        fw_action = self.get('fw_action', {})

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


# region * Token FireTasks for testing

class MiddleTask(FiretaskBase):
    required_params = ["message"]
    option_params = ["fw_action"]
    _fw_name = "{{pybat.workflow.firetask.MiddleTask}}"

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

# endregion
