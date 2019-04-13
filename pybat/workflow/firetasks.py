# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import subprocess

import numpy as np
from custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler, UnconvergedErrorHandler
from custodian.vasp.jobs import VaspJob
from fireworks import FiretaskBase
from fireworks import Firework, FWAction, ScriptTask
from pymatgen import Structure

from pybat.core import Cathode

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2019, Marnik Bercx, University of Antwerp"
__version__ = "pre-alpha"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Mar 2019"


class VaspTask(FiretaskBase):
    """
    Firetask that represents a VASP calculation run.

    Required params:
        directory (str): Directory in which the VASP calculation should be run.

    """
    required_params = ["directory"]
    _fw_name = "{{pybat.workflow.firetasks.VaspTask}}"

    def run_task(self, fw_spec):
        os.chdir(self["directory"])
        subprocess.run(fw_spec["_fw_env"]["vasp_cmd"], shell=True)


class CustodianTask(FiretaskBase):
    """
    Firetask that represents a calculation run inside a Custodian.

    Required params:
        directory (str): Directory in which the VASP calculation should be run.

    """
    required_params = ["directory"]
    _fw_name = "{{pybat.workflow.firetasks.CustodianTask}}"

    def run_task(self, fw_spec):
        directory = os.path.abspath(self["directory"])
        os.chdir(directory)

        output = os.path.join(directory, "out")
        # TODO Make the output file more general
        vasp_cmd = fw_spec["_fw_env"]["vasp_cmd"].split(" ")

        handlers = [VaspErrorHandler(output_filename=output),
                    UnconvergedErrorHandler(output_filename=output)]

        jobs = [VaspJob(vasp_cmd=vasp_cmd,
                        output_file=output,
                        stderr_file=output,
                        auto_npar=False)]

        c = Custodian(handlers, jobs, max_errors=10)
        c.run()


class PassCathodeTask(FiretaskBase):
    """
    Pass a final cathode structure to children Fireworks in the Workflow through the
    fw_spec and write the final cathode to the directory specified.

    """

    required_params = ["directory"]
    optional_params = ["ignore_magmom"]

    def run_task(self, fw_spec):
        directory = os.path.abspath(self["directory"])
        cathode = Cathode.from_file(os.path.join(directory,
                                                 "initial_cathode.json"))
        cathode.update_sites(directory, ignore_magmom=self.get("ignore_magmom", False))

        cathode.to("json", os.path.join(directory, "final_cathode.json"))

        FWAction()

        return FWAction(update_spec={"structure": cathode.as_dict()})


class PulayTask(FiretaskBase):
    """
    Check if the lattice vectors of a structure have changed significantly during
    the geometry optimization, which could indicate that there where Pulay stresses
    present. If so, start a new geometry optimization with the final structure.

    Required params:
        directory (str): Directory in which the geometry optimization calculation
            was run.

    Optional params:
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
            if fw_action:
                return FWAction.from_dict(fw_action)
            else:
                return FWAction()

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
            firework_spec = {}
            if number_nodes is None or number_nodes == 0:
                firework_spec.update({"_category": "none"})
            else:
                firework_spec.update({"_category": str(number_nodes) + "nodes"})

            # Combine the two FireTasks into one FireWork
            optimize_fw = Firework(tasks=[copy_contcar, vasprun, pulay_task],
                                   name="Pulay Step",
                                   spec=firework_spec)

            return FWAction(additions=optimize_fw)


class ConfigurationTask(FiretaskBase):
    """
    Construct a set of configurations based on the specified parameters.

    Required params:


    """

    required_params = ["structure", "directory", "substitution_sites", "element_list",
                       "sizes"]
    optional_params = ["concentration_restrictions", "max_configurations",
                       "include_existing", "configuration_list"]
    _fw_name = "{{pybat.workflow.firetasks.ConfigurationTask}}"

    def run_task(self, fw_spec):

        # If requested, check for existing configurations in the directory tree
        current_conf_dict = find_configuration_dict(self["directory"])

        # Adjust max_configurations based on existing configurations
        if self.get("max_configurations", None):
            if self.get("include_existing", True):
                max_conf_to_generate = self.get("max_configurations")
            else:
                max_conf_to_generate = self.get("max_configurations") \
                                       + len(current_conf_dict)
        else:
            max_conf_to_generate = None

        if not self.get("configuration_list", None):
            configurations = self["structure"].get_cation_configurations(
                substitution_sites=self["substitution_sites"],
                cation_list=self["element_list"],
                sizes=self["sizes"],
                concentration_restrictions=self.get("configuration_restrictions", None),
                max_configurations=max_conf_to_generate
            )
        else:
            configurations = self.get("configuration_list", None)

        configuration_dict = {}
        conf_number = 0
        for configuration in configurations:

            conf_hash = configuration.__hash__()

            # If the configuration is not found in the directory tree
            if conf_hash not in current_conf_dict.keys():

                # Make sure the configuration directory number is new
                while "conf_" + str(conf_number) in [
                    e for l in [v["directory"].split("/") for v in
                                current_conf_dict.values()] for e in l if "conf" in e]:
                    conf_number += 1

                configuration_dir = os.path.join(self["directory"], "conf_" +
                                                 str(conf_number),
                                                 "prim")

                # Create the configuration directory and write the structure to a file
                os.makedirs(configuration_dir)
                configuration.to("json", os.path.join(configuration_dir,
                                                      "configuration.json"))
                configuration_dict[str(conf_hash)] = {
                    "structure": configuration.as_dict(),
                    "directory": configuration_dir
                }
                conf_number += 1

            elif self.get("include_existing", True):

                configuration_dict[str(conf_hash)] = {
                    "structure": configuration.as_dict(),
                    "directory": os.path.join(
                        self["directory"],
                        current_conf_dict[conf_hash]["directory"]
                    )
                }

            # Break out of the loop if we have enough configurations
            if len(configuration_dict) == self.get("max_configurations", None):
                break

        return FWAction(update_spec={"configuration_dict": configuration_dict})


class EnergyConfTask(FiretaskBase):
    """
    Add a list of FireWorks to the workflow that first optimize the geometry and
    then perform an static calculation for all the configurations in the
    configuration_dict in the fw_spec.

    """

    required_params = ["functional"]
    optional_params = ["in_custodian", "number_nodes"]
    _fw_name = "{{pybat.workflow.firetasks.EnergyConfTask}}"

    def run_task(self, fw_spec):

        try:
            configuration_dict = fw_spec["configuration_dict"]
        except KeyError:
            raise KeyError("Can not find the configuration_dict in the fw_spec.")

        functional_dir = self["functional"][0]
        if self["functional"][0] == "pbeu":
            functional_dir += "_" + "".join(
                k + str(self["functional"][1]["LDAUU"][k])
                for k in self["functional"][1]["LDAUU"].keys()
            )

        firework_list = []

        for configuration in configuration_dict.values():
            optimize_dir = os.path.join(
                configuration["directory"], functional_dir + "_optimize"
            )
            static_dir = os.path.join(
                configuration["directory"], functional_dir + "_static"
            )

            # This import needs to happen here because the Fireworks depend on the
            # firetasks in this module. # TODO fix this?
            from pybat.workflow.fireworks import PybatStaticFW, PybatOptimizeFW

            static_fw = PybatStaticFW(
                structure=os.path.join(optimize_dir, "final_cathode.json"),
                functional=self["functional"],
                directory=static_dir,
                write_chgcar=False,
                in_custodian=self.get("in_custodian", False),
                number_nodes=self.get("number_nodes", None)
            )

            fw_action = FWAction(additions=static_fw)

            optimize_fw = PybatOptimizeFW(
                structure=configuration["structure"],
                functional=self["functional"],
                directory=optimize_dir,
                in_custodian=self.get("in_custodian", False),
                number_nodes=self.get("number_nodes", None),
                fw_action=fw_action
            )

            firework_list.append(optimize_fw)

        return FWAction(additions=firework_list)


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

# region * Utility scripts

def find_all(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            result.append(os.path.join(root, name))
    return result


def find_configuration_dict(path):
    path = os.path.abspath(path)
    return {Cathode.from_file(file).__hash__(): {
        "structure": Cathode.from_file(file).as_dict(),
        "directory": file.replace(path, "").replace("configuration.json", "").strip("/")
    } for file in find_all("configuration.json", path)}

# endregion
