# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os

from fireworks import Firework, FWAction, ScriptTask, PyTask

from pybat.workflow.firetasks import VaspTask, CustodianTask, PulayTask, MiddleTask

"""
Package that contains all the fireworks used by pybat to construct Workflows.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2019, Marnik Bercx, University of Antwerp"
__version__ = "pre-alpha"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Mar 2019"


class PybatStaticFW(Firework):

    def __init__(self, structure, functional, directory, write_chgcar=False,
                 in_custodian=False, number_nodes=None):
        """
        Create a FireWork for performing a static calculation.

        Args:
            structure: pymatgen.Structure OR path to structure file for which to run
                the static calculation.
            functional (tuple): Tuple with the functional choices. The first element
                contains a string that indicates the functional used ("pbe", "hse", ...),
                whereas the second element contains a dictionary that allows the user
                to specify the various functional tags.
            directory (str): Directory in which the static calculation should be
                performed.
            write_chgcar (bool): Flag that indicates whether the CHGCAR file should
                be written.
            in_custodian (bool): Flag that indicates whether the calculation should be
                run inside a Custodian.
            number_nodes (int): Number of nodes that should be used for the calculations.
                Is required to add the proper `_category` to the Firework generated, so
                it is picked up by the right Fireworker.

        Returns:
            Firework: A firework that represents a static calculation.

        """

        # Create the PyTask that sets up the calculation
        setup_static = PyTask(
            func="pybat.cli.commands.setup.static",
            kwargs={"structure": structure,
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
        firework_spec = {}
        if number_nodes is None or number_nodes == 0:
            firework_spec.update({"_category": "none"})
        else:
            firework_spec.update({"_category": str(number_nodes) + "nodes"})

        # Combine the two FireTasks into one FireWork
        super(PybatStaticFW, self).__init__(
            tasks=[setup_static, vasprun], name="Static calculation", spec=firework_spec
        )


class PybatOptimizeFW(Firework):

    def __init__(self, structure, functional, directory, is_metal=False,
                 in_custodian=False, number_nodes=None, fw_action=None):
        """
        Initialize a Firework for a geometry optimization.

        Args:
            structure: pymatgen.Structure OR path to structure file for which to run
                the geometry optimization.
            functional (tuple): Tuple with the functional choices. The first element
                contains a string that indicates the functional used ("pbe", "hse", ...),
                whereas the second element contains a dictionary that allows the user
                to specify the various functional tags.
            directory (str): Directory in which the geometry optimization should be
                performed.
            is_metal (bool): Flag that indicates the material being studied is a
                metal, which changes the smearing from Gaussian to second order
                Methfessel-Paxton of 0.2 eV.
            in_custodian (bool): Flag that indicates whether the calculation should be
                run inside a Custodian.
            number_nodes (int): Number of nodes that should be used for the calculations.
                Is required to add the proper `_category` to the Firework generated, so
                it is picked up by the right Fireworker.
            fw_action (fireworks.FWAction): FWAction to return after the final
                PulayTask is completed.
        """

        # Create the PyTask that sets up the calculation
        setup_optimize = PyTask(
            func="pybat.cli.commands.setup.optimize",
            kwargs={"structure": structure,
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
                               fw_action=fw_action)

        # Only add number of nodes to spec if specified
        firework_spec = {}
        if number_nodes is None or number_nodes == 0:
            firework_spec.update({"_category": "none"})
        else:
            firework_spec.update({"_category": str(number_nodes) + "nodes"})

        # Combine the FireTasks into one FireWork
        super(PybatOptimizeFW, self).__init__(
            tasks=[setup_optimize, vasprun, get_cathode, pulay_task],
            name="Geometry optimization", spec=firework_spec
        )


class NebFirework(Firework):

    def __init__(self, directory, nimages, functional, is_metal=False, is_migration=False,
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
        firework_spec = {}
        if number_nodes is None or number_nodes == 0:
            firework_spec.update({"_category": "none"})
        else:
            firework_spec.update({"_category": str(number_nodes) + "nodes"})

        # Combine the two FireTasks into one FireWork
        super(NebFirework, self).__init__(
            tasks=[setup_neb, vasprun], name="NEB calculation", spec=firework_spec
        )


# region * Token Fireworks for testing


class FirstFirework(Firework):

    def __init__(self, final_message):
        super(FirstFirework, self).__init__(
            tasks=[ScriptTask.from_str("echo 'Here we go!'"),
                   MiddleTask(message="next",
                              fw_action=FWAction(additions=Firework([ScriptTask.from_str(
                                  "echo '" + final_message + "'")])))])

# endregion
