# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License


import os
import shutil
import warnings

from fireworks import LaunchPad, FWorker
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

"""
Set of scripts for configuring the pybat package.
"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "pre-alpha"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Mar 2019"


def launchpad(launchpad_file=None, database="base"):
    """
    Script to set up the configuration of the launchpad for accessing the workflow
    server.

    Args:
        launchpad_file (str): my_launchpad.yaml file from which to load the mongoDB database
            details.

    Returns:
        None

    """

    if launchpad_file:
        lpad = LaunchPad.from_file(launchpad_file)
    else:
        host = input("Please provide the server host: ")
        port = int(input("Please provide the port number: "))
        name = input("Please provide the server name: ")
        username = input("Please provide your username: ")
        password = input("Please provide your password: ")

        lpad = LaunchPad(
            host=host, port=port, name=name, username=username, password=password,
            ssl=True, authsource="admin"
        )

    # Test the launchpad
    print("\nAttempting connection to mongoDB database...")
    _ = lpad.get_fw_ids()
    print("Connection successful!\n")

    config_lpad_file = os.path.join(os.path.expanduser("~"), ".pybat_config",
                                    "launchpad", database + "_launchpad.yaml")
    try:
        os.makedirs(
            os.path.join(os.path.expanduser("~"), ".pybat_config", "launchpad")
        )
    except FileExistsError:
        pass

    lpad.to_file(config_lpad_file)
    print("Launchpad file written to " + config_lpad_file + "\n")


def fworker(fireworker_file=None, fworker_name="base"):
    """
    Script to set up the configuration of the fireworker.

    Returns:
        None

    """

    if fireworker_file:
        fireworker = FWorker.from_file(fireworker_file)
    else:
        name = input("Please provide the fireworker name: ")
        vasp_cmd = input("Please provide the full vasp command: ")
        fireworker = FWorker(name=name, env={"vasp_cmd": vasp_cmd})

    # Make sure the fireworker has the required node categories.
    if fireworker.category:
        print("\nNote: Pybat will overwrite the category in the fireworker file to make "
              "sure that the jobs submitted on this fireworker only pick up the "
              "Fireworks with the correct category setting, i.e. corresponding to a "
              "number of nodes.\n")
    fireworker.category = ["none", "1nodes"]

    try:
        os.makedirs(
            os.path.join(os.path.expanduser("~"), ".pybat_config", "fworker")
        )
    except FileExistsError:
        pass

    config_fw_file = os.path.join(os.path.expanduser("~"), ".pybat_config",
                                  "fworker", fworker_name + "_fworker.yaml")

    fireworker.to_file(config_fw_file)
    print("Fireworker file written to " + config_fw_file + "\n")


def qadapter(qadapter_file=None, fworker_name="base"):
    """
    Script to set up the configuration of the queueing system. Note that we store the
    queue_adapter in the same configuration directory as the Fireworker, i.e. fworker.
    This is in the assumption that each worker has one preferred queueing system.

    Returns:
        None

    """
    if qadapter_file:
        queue_adapter = CommonAdapter.from_file(qadapter_file)
    else:
        logdir = input("Please provide the path to the directory where the log "
                       "should be stored. \n"
                       "(Note: if the directory does not exist, it will be created): ")
        logdir = os.path.abspath(logdir)
        if not os.path.exists(logdir):
            os.makedirs(logdir)
        queue_adapter = CommonAdapter.from_file(os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "examples", "config",
            "fworker", "example_qadapter.yaml"
        ))
        print("\nNote: 'rocket_launch' has been set to an infinite rapidfire mode.")

    try:
        os.makedirs(
            os.path.join(os.path.expanduser("~"), ".pybat_config", "fworker")
        )
    except FileExistsError:
        pass

    config_q_file = os.path.join(os.path.expanduser("~"), ".pybat_config",
                                 "fworker", fworker_name + "_qadapter.yaml")

    template_file = os.path.join(
        os.path.expanduser("~"), ".pybat_config", "fworker",
        fworker_name + "_job_template.sh"
    )
    if not os.path.exists(template_file):
        print()
        warn_message = "No corresponding template file found! Don't forget to use " \
                       "pybat config jobscript to add the template file of the '" \
                       + fworker_name + "' fireworker."
        warnings.warn(warn_message)

    queue_adapter.template_file = template_file
    queue_adapter.to_file(config_q_file)
    print("\nQueue adapter file written to " + config_q_file)


def jobscript(template_file, fworker_name="base"):
    """
    Add a template jobscript for submitting jobs to a queueing system. Should be
    adaptable by a queue_adapter, i.e. it should contain variable such as $${nnodes},
    $${rocket_launch}, etc... Again, we'll assume for now that one template per
    fireworker is sufficient.

    Args:
        template_file:
        fworker_name:

    Returns:

    """

    try:
        os.makedirs(
            os.path.join(os.path.expanduser("~"), ".pybat_config", "fworker")
        )
    except FileExistsError:
        pass

    config_template_file = os.path.join(os.path.expanduser("~"), ".pybat_config",
                                        "fworker", fworker_name + "_job_template.sh")

    shutil.copy(template_file, config_template_file)
    print("\nCopied job template to " + config_template_file + "\n")


def check():
    config_dir = os.path.join(os.path.expanduser("~"), ".pybat_config")

    if not os.path.exists(config_dir):
        print("\nNo configuration directory found! Use 'pybat config' to set up the "
              "pybat configuration.")
        return

    if not os.path.exists(os.path.join(config_dir, "launchpad")):
        print("\nNo launchpad directory found. Use 'pybat config launchpad' to set up a "
              "launchpad configuration.")
    else:
        print("\nConfig launchpads" + "\n-----------------")
        lpads = os.listdir(os.path.join(config_dir, "launchpad"))

        for lpad_file in lpads:
            name = lpad_file.split("_")[0]
            lpad = load_config("launchpad", name)
            print(name + ": (host) " + lpad.host)
            print(len(name) * " " + "  (name) " + lpad.name)
            print(len(name) * " " + "  (user) " + lpad.username)

    if not os.path.exists(os.path.join(config_dir, "fworker")):
        print("\nNo fworker directory found. Use 'pybat config fworker' to set up a "
              "fworker configuration.")
    else:
        print("\nConfig fireworkers" + "\n------------------")
        fworker_files = [f.split("_") for f
                         in os.listdir(os.path.join(config_dir, "fworker"))]
        fworker_dict = {}

        for file in fworker_files:
            if file[0] not in fworker_dict.keys():
                fworker_dict[file[0]] = []
            fworker_dict[file[0]].append("_".join(file[1:]))

        for f, l in fworker_dict.items():
            needed_files = ["fworker.yaml", "qadapter.yaml", "job_template.sh"]
            t = [i if i in needed_files else "NONE" for i in needed_files]
            print(f + ": " + t[0])
            for el in t[1:]:
                print(len(f) * " " + "  " + el)
        print()







def load_config(config, name="base"):
    """

    Args:
        config:
        name:

    Returns:

    """
    try:
        if config == "launchpad":
            return LaunchPad.from_file(
                os.path.join(os.path.expanduser("~"), ".pybat_config",
                             "launchpad", name + "_launchpad.yaml")
            )
        if config == "fworker":
            return FWorker.from_file(
                os.path.join(os.path.expanduser("~"), ".pybat_config",
                             "fworker", name + "_fworker.yaml")
            )
        if config == "qadapter":
            return CommonAdapter.from_file(
                os.path.join(os.path.expanduser("~"), ".pybat_config",
                             "fworker", name + "_qadapter.yaml")
            )
    except FileNotFoundError:
        raise FileNotFoundError(
            "Did not find the corresponding configuration file in " \
            + os.path.join(os.path.expanduser("~"), ".pybat_config") \
            + ". Use 'pybat config " + config + "' to set up the " + name + \
            " configuration for the " + config + "."
        )
