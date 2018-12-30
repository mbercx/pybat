# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License


import os

from pathlib import Path
from ruamel.yaml import YAML

"""
Set of scripts for configuring the pybat package.
"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "May 2018"

# Load the workflow configuration
CONFIG_FILE = os.path.join(os.path.expanduser("~"), ".pybat_wf_config.yaml")


def lpad(launchpad_file=None):
    """
    Script to set up the configuration of the launchpad for accessing the workflow
    server.

    Args:
        launchpad_file (str): my_launchpad.yaml file from which to load the mongoDB database
            details.

    Returns:
        None

    """
    yaml = YAML()
    yaml.default_flow_style = False

    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, 'r') as config_file:
            config_dict = yaml.load(config_file.read())
    else:
        config_dict = {"SERVER": {}, "WORKFLOW": {}}

    if launchpad_file:
        with open(launchpad_file, 'r') as launchpad_file:
            config_dict["SERVER"].update(yaml.load(launchpad_file.read()))
    else:
        config_dict["SERVER"]["host"] = input("Please provide the server "
                                              "host: ")
        config_dict["SERVER"]["port"] = input("Please provide the port "
                                              "number: ")
        config_dict["SERVER"]["name"] = input("Please provide the server "
                                              "name: ")
        config_dict["SERVER"]["username"] = input("Please provide your "
                                                  "username: ")
        config_dict["SERVER"]["password"] = input("Please provide your "
                                                  "password: ")
        # TODO Add server check

    with Path(CONFIG_FILE) as config_file:
        yaml.dump(config_dict, config_file)


def script(script_path=None):
    """
    Script to set up the configuration of the workflow jobscript.

    Returns:
        None

    """
    yaml = YAML()
    yaml.default_flow_style = False

    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, 'r') as config_file:
            config_dict = yaml.load(config_file.read())
    else:
        config_dict = {"SERVER": {}, "WORKFLOW": {}}

    if not script_path:
        script_path = input(
            "Please provide the full path to the workflow script: "
        )
        while not os.path.exists(script_path):

            script_path = input(
                "Provided path does not exist. Please provide the full path to the "
                "workflow script again: "
            )

        if not os.path.isabs(script_path):

            print("Provided path is not an absolute path. Finding absolute "
                  "path for proper configuration of the workflows...")
            script_path = os.path.abspath(script_path)

    else:
        script_path = os.path.abspath(script_path)

    config_dict["WORKFLOW"]["script_path"] = script_path

    with Path(CONFIG_FILE) as config_file:
        yaml.dump(config_dict, config_file)


def base(settings="all"):
    """
    Script to set up the configuration of the workflow server and jobscripts.

    Returns:
        None

    """
    yaml = YAML()
    yaml.default_flow_style = False

    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, 'r') as config_file:
            config_dict = yaml.load(config_file.read())
    else:
        config_dict = {"SERVER": {}, "WORKFLOW": {}}

    if settings in ["server", "all"]:
        config_dict["SERVER"]["host"] = input("Please provide the server "
                                              "host: ")
        config_dict["SERVER"]["port"] = input("Please provide the port "
                                              "number: ")
        config_dict["SERVER"]["name"] = input("Please provide the server "
                                              "name: ")
        config_dict["SERVER"]["username"] = input("Please provide your "
                                                  "username: ")
        config_dict["SERVER"]["password"] = input("Please provide your "
                                                  "password: ")
        # TODO Add server check

    if settings in ["workflow", "all"]:

        script_path = input(
            "Please provide the full path to the workflow script: "
        )
        while not os.path.exists(script_path):

            script_path = input(
                "Provided path does not exist. Please provide the full path to the "
                "workflow script again: "
            )

        if not os.path.isabs(script_path):

            print("Provided path is not an absolute path. Finding absolute "
                  "path for proper configuration of the workflows...")
            script_path = os.path.abspath(script_path)

        config_dict["WORKFLOW"]["script_path"] = script_path

    with Path(CONFIG_FILE) as config_file:
        yaml.dump(config_dict, config_file)