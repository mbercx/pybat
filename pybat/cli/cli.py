# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os

import click
from fireworks import LaunchPad
from pymongo.errors import ServerSelectionTimeoutError
from ruamel.yaml import YAML

from pybat.core import Cathode
from pybat.core import LiRichCathode

"""
Command line interface for the pybat package.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "pre-alpha"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Mar 2019"

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

# This is used to make '-h' a shorter way to access the CLI help
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

# region * Help strings for options
# Because several options have the same help string, it's easier to gather those here,
# so that in case we adjust them it only has to be done once.

DISTANCE_HELP = "Distance between the oxygen pairs for the intialization of the dimer " \
                "structure. If not specified, the user will be requested to specify " \
                "the value when running the script."

IN_CUSTODIAN_HELP = "Run the calculations in a Custodian for automatic error handling."

IS_METAL_HELP = "Flag to indicate that the structure is metallic. This will make the " \
                "algorithm choose Methfessel-Paxton smearing of 0.2 eV."

FUNCTIONAL_HELP = "Option for configuring the functional used in the calculation. User " \
                  "must provide the functional information in the form of a single " \
                  "string, starting with the string that determines the functional, " \
                  "then with string/float pairs for specifying further settings. " \
                  "Defaults to 'pbe'. \nExamples:\n" \
                  "* 'pbeu Mn\xa03.9 V 3.1' ~ PBE+U Dudarev approach) with effective U " \
                  "equal to 3.9 for Mn and 3.1 for V.\n" \
                  "* 'scan' ~ SCAN\n" \
                  "* 'hse' ~ HSE06\n" \
                  "*\xa0'hse\xa0hfscreen\xa00.3'\xa0~\xa0HSE03\n"

MIGRATION_INDICES_HELP = "Starting and final indices of a migration. Provided with two " \
                         "integers when using the option, i.e. -I 4 12."

NUMBER_NODES_HELP = "Number of nodes that should be used for the calculations. Is " \
                    "required to add the proper `_category` to the Firework generated, " \
                    "so it is picked up by the right Fireworker."

WRITE_CIF_HELP = "Flag that indicates that the structure(s) should also be written as " \
                 "a .cif file."


# endregion

@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    """
    CLI tools for performing calculations for studying batteries.
    """
    pass


@main.command(context_settings=CONTEXT_SETTINGS)
@click.option("-L", "-lpad_name", default="base")
@click.option("-F", "--fworker_name", default="base")
def qlaunch(lpad_name, fworker_name):
    """
    Launch jobs to the queue that will accept Fireworks.

    Returns:

    """
    from fireworks.queue.queue_launcher import launch_rocket_to_queue
    from pybat.cli.commands.config import load_config

    fireworker = load_config("fworker", fworker_name)
    qadapter = load_config("qadapter", fworker_name)

    launch_rocket_to_queue(
        launchpad=load_config("launchpad", lpad_name), fworker=fireworker,
        qadapter=qadapter, launcher_dir='.', create_launcher_dir=False, fill_mode=True
    )

# TODO Add checks for U-value input

# region * Config


@main.group(context_settings=CONTEXT_SETTINGS)
def config():
    """
    Configure the Workflows setup.

    """
    pass


@config.command(context_settings=CONTEXT_SETTINGS)
@click.option("-l", "--launchpad_file", default="")
@click.option("-N", "--name", default="base")
def launchpad(launchpad_file, name):
    """
    Configure a Workflows server or launchpad.

    """
    if launchpad_file == "":
        launchpad_file = None
    from pybat.cli.commands.config import launchpad
    launchpad(launchpad_file=launchpad_file, database=name)


@config.command(context_settings=CONTEXT_SETTINGS)
@click.option("-F", "--fworker_file", default="")
@click.option("-N", "--name", default="base")
def fworker(fworker_file, name):
    """
    Configure the base settings of a fireworker.

    """
    if fworker_file == "":
        fworker_file = None
    from pybat.cli.commands.config import fworker
    fworker(fireworker_file=fworker_file, fworker_name=name)


@config.command(context_settings=CONTEXT_SETTINGS)
@click.argument("qadapter_file", nargs=1)
@click.option("-N", "--name", default="base")
def queue(qadapter_file, name):
    """
    Configure the standard queue adapter of a fireworker.

    """
    from pybat.cli.commands.config import queue

    queue(qadapter_file=qadapter_file, fworker_name=name)


@config.command(context_settings=CONTEXT_SETTINGS)
@click.argument("template_file", nargs=1)
@click.option("-N", "--name", default="base")
def jobscript(template_file, name):
    """
    Configure the standard queue adapter of a fireworker.

    """
    from pybat.cli.commands.config import jobscript

    jobscript(template_file=template_file, fworker_name=name)


# endregion

# region * Define


@main.group(context_settings=CONTEXT_SETTINGS)
def define():
    """
    Define a structural change.

    """
    pass


@define.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--migration_indices", "-I", default=(0, 0),
              help=MIGRATION_INDICES_HELP)
@click.option("--write_cif", "-w", is_flag=True,
              help=WRITE_CIF_HELP)
def migration(structure_file, migration_indices, write_cif):
    """
    Define a migration of an ion in a structure.

    """
    from pybat.cli.commands.define import define_migration

    cat = Cathode.from_file(structure_file)

    define_migration(structure=cat,
                     migration_indices=migration_indices,
                     write_cif=write_cif)


@define.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--dimer_indices", "-i", default=(0, 0))
@click.option("--distance", "-d", default=float(0))
@click.option("--remove_cations", "-r", is_flag=True)
@click.option("--write_cif", "-w", is_flag=True,
              help=WRITE_CIF_HELP)
def dimer(structure_file, dimer_indices, distance, remove_cations, write_cif):
    """
    Define the formation of a dimer in a structure.

    """
    from pybat.cli.commands.define import define_dimer

    cat = LiRichCathode.from_file(structure_file)

    define_dimer(structure=cat,
                 dimer_indices=dimer_indices,
                 distance=distance,
                 remove_cations=remove_cations,
                 write_cif=write_cif)


# endregion

# region * Get


@main.group(context_settings=CONTEXT_SETTINGS)
def get():
    """
    Obtain data from output files.

    """
    pass


@get.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".")
@click.option("--write_cif", "-w", is_flag=True,
              help=WRITE_CIF_HELP)
def structure(directory, write_cif):
    """
    Obtain the structure with its magnetic configuration.

    """
    from pybat.cli.commands.get import get_structure

    get_structure(directory=directory,
                  write_cif=write_cif)


@get.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".",
              help="The directory which contains the required data for the "
                   "final cathode JSON file input. This includes initial_cathode.json "
                   "and the output of a VASP calculation, i.e. the CONTCAR and OUTCAR "
                   "files. ")
@click.option("--to_current_dir", "-c", is_flag=True,
              help="Flag to indicate that the final cathode file should be "
                   "witten to the current directory instead of the directory "
                   "which contains the data.")
@click.option("--ignore_magmom", "-i", is_flag=True,
              help="Flag that indicates that the final magnetic moments should be "
                   "ignored, i.e. that the magnetic moments should be kept as the ones "
                   "in initial_cathode.json.")
@click.option("--write_cif", "-w", is_flag=True,
              help=WRITE_CIF_HELP)
def cathode(directory, to_current_dir, ignore_magmom, write_cif):
    """
    Obtain the Cathode with its magnetic configuration and vacancies.

    """
    from pybat.cli.commands.get import get_cathode

    get_cathode(directory=directory,
                to_current_dir=to_current_dir,
                ignore_magmom=ignore_magmom,
                write_cif=write_cif)


@get.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".")
@click.option("--method", "-M", default="pymatgen")
def barrier(directory, method):
    """
    Combine the images of a NEB calculation to show the transition.
    """
    from pybat.cli.commands.get import get_barrier

    get_barrier(directory=directory,
                method=method)


@get.command(context_settings=CONTEXT_SETTINGS)
def voltage():
    """
    Calculate the voltage for a battery cathode versus a Li anode.
    """
    pass
    # TODO


@get.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".")
def endiff(directory):
    """
    Calculate the energy difference for a transition.
    """
    from pybat.cli.commands.get import get_endiff

    get_endiff(directory)


# endregion

# region * Setup


@main.group(context_settings=CONTEXT_SETTINGS)
def setup():
    """
    Set up calculations.
    """
    pass


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.option("--functional", "-f", default="pbe",
              help="Option for configuring the functional used in the calculation. "
                   "User must provide the functional information in the form of a "
                   "single string, starting with the string that determines the "
                   "functional, then with string/float pairs for specifying further "
                   "settings. Defaults to 'pbe'. Examples:\n"
                   "* 'pbeu Mn\xa03.9 V 3.1' ~ PBE+U (Dudarev approach) with effective "
                   "U equal to 3.9 for Mn and 3.1 for V.\n"
                   "* 'hse' ~ HSE06\n"
                   "*\xa0'hse\xa0hfscreen\xa00.3'\xa0~\xa0HSE03\n"
              )
@click.argument("structure_file", nargs=1)
@click.option("--calculation_dir", "-d", default="",
              help="The directory in which to set up the calculation. "
                   "Default is FUNCTIONAL_static.")
@click.option("--write_chgcar", "-C", is_flag=True)
def static(structure_file, functional, calculation_dir, write_chgcar):
    """
    Set up a static calculation for a structure.
    """
    from pybat.cli.commands.setup import static

    cat = Cathode.from_file(structure_file)

    static(structure=cat,
           functional=string_to_functional(functional),
           calculation_dir=calculation_dir,
           write_chgcar=write_chgcar)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--functional", "-f", default="pbe",
              help="Option for configuring the functional used in the calculation. "
                   "User must provide the functional information in the form of a "
                   "single string, starting with the string that determines the "
                   "functional, then with string/float pairs for specifying further "
                   "settings. Defaults to 'pbe'. Examples:\n"
                   "* 'pbeu Mn\xa03.9 V 3.1' ~ PBE+U (Dudarev approach) with effective "
                   "U equal to 3.9 for Mn and 3.1 for V.\n"
                   "* 'hse' ~ HSE06\n"
                   "*\xa0'hse\xa0hfscreen\xa00.3'\xa0~\xa0HSE03\n"
              )
@click.option("--calculation_dir", "-d", default="",
              help="The directory in which to set up the calculation. "
                   "Default is FUNCTIONAL_.")
@click.option("--is_metal", "-m", is_flag=True,
              help="Flag to indicate that the structure is metallic. This "
                   "will make the algorithm choose Methfessel-Paxton "
                   "smearing of 0.2 eV.")
def optimize(structure_file, functional, calculation_dir, is_metal):
    """
    Set up a geometry optimization for a structure.
    """
    from pybat.cli.commands.setup import optimize

    cat = Cathode.from_file(structure_file)

    optimize(structure=cat,
             functional=string_to_functional(functional),
             calculation_dir=calculation_dir,
             is_metal=is_metal)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".",
              help="Directory in which to set up the calculations for the "
                   "first step in the transition path determination. Note "
                   "that this directory has to contain the structure files "
                   "for the initial and final state. ")
@click.option("--functional", "-f", default="pbe",
              help="Option for configuring the functional used in the calculation. "
                   "User must provide the functional information in the form of a "
                   "single string, starting with the string that determines the "
                   "functional, then with string/float pairs for specifying further "
                   "settings. Defaults to 'pbe'. Examples:\n"
                   "* 'pbeu Mn\xa03.9 V 3.1' ~ PBE+U (Dudarev approach) with effective "
                   "U equal to 3.9 for Mn and 3.1 for V.\n"
                   "* 'hse' ~ HSE06\n"
                   "*\xa0'hse\xa0hfscreen\xa00.3'\xa0~\xa0HSE03\n"
              )
@click.option("--is_metal", "-m", is_flag=True,
              help="Flag to indicate that the structure is metallic. This "
                   "will make the algorithm choose Methfessel-Paxton "
                   "smearing of 0.2 eV.")
@click.option("--is_migration", "-M", is_flag=True,
              help="Flag to designate the transition as a migration. "
                   "Activating this flag means that a static calculation will "
                   "be set up to determine the charge density of the host "
                   "structure, i.e. without the migrating ion. This charge "
                   "density can then be used to find a good initial guess "
                   "for the migration pathway.")
@click.option("--optimize_initial", "-O", is_flag=True,
              help="Flag to indicate that the initial structure should also be "
                   "optimized.")
def transition(directory, functional, is_metal, is_migration, optimize_initial):
    """
    Set up a the geometry optimizations for the initial and final state of a
    transition.
    """
    from pybat.cli.commands.setup import transition

    transition(directory=directory,
               functional=string_to_functional(functional),
               is_metal=is_metal,
               is_migration=is_migration,
               optimize_initial=optimize_initial)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".",
              help="Directory in which to set up the calculations for the "
                   "first step in the transition path determination. Note "
                   "that this directory has to contain the structure files "
                   "for the initial and final state.")
@click.option("--functional", "-f", default="pbe",
              help="Option for configuring the functional used in the calculation. "
                   "User must provide the functional information in the form of a "
                   "single string, starting with the string that determines the "
                   "functional, then with string/float pairs for specifying further "
                   "settings. Defaults to 'pbe'. Examples:\n"
                   "* 'pbeu Mn\xa03.9 V 3.1' ~ PBE+U (Dudarev approach) with effective "
                   "U equal to 3.9 for Mn and 3.1 for V.\n"
                   "* 'hse' ~ HSE06\n"
                   "*\xa0'hse\xa0hfscreen\xa00.3'\xa0~\xa0HSE03\n"
              )
@click.option("--nimages", "-N", default=7,
              help="Number of images. Defaults to 7.")
@click.option("--is_metal", "-m", is_flag=True,
              help="Flag to indicate that the structure is metallic. This "
                   "will make the algorithm choose Methfessel-Paxton "
                   "smearing of 0.2 eV.")
@click.option("--is_migration", "-M", is_flag=True,
              help="Flag to designate the transition as a migration. "
                   "Activating this flag means that a static calculation will "
                   "be set up to determine the charge density of the host "
                   "structure, i.e. without the migrating ion. This charge "
                   "density can then be used to find a good initial guess "
                   "for the migration pathway.")
def neb(directory, functional, nimages, is_metal, is_migration):
    """
    Set up the Nudged Elastic Band calculation based on the output in the
    initial and final directory.

    Returns:

    """
    from pybat.cli.commands.setup import neb

    neb(directory=directory,
        functional=string_to_functional(functional),
        nimages=nimages,
        is_metal=is_metal,
        is_migration=is_migration)


# endregion

# region * Utility


@main.group(context_settings=CONTEXT_SETTINGS)
def util():
    """
    Utility command for the pybat package.

    """
    pass


@util.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".")
@click.option("--filename", "-f", default="neb_result.cif")
def showpath(directory, filename):
    """
    Combine the images of a NEB calculation to show the transition.

    """
    from pybat.cli.commands.util import show_path

    show_path(directory=directory,
              filename=filename)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--file_format", "-F", default="json")
def conv(structure_file, file_format):
    """
    Convert a structure into the conventional unit cell.

    """
    from pybat.cli.commands.util import conventional_structure

    conventional_structure(structure_file=structure_file,
                           fmt=file_format)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument("vasprun_file", nargs=1)
def data(vasprun_file):
    """
    Compress the data of the vasprun.xml file to a JSON file.

    """
    from pybat.cli.commands.util import data

    data(vasprun_file=vasprun_file)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--file_format", "-F", default="json")
def prim(structure_file, file_format):
    """
    Convert a structure into the primitive unit cell.

    """
    from pybat.cli.commands.util import primitive_structure

    primitive_structure(structure_file=structure_file,
                        fmt=file_format)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument("cell", nargs=1)
@click.argument("structure_file", nargs=1)
@click.option("--file_format", "-F", default="json")
def supercell(cell, structure_file, file_format):
    """
    Convert a structure to a supercell.

    """
    from pybat.cli.commands.util import make_supercell

    make_supercell(structure_file=structure_file,
                   supercell=cell,
                   fmt=file_format)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
def show(structure_file):
    """
    Print the cathode details to the screen.

    """
    from pybat.cli.commands.util import print_structure

    print_structure(structure_file=structure_file)


# endregion

# region * Workflow

@main.group(context_settings=CONTEXT_SETTINGS)
def workflow():
    """
    Scripts for setting up workflows and submitting them to the server.
    """
    pass


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--functional", "-f", default="pbe",
              help=FUNCTIONAL_HELP
              )
@click.option("--directory", "-d", default="",
              help="Directory in which the static calculation will be performed.")
@click.option("--write_chgcar", "-C", is_flag=True,
              help="Flag that indicates that the CHGCAR should be written, which "
                   "is not done by default to conserve space.")
@click.option("--in_custodian", "-c", is_flag=True,
              help=IN_CUSTODIAN_HELP)
@click.option("--number_nodes", "-n", default=0,
              help=NUMBER_NODES_HELP)
def static(structure_file, functional, directory, write_chgcar, in_custodian,
           number_nodes):
    """
    Set up an static calculation workflow.
    """
    from pybat.workflow.workflows import get_wf_static

    # Process the input options
    if number_nodes == 0:
        number_nodes = None
    functional = string_to_functional(functional)

    # Set up the calculation directory
    directory = set_up_directory(directory, functional, "static")

    cat = Cathode.from_file(structure_file)

    LAUNCHPAD.add_wf(
        get_wf_static(structure=cat,
                      functional=functional,
                      directory=directory,
                      write_chgcar=write_chgcar,
                      in_custodian=in_custodian,
                      number_nodes=number_nodes)
    )


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--functional", "-f", default="pbe", help=FUNCTIONAL_HELP)
@click.option("--directory", "-d", default="")
@click.option("--is_metal", "-m", is_flag=True, help=IS_METAL_HELP)
@click.option("--in_custodian", "-c", is_flag=True)
@click.option("--number_nodes", "-n", default=0, help=NUMBER_NODES_HELP)
def optimize(structure_file, functional, directory, is_metal, in_custodian, number_nodes):
    """
    Set up a geometry optimization workflow.
    """
    from pybat.workflow.workflows import get_wf_optimize

    cat = Cathode.from_file(structure_file)

    # Process the input options
    if number_nodes == 0:
        number_nodes = None
    functional = string_to_functional(functional)

    # Set up the calculation directory
    directory = set_up_directory(directory, functional, "optimize")

    LAUNCHPAD.add_wf(
        get_wf_optimize(structure=cat,
                        functional=functional,
                        directory=directory,
                        is_metal=is_metal,
                        in_custodian=in_custodian,
                        number_nodes=number_nodes)
    )


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--sub_sites", "-s", default="None",
              help="Indices of the sites that should be substituted to generate the "
                   "possible configurations. Can either be a List of integers or a "
                   "range()."
                   "\nExamples:\n"
                   "\"range(4)\"\n"
                   "\"[0, 2, 4]\"")
@click.option("--element_list", "-E", default="None",
              help="List of elements that should be placed on the substitution sites to "
                   "generate the configurations."
                   "\nExamples:\n"
                   "\"['Ni', 'Mn', 'Co']\"\n"
                   "\"['Li', 'Vac']\"")
@click.option("--sizes", "-S", default="None",
              help="Allowed unit cell sizes for the generation of configurations. Can "
                   "be either a List of numbers or a range(). \nExamples:\n"
                   "\"range(1, 7)\" \n"
                   "\"[1,2,5]\"")
@click.option("--conc_restrict", "-R", default="None",
              help="Dictionary of the allowed concentration ranges for each element. "
                   "Note that the concentration is defined versus the total amount of "
                   "atoms in the unit cell. \n Examples:\n \"{'Li':(0.2, 0.3)}\"\n "
                   "\"{'Ni':(0.1, 0.2), 'Mn':(0, 0.1)}\"")
@click.option("--max_conf", "-M", default=0,
              help="Maximum number of configurations to generate. Note that this "
                   "is on top of the number of already existing configurations in the "
                   "directory tree in case the existing configurations "
                   "are not included using the include_existing option.")
@click.option("--functional", "-f", default="pbe", help=FUNCTIONAL_HELP)
@click.option("--directory", "-d", default="",
              help="Directory in which to set up the configuration workflow.")
@click.option("--include_existing", "-X", is_flag=True,
              help="Include the existing configurations in the directory tree for the "
                   "calculations. This means that the energy will also be calculated "
                   "for configuration.json files found in the directory tree.")
@click.option("--in_custodian", "-c", is_flag=True, help=IN_CUSTODIAN_HELP)
@click.option("--number_nodes", "-n", default=0, help=NUMBER_NODES_HELP)
def configuration(structure_file, functional, sub_sites, element_list, sizes,
                  directory, include_existing, conc_restrict, max_conf, in_custodian,
                  number_nodes):
    """
    Set up a workflow for a set of configurations.

    Calculate the geometry and total energy of a set of configurations. The
    configurations can be specified using the various options. In case some required
    input is not specified during command execution, it will be requested from the user.

    """
    from pybat.workflow.workflows import get_wf_configurations

    cat = Cathode.from_file(structure_file)

    # Set up the directory
    if directory == "":
        directory = os.getcwd()
    directory = os.path.abspath(directory)

    # Process the option input
    try:
        substitution_sites = eval(sub_sites)
    except SyntaxError:
        substitution_sites = [int(site) for site in sub_sites.split(" ")]
    try:
        element_list = eval(element_list)
    except SyntaxError:
        element_list = [el for el in element_list.split(" ")]
    try:
        sizes = [int(i) for i in sizes.strip("[]").split(",")]
    except ValueError:
        sizes = eval(sizes)
    try:
        conc_restrict = eval(conc_restrict)
    except TypeError:
        conc_restrict = None

    # In case some required settings are missing, request them
    if not substitution_sites:
        print(cat)
        print()
        substitution_sites = [int(i) for i in input(
            "Please provide the substitution site indices, separated by a space: "
        ).split(" ")]
    if not element_list:
        element_list = [i for i in input(
            "Please provide the substitution elements, separated by a space: "
        ).split(" ")]
    if not sizes:
        sizes = [int(i) for i in input(
            "Please provide the possible unit cell sizes, separated by a space: "
        ).split(" ")]

    LAUNCHPAD.add_wf(
        get_wf_configurations(structure=cat,
                              substitution_sites=substitution_sites,
                              element_list=element_list,
                              sizes=sizes,
                              concentration_restrictions=conc_restrict,
                              max_configurations=max_conf,
                              functional=string_to_functional(functional),
                              directory=directory,
                              include_existing=include_existing,
                              in_custodian=in_custodian,
                              number_nodes=number_nodes)
    )


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--dimer_indices", "-i", default=(0, 0))
@click.option("--distance", "-d", default=float(0), help=DISTANCE_HELP)
@click.option("--functional", "-f", default="pbe", help=FUNCTIONAL_HELP)
@click.option("--is_metal", "-m", is_flag=True, help=IS_METAL_HELP)
@click.option("--in_custodian", "-c", is_flag=True, help=IN_CUSTODIAN_HELP)
@click.option("--number_nodes", "-n", default=0, help=NUMBER_NODES_HELP)
def dimer(structure_file, dimer_indices, distance, functional, is_metal,
          in_custodian, number_nodes):
    """
    Set up dimer calculation workflows.
    """
    from pybat.workflow.workflows import get_wf_dimer

    cat = LiRichCathode.from_file(structure_file)

    LAUNCHPAD.add_wf(
        get_wf_dimer(structure=cat,
                     dimer_indices=dimer_indices,
                     distance=distance,
                     functional=string_to_functional(functional),
                     is_metal=is_metal,
                     in_custodian=in_custodian,
                     number_nodes=number_nodes)
    )


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument("directory", nargs=1)
@click.option("--nimages", "-N", default=7, show_default=True,
              help="Number of images for the NEB calculation. Defaults to 7.")
@click.option("--functional", "-f", default="pbe",
              help=FUNCTIONAL_HELP)
@click.option("--is_metal", "-m", is_flag=True, help=IS_METAL_HELP)
@click.option("--is_migration", "-M", is_flag=True,
              help="Flag to designate the transition as a migration. "
                   "Activating this flag means that a static calculation will "
                   "be set up to determine the charge density of the host "
                   "structure, i.e. without the migrating ion. This charge "
                   "density can then be used to find a good initial guess "
                   "for the migration pathway.")
@click.option("--in_custodian", "-c", is_flag=True, help=IN_CUSTODIAN_HELP)
@click.option("--number_nodes", "-n", default=0, help=NUMBER_NODES_HELP)
def neb(directory, nimages, functional, is_metal, is_migration, in_custodian,
        number_nodes):
    """
    Set up dimer calculation workflows.
    """
    from pybat.workflow.workflows import get_wf_neb

    LAUNCHPAD.add_wf(
        get_wf_neb(directory=directory,
                   nimages=nimages,
                   functional=string_to_functional(functional),
                   is_metal=is_metal,
                   is_migration=is_migration,
                   in_custodian=in_custodian,
                   number_nodes=number_nodes)
    )


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--distance", "-d", default=float(0), help=DISTANCE_HELP)
@click.option("--functional", "-f", default="pbe", help=FUNCTIONAL_HELP)
@click.option("--is_metal", "-m", is_flag=True, help=IS_METAL_HELP)
@click.option("--in_custodian", "-c", is_flag=True, help=IN_CUSTODIAN_HELP)
@click.option("--number_nodes", "-n", default=0, help=NUMBER_NODES_HELP)
def noneq_dimers(structure_file, distance, functional, is_metal, in_custodian,
                 number_nodes):
    """
    Set up dimer calculations for all nonequivalent dimers in a structure.
    """
    from pybat.workflow.workflows import get_wfs_noneq_dimers

    cat = LiRichCathode.from_file(structure_file)

    for wf in get_wfs_noneq_dimers(structure=cat,
                                   distance=distance,
                                   functional=string_to_functional(functional),
                                   is_metal=is_metal,
                                   in_custodian=in_custodian,
                                   number_nodes=number_nodes):
        LAUNCHPAD.add_wf(wf)


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument("site_index", default=None)
@click.argument("structure_file", nargs=1)
@click.option("--distance", "-d", default=float(0), help=DISTANCE_HELP)
@click.option("--functional", "-f", default="pbe", help=FUNCTIONAL_HELP)
@click.option("--is_metal", "-m", is_flag=True, help=IS_METAL_HELP)
@click.option("--in_custodian", "-c", is_flag=True, help=IN_CUSTODIAN_HELP)
@click.option("--number_nodes", "-n", default=0, help=NUMBER_NODES_HELP)
def site_dimers(site_index, structure_file, distance, functional, is_metal, in_custodian,
                number_nodes):
    """
    Set up dimer calculations for all nonequivalent dimers in a structure.
    """
    from pybat.workflow.workflows import get_wfs_site_dimers

    cat = LiRichCathode.from_file(structure_file)

    for wf in get_wfs_site_dimers(structure=cat,
                                  site_index=site_index,
                                  distance=distance,
                                  functional=string_to_functional(functional),
                                  is_metal=is_metal,
                                  in_custodian=in_custodian,
                                  number_nodes=number_nodes):
        LAUNCHPAD.add_wf(wf)


# endregion

# region * Test

@main.group(context_settings=CONTEXT_SETTINGS)
def test():
    """
    Scripts or methods in the testing phase!
    """
    pass


# endregion


# region * Local scripts


def string_to_functional(dict_string):
    """
    Turn input string for function option into a (string, dictionary) representing
    the functional to be used in calculations.

    Args:
        dict_string: String that provides the string/float pairs, separated by a
            space.

    Returns:
        tuple: Tuple of (String that represents the functional, Dictionary of
            string/float pairs).

    """
    dict_pairs = dict_string.split(" ")

    functional = dict_pairs[0]
    functional_dict = dict(
        zip(dict_pairs[1::2], [float(number) for number in dict_pairs[2::2]])
    )

    # Turn functional string into dict
    if functional == "pbeu":
        pbeu_dict = dict()
        pbeu_dict["LDAUU"] = functional_dict
        return functional, pbeu_dict
    else:
        return functional, functional_dict


def set_up_directory(directory, functional, calculation):
    # Set up the calculation directory
    if directory == "":
        directory = os.path.join(os.getcwd(), functional[0])
        if functional[0] == "pbeu":
            directory += "_" + "".join(k + str(functional[1]["LDAUU"][k]) for k
                                       in functional[1]["LDAUU"].keys())
        directory += "_" + calculation
    else:
        directory = os.path.abspath(directory)

    return directory

# endregion
