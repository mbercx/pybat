# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import click

"""
Command line interface for the pybat package.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "May 2018"

# This is used to make '-h' a shorter way to access the CLI help
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    """
    CLI tools for performing calculations for studying batteries.
    """
    pass


##########
# DEFINE #
##########


@main.group(context_settings=CONTEXT_SETTINGS)
def define():
    """
    Define a structural change.

    """
    pass


@define.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--migration_indices", "-i", default=(0, 0))
@click.option("--write_cif", "-w", is_flag=True)
def migration(structure_file, migration_indices, write_cif):
    """
    Define a migration of an ion in a structure.

    """
    from pybat.cli.commands.define import define_migration

    define_migration(structure_file=structure_file,
                     migration_indices=migration_indices,
                     write_cif=write_cif)


@define.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--dimer_indices", "-i", default=(0, 0))
@click.option("--distance", "-d", default=float(0))
@click.option("--remove_cations", "-R", is_flag=True)
@click.option("--write_cif", "-w", is_flag=True)
def dimer(structure_file, dimer_indices, distance, remove_cations, write_cif):
    """
    Define the formation of a dimer in a structure.

    """
    from pybat.cli.commands.define import define_dimer

    define_dimer(structure_file=structure_file,
                 dimer_indices=dimer_indices,
                 distance=distance,
                 remove_cations=remove_cations,
                 write_cif=write_cif)


#########
# SETUP #
#########


@main.group(context_settings=CONTEXT_SETTINGS)
def setup():
    """
    Set up calculations.
    """
    pass


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--calculation_dir", "-d", default="",
              help="The directory in which to set up the calculation. "
                   "Default is FUNCTIONAL_relax.")
@click.option("--is_metal", "-m", is_flag=True,
              help="Flag to indicate that the structure is metallic. This "
                   "will make the algorithm choose Methfessel-Paxton "
                   "smearing of 0.2 eV.")
@click.option("--hse_calculation", "-H", is_flag=True)
def relax(structure_file, calculation_dir, is_metal, hse_calculation):
    """
    Set up a geometry optimization for a structure.
    """
    from pybat.cli.commands.setup import relax

    relax(structure_file=structure_file,
          calculation_dir=calculation_dir,
          is_metal=is_metal,
          hse_calculation=hse_calculation)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".",
              help="Directory in which to set up the calculations for the "
                   "first step in the transition path determination. Note "
                   "that this directory has to contain the structure files "
                   "for the initial and final state. ")
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
@click.option("--hse_calculation", "-H", is_flag=True)
@click.option("--optimize_initial", "-I", is_flag=True)
def transition(directory, is_metal, is_migration, hse_calculation,
               optimize_initial):
    """
    Set up a the geometry optimizations for the initial and final state of a
    transition.
    """

    from pybat.cli.commands.setup import transition

    transition(directory=directory,
               is_metal=is_metal,
               is_migration=is_migration,
               hse_calculation=hse_calculation,
               optimize_initial=optimize_initial)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("dimer_distance", "-D", default=1.4)
@click.option("--is_metal", "-m", is_flag=True,
              help="Flag to indicate that the structure is metallic. This "
                   "will make the algorithm choose Methfessel-Paxton "
                   "smearing of 0.2 eV.")
@click.option("hse_calculation", "-H", is_flag=True)
def dimers(structure_file, dimer_distance, is_metal, hse_calculation):
    """
    Set up the dimer formation calculations for all nonequivalent dimers. Will
    start from the initial cathode structure and set up the geometry
    optimizations for the various dimer formations in a structure, each in its
    own directory. Will also set up a geometry optimization for the initial
    structure in the directory 'initial'.
    """
    from pybat.cli.commands.setup import dimers

    dimers(structure_file=structure_file,
           dimer_distance=dimer_distance,
           is_metal=is_metal,
           hse_calculation=hse_calculation)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".",
              help="Directory in which to set up the calculations for the "
                   "first step in the transition path determination. Note "
                   "that this directory has to contain the structure files "
                   "for the initial and final state.")
@click.option("--nimages", "-n", default=8,
              help="Number of images.")
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
@click.option("--hse_calculation", "-H", is_flag=True)
def neb(directory, nimages, is_metal, is_migration, hse_calculation):
    """
    Set up the Nudged Elastic Band calculation based on the output in the
    initial and final directory.

    Returns:

    """
    from pybat.cli.commands.setup import neb

    neb(directory=directory,
        nimages=nimages,
        is_metal=is_metal,
        is_migration=is_migration,
        hse_calculation=hse_calculation)


###########
# UTILITY #
###########


@main.group(context_settings=CONTEXT_SETTINGS)
def util():
    """
    Utility command for the pybat package.

    """
    pass


@util.command(context_settings=CONTEXT_SETTINGS)
@click.option("--settings", "-s", default="all",
              help="Choice of settings to adjust. Defaults to \"all\", "
                   "but the user can also choose to only change the settings "
                   "for the \"server\", or the \"workflow\".")
def config(settings):
    """
    Configure the workflows server and script.

    """
    from pybat.cli.commands.util import config
    config(settings=settings)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".")
@click.option("--filename", "-f", default="neb_result")
def showpath(directory, filename):
    """
    Combine the images of a NEB calculation to show the transition.

    """
    from pybat.cli.commands.util import show_path

    show_path(directory=directory,
              filename=filename)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--file_format", "-F", default="cif")
def conv(structure_file, file_format):
    """
    Convert a structure into the conventional unit cell.

    """
    from pybat.cli.commands.util import conventional_structure

    conventional_structure(structure_file=structure_file,
                           fmt=file_format)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--file_format", "-F", default="cif")
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
@click.option("--file_format", "-F", default="cif")
def supercell(cell, structure_file, file_format):
    """
    Convert a structure to a supercell.

    """
    from pybat.cli.commands.util import make_supercell

    make_supercell(structure_file=structure_file,
                   supercell=cell,
                   fmt=file_format)


#######
# GET #
#######v


@main.group(context_settings=CONTEXT_SETTINGS)
def get():
    """
    Obtain data from output files.

    """
    pass


@get.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".")
@click.option("--write_cif", "-w", is_flag=True)
def structure(directory, write_cif):
    """
    Obtain the structure with its magnetic configuration.

    """
    from pybat.cli.commands.get import get_structure

    get_structure(directory=directory,
                  write_cif=write_cif)


@get.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".")
@click.option("--write_cif", "-w", is_flag=True)
def cathode(directory, write_cif):
    """
    Obtain the Cathode with its magnetic configuration and vacancies.

    """
    from pybat.cli.commands.get import get_cathode

    get_cathode(directory=directory,
                write_cif=write_cif)


@get.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".")
def barrier(directory):
    """
    Combine the images of a NEB calculation to show the transition.

    """
    from pybat.cli.commands.get import get_barrier

    get_barrier(directory=directory)


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


############
# WORKFLOW #
############

@main.group(context_settings=CONTEXT_SETTINGS)
def workflow():
    """
    Scripts for setting up workflows and submitting them to the server.
    """
    pass


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--directory", "-d", default="")
@click.option("--is_metal", "-m", is_flag=True,
              help="Flag to indicate that the structure is metallic. This "
                   "will make the algorithm choose Methfessel-Paxton "
                   "smearing of 0.2 eV.")
@click.option("--hse_calculation", "-H", is_flag=True)
@click.option("--in_custodian", "-C", is_flag=True)
def relax(structure_file, directory, is_metal, hse_calculation,
          in_custodian):
    """
    Set up a geometry optimization workflow.
    """
    from pybat.workflow import relax_workflow

    relax_workflow(structure_file=structure_file,
                   is_metal=is_metal,
                   directory=directory,
                   hse_calculation=hse_calculation,
                   in_custodian=in_custodian)


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--dimer_indices", "-i", default=(0, 0))
@click.option("--distance", "-d", default=float(0))
@click.option("--is_metal", "-m", is_flag=True,
              help="Flag to indicate that the structure is metallic. This "
                   "will make the algorithm choose Methfessel-Paxton "
                   "smearing of 0.2 eV.")
@click.option("--hse_calculation", "-H", is_flag=True)
@click.option("--in_custodian", "-C", is_flag=True)
def dimer(structure_file, dimer_indices, distance, is_metal,
          hse_calculation, in_custodian):
    """
    Set up dimer calculation workflows.
    """
    from pybat.workflow import dimer_workflow

    dimer_workflow(structure_file=structure_file,
                   dimer_indices=dimer_indices,
                   distance=distance,
                   is_metal=is_metal,
                   hse_calculation=hse_calculation,
                   in_custodian=in_custodian)


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--distance", "-d", default=float(1.4))
@click.option("--is_metal", "-m", is_flag=True,
              help="Flag to indicate that the structure is metallic. This "
                   "will make the algorithm choose Methfessel-Paxton "
                   "smearing of 0.2 eV.")
@click.option("--hse_calculation", "-H", is_flag=True)
@click.option("--in_custodian", "-C", is_flag=True)
def noneq_dimers(structure_file, distance, is_metal, hse_calculation,
                 in_custodian):
    """
    Set up dimer calculations for all nonequivalent dimers in a structure.
    """
    from pybat.workflow import noneq_dimers_workflow

    noneq_dimers_workflow(structure_file=structure_file,
                          is_metal=is_metal,
                          distance=distance,
                          hse_calculation=hse_calculation,
                          in_custodian=in_custodian)


########
# TEST #
########

@main.group(context_settings=CONTEXT_SETTINGS)
def test():
    """
    Scripts or methods in the testing phase!
    """
    pass


@test.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--directory", "-d", default="")
@click.option("--hse_calculation", "-H", is_flag=True)
@click.option("--in_custodian", "-C", is_flag=True)
def workflow(structure_file, directory, hse_calculation, in_custodian):
    """
    Testing for the workflow scripts.
    """
    from pybat.workflow import scf_workflow

    scf_workflow(structure_file=structure_file,
                 directory=directory,
                 hse_calculation=hse_calculation,
                 in_custodian=in_custodian)


@test.command(context_settings=CONTEXT_SETTINGS)
@click.argument("site_index", nargs=1)
@click.argument("structure_file", nargs=1)
@click.option("--distance", "-d", default=float(1.4))
@click.option("--is_metal", "-m", is_flag=True,
              help="Flag to indicate that the structure is metallic. This "
                   "will make the algorithm choose Methfessel-Paxton "
                   "smearing of 0.2 eV.")
@click.option("--hse_calculation", "-H", is_flag=True)
@click.option("--in_custodian", "-C", is_flag=True)
def dimers(site_index, structure_file, distance, is_metal, hse_calculation,
           in_custodian):
    """
    Testing for the workflow scripts.

    """
    from pybat.workflow import site_dimers_workflow

    site_dimers_workflow(structure_file=structure_file,
                         site_index=site_index,
                         distance=distance,
                         is_metal=is_metal,
                         hse_calculation=hse_calculation,
                         in_custodian=in_custodian)
