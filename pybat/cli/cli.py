
import click

"""
Command line interface for the pybat package.

"""

# This is used to make '-h' a shorter way to access the CLI help
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    """
    CLI tools for performing calculations for studying batteries.
    """
    pass


@main.group(context_settings=CONTEXT_SETTINGS)
def define():
    """
    Set up calculations.

    """
    pass


@define.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--provide_coords", "-C", is_flag=True)
@click.option("--write_cif", "-w", is_flag=True)
def migration(structure_file, provide_coords, write_cif):
    """
    Define a migration of an ion in a structure.

    """
    from pybat.cli.commands.define import define_migration

    define_migration(structure_file=structure_file,
                     provide_coords=provide_coords,
                     write_cif=write_cif)


@define.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--dimer_indices", "-i", default=(0,0))
@click.option("--distance", "-d", default=float(0))
@click.option("--remove_cations", "-R", is_flag=True)
def dimer(structure_file, dimer_indices, distance, remove_cations):
    """
    Define the formation of a dimer in a structure.

    """
    from pybat.cli.commands.define import define_dimer

    define_dimer(structure_file=structure_file,
                 dimer_indices=dimer_indices,
                 distance=distance,
                 remove_cations=remove_cations)


@main.group(context_settings=CONTEXT_SETTINGS)
def setup():
    """
    Set up calculations.
    """
    pass


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--calculation_dir", "-d", default="relax",
              help="The directory in which to set up the calculation.")
@click.option("--hse", "-H", is_flag=True)
def relax(structure_file, calculation_dir, hse):
    """
    Set up a geometry optimization for a structure.
    """
    from pybat.cli.commands.setup import set_up_relaxation

    set_up_relaxation(structure_file=structure_file,
                      calculation_dir=calculation_dir,
                      hse_calculation=hse)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".",
              help="Directory in which to set up the calculations for the "
                   "first step in the transition path determination. Note "
                   "that this directory has to contain the structure files "
                   "for the initial and final state. ")
@click.option("--is_migration", "-m", is_flag=True,
              help="Flag to designate the transition as a migration. "
                   "Activating this flag means that a static calculation will "
                   "be set up to determine the charge density of the host "
                   "structure, i.e. without the migrating ion. This charge "
                   "density can then be used to find a good initial guess "
                   "for the migration pathway.")
@click.option("--hse", "-H", is_flag=True)
def transition(directory, is_migration, hse):
    """
    Set up a the geometry optimizations for the initial and final state of a
    transition.
    """

    from pybat.cli.commands.setup import find_transition_structures
    from pybat.cli.commands.setup import set_up_transition

    (initial_structure,
     final_structure) = find_transition_structures(directory)

    set_up_transition(directory=directory,
                      initial_structure=initial_structure,
                      final_structure=final_structure,
                      is_migration=is_migration)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".",
              help="Directory in which to set up the calculations for the "
                   "first step in the transition path determination. Note "
                   "that this directory has to contain the structure files "
                   "for the initial and final state.")
@click.option("--nimages", "-n", default=8,
              help="Number of images.")
@click.option("--is_migration", "-m", is_flag=True,
              help="Flag to designate the transition as a migration. "
                   "Activating this flag means that a static calculation will "
                   "be set up to determine the charge density of the host "
                   "structure, i.e. without the migrating ion. This charge "
                   "density can then be used to find a good initial guess "
                   "for the migration pathway.")
@click.option("--hse", "-H", is_flag=True)
def neb(directory, nimages, is_migration, hse):
    """
    Set up the Nudged Elastic Band calculation based on the output in the
    initial and final directory.

    Returns:

    """
    from pybat.cli.commands.setup import set_up_NEB

    set_up_NEB(directory=directory,
               nimages=nimages,
               is_migration=is_migration,
               hse_calculation=hse)


@main.group(context_settings=CONTEXT_SETTINGS)
def util():
    """
    Utility command for the pybat package.

    """
    pass

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
@click.option("--format", "-F", default="cif")
def conv(structure_file, format):
    """
    Convert a structure into the conventional unit cell.

    """
    from pybat.cli.commands.util import conventional_structure

    conventional_structure(structure_file=structure_file,
                           fmt=format)

@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument("supercell", nargs=1)
@click.argument("structure_file", nargs=1)
@click.option("--format", "-F", default="cif")
def super(supercell, structure_file, format):
    """
    Convert a structure to a supercell.

    """
    from pybat.cli.commands.util import make_supercell

    make_supercell(structure_file=structure_file,
                   supercell=supercell,
                   fmt=format)


@main.group(context_settings=CONTEXT_SETTINGS)
def get():
    """
    Obtain data from output files.

    """
    pass


@get.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".")
def structure(directory):
    """
    Obtain the structure with its magnetic configuration.

    """
    from pybat.cli.commands.get import get_structure

    get_structure(directory=directory)


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

@main.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file")
def test(structure_file):
    """
    Scripts or methods in the testing phase!
    """
    from pybat.core import test_script

    test_script(structure_file)
