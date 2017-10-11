
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
def migration():
    """
    Set up a NEB calculation to study the migration of elements in the
    structure.
    """
    pass


@main.group(context_settings=CONTEXT_SETTINGS)
def setup():
    """
    Set up calculations.
    """
    pass


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
def transition(directory, is_migration):
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
                   "for the initial and final state. ")
@click.option("--nimages", "-n", default=8,
              help="Number of images.")
@click.option("--is_migration", "-m", is_flag=True,
              help="Flag to designate the transition as a migration. "
                   "Activating this flag means that a static calculation will "
                   "be set up to determine the charge density of the host "
                   "structure, i.e. without the migrating ion. This charge "
                   "density can then be used to find a good initial guess "
                   "for the migration pathway.")
def neb(directory, nimages, is_migration):
    """
    Set up the Nudged Elastic Band calculation based on the output in the
    initial and final directory.

    Returns:

    """
    from pybat.cli.commands.setup import set_up_NEB

    set_up_NEB(directory=directory,
               nimages=nimages,
               is_migration=is_migration)
