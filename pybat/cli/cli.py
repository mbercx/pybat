
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
def setup():
    """
    Set up calculations.
    """
    pass

@setup.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".",
              help="Directory in which to set up the calculations for the "
                   "first step in the transition path determination. Note "
                   "that this directory has to contain the structure files")
def transition(directory, initial_):
    """
    Set up a NEB calculation to study a transition of state in the structure.
    """

    from pybat.cli.commands.setup import set_up_transition




@setup.command(context_settings=CONTEXT_SETTINGS)
def migration():
    """
    Set up a NEB calculation to study the migration of elements in the
    structure.
    """
    pass

@setup.command(context_settings=CONTEXT_SETTINGS)
def neb():
    """
    Set up the Nudged Elastic Band calculation based on the output in the
    initial and final directory.

    Returns:

    """

