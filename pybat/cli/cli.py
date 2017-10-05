
import click

"""
Command line interface for the pybat package.

"""

@click.group()
def main():
    """
    Tools for performing calculations for studying batteries.
    """
    pass

@main.group()
def setup():
    """
    Set up calculations for VASP.
    """
    pass

@setup.command()
def migrate():
    """
    Set up a NEB calculation to study the migration of elements in the
    structure.
    """

