
import os

from pymatgen.analysis.transition_state import NEBAnalysis

"""
Utility command for the pybat package.

"""

def show_path(directory, filename):
    """
    Show the final migration for a NEB calculation.

    Returns:

    """
    # TODO This is quite inefficient, since the NEBAnalysis script also parses the OUTCAR files. However, this was the fastest solution implementation wise.
    neb = NEBAnalysis.from_dir(directory)

    transition_structure = neb.structures[0].copy()
    for structure in neb.structures[1:]:
        for site in structure:
            transition_structure.append(site.specie, site.frac_coords)

    transition_structure.to("cif", filename + ".cif")

def plot_barrier(directory):
    """
    Plot the migration barrier of a transition in a directory.
    Args:
        directory:

    Returns:

    """
    neb = NEBAnalysis.from_dir(directory)
    neb.get_plot().show()