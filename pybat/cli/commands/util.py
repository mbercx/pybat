
import os

from pymatgen.core import Structure
from pymatgen.analysis.transition_state import NEBAnalysis
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

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


def conventional_structure(structure_file, fmt="cif"):
    """

    Args:
        structure_file:

    Returns:

    """
    structure = Structure.from_file(structure_file)
    spg = SpacegroupAnalyzer(structure)

    conv_structure_file = structure_file.split(".")[0] + "_conv" + ".json"
    spg.get_conventional_standard_structure().to(fmt, conv_structure_file)

