
import os

from pymatgen import Structure
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.analysis.transition_state import NEBAnalysis


def get_structure(directory):
    """
    Construct a .json file with the structure and magnetic moment from the
    output of a VASP calculation, i.e. the CONTCAR and OUTCAR file.

    Args:
        directory:

    Returns:

    """
    directory = os.path.abspath(directory)
    structure = Structure.from_file(os.path.join(directory, "CONTCAR"))
    out = Outcar(os.path.join(directory, "OUTCAR"))

    magmom = [site["tot"] for site in out.magnetization]

    structure.add_site_property("magmom", magmom)
    structure.to("json", "structure.json")


def get_barrier(directory):
    """
    Plot the migration barrier of a transition in a directory.
    Args:
        directory:

    Returns:

    """
    neb = NEBAnalysis.from_dir(directory)
    neb.get_plot().show()
