# coding: utf8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os

import numpy as np

from itertools import chain
from monty.serialization import loadfn
from pymatgen.core import Structure, PeriodicSite
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.sets import DictSet

"""
Package that described the various calculation sets used for the setup scripts.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "pre-alpha"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Mar 2019"

MODULE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "set_configs/")


def _load_yaml_config(fname):
    config = loadfn(os.path.join(MODULE_DIR, "%s.yaml" % fname))
    return config


class BulkRelaxSet(DictSet):
    """
    VASP input set for the bulk relaxation.

    """
    CONFIG = _load_yaml_config("bulkRelaxSet")

    def __init__(self, structure, **kwargs):
        super(BulkRelaxSet, self).__init__(
            structure, BulkRelaxSet.CONFIG, **kwargs)
        self.kwargs = kwargs


class BulkSCFSet(DictSet):
    """
    VASP input set for the bulk SCF calculation.

    """
    CONFIG = _load_yaml_config("bulkSCFSet")

    def __init__(self, structure, **kwargs):
        super(BulkSCFSet, self).__init__(
            structure, BulkSCFSet.CONFIG, **kwargs)
        self.kwargs = kwargs


class PybatNEBSet(BulkRelaxSet):
    """
    Class for writing NEB inputs, based on the settings of BulkRelaxSet.

    This class was largely copied from the MITNEBSet in pymatgen.io.vasp.sets,
    But then using the defaults specified in BulkRelaxSet. We've also added a new
    method for visualizing the transition.

    Args:
        \\*\\*kwargs: Other kwargs supported by :class:`DictSet`.

    """

    def __init__(self, structures, **kwargs):
        if len(structures) < 3:
            raise ValueError("You need at least 3 structures for an NEB.")
        kwargs["sort_structure"] = False
        super(PybatNEBSet, self).__init__(structures[0], **kwargs)
        self._structures = self._process_structures(structures)

        if "EDIFF" not in self._config_dict["INCAR"]:
            self._config_dict["INCAR"]["EDIFF"] = self._config_dict[
                "INCAR"].pop("EDIFF_PER_ATOM")

        # NEB specific defaults
        defaults = {'IMAGES': len(structures) - 2, 'IBRION': 1, 'ISYM': 0,
                    'LCHARG': False, 'LCLIMB': True}
        self._config_dict["INCAR"].update(defaults)

    @property
    def poscar(self):
        return Poscar(self.structures[0])

    @property
    def poscars(self):
        return [Poscar(s) for s in self.structures]

    @property
    def structures(self):
        return self._structures

    def _process_structures(self, structures):
        """
        Remove any atom jumps across the cell

        """
        input_structures = structures
        structures = [input_structures[0]]
        for s in input_structures[1:]:
            prev = structures[-1]
            for i in range(len(s)):
                t = np.round(prev[i].frac_coords - s[i].frac_coords)
                if np.any(np.abs(t) > 0.5):
                    s.translate_sites([i], t, to_unit_cell=False)
            structures.append(s)
        return structures

    def write_input(self, output_dir, make_dir_if_not_present=True,
                    write_cif=False, write_path_cif=False,
                    write_endpoint_inputs=False):
        """
        NEB inputs has a special directory structure where inputs are in 00,
        01, 02, ....

        Args:
            output_dir (str): Directory to output the VASP input files
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            write_cif (bool): If true, writes a cif along with each POSCAR.
            write_path_cif (bool): If true, writes a cif for each image.
            write_endpoint_inputs (bool): If true, writes input files for
                running endpoint calculations.

        """

        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        self.incar.write_file(os.path.join(output_dir, 'INCAR'))
        self.kpoints.write_file(os.path.join(output_dir, 'KPOINTS'))
        self.potcar.write_file(os.path.join(output_dir, 'POTCAR'))

        for i, p in enumerate(self.poscars):
            d = os.path.join(output_dir, str(i).zfill(2))
            if not os.path.exists(d):
                os.makedirs(d)
            p.write_file(os.path.join(d, 'POSCAR'))
            if write_cif:
                p.structure.to(filename=os.path.join(d, '{}.cif'.format(i)))
        if write_endpoint_inputs:
            end_point_param = BulkRelaxSet(
                self.structures[0],
                user_incar_settings=self.user_incar_settings)

            for image in ['00', str(len(self.structures) - 1).zfill(2)]:
                end_point_param.incar.write_file(
                    os.path.join(output_dir, image, 'INCAR')
                )
                end_point_param.kpoints.write_file(
                    os.path.join(output_dir, image, 'KPOINTS')
                )
                end_point_param.potcar.write_file(
                    os.path.join(output_dir, image, 'POTCAR')
                )
        if write_path_cif:
            sites = set()
            lattice = self.structures[0].lattice
            for site in chain(*(s.sites for s in self.structures)):
                sites.add(PeriodicSite(site.species_and_occu, site.frac_coords,
                                       lattice))
            nebpath = Structure.from_sites(sorted(sites))
            nebpath.to(filename=os.path.join(output_dir, 'path.cif'))

    def visualize_transition(self, filename="transition.cif"):
        """
        Write a file to show the transition by simply adding all images.

        """

        transition_structure = self.structures[0].copy()
        for structure in self.structures[1:]:
            for site in structure:
                transition_structure.append(site.specie, site.frac_coords)

        transition_structure.to("cif", filename)
