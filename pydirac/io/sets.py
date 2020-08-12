# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Pydirac: PYthon tool for DIRAC software.
#  Copyright (C) 2020-2020 The Pydirac Development Team
#
#  This file is part of Pydirac.
#
#  Pydirac is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 3
#  of the License, or (at your option) any later version.
#
#  Pydirac is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/>
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

"""
This module defines the DiracInputSet abstract base class and a concrete
implementation for the parameters developed and tested by the author.
The basic concept behind an input set is to specify a scheme to generate
a consistent set of DIRAC inputs from a structure without further user
intervention. This ensures comparability across runs.
"""


import abc
from pathlib import Path
from copy import deepcopy
import shutil

from monty.json import MSONable, jsanitize
from monty.serialization import loadfn

from pydirac.io.inputs import DiracInput, Mol, Inp
from pydirac.core.molecule import Molecule
from collections import OrderedDict
import re
from monty.json import MSONable
from os.path import join as opj
# from mendeleev import element
from pydirac.core.periodic_table import Element
from pydirac.core.settings import Settings
from pydirac.utility.config import *
from pydirac.io.inputs import Mol
from pydirac.io.basic import input_from_calctype


dir_path = os.path.dirname(os.path.abspath(__file__))


MODULE_DIR = Path(__file__).resolve().parent


class DiracInputSet(MSONable, metaclass=abc.ABCMeta):
    """
    Base class representing a set of Dirac input parameters with a structure
    supplied as init parameters. Typically, you should not inherit from this
    class. Start from DictSet or MPRelaxSet or MITRelaxSet.
    """

    @property
    @abc.abstractmethod
    def inp(self):
        """Inp object"""
        pass

    @property
    @abc.abstractmethod
    def mol(self):
        """Mol object
        Returns:
        """
        pass

    def get_dirac_input(self):
        return DiracInput(inp=self.inp, mol=self.mol)

    def write_input(self, output_dir, make_dir_if_not_present=True):
        """
        Writes a set of DIRAC input to a directory.
        Returns:

        """
        dinput = self.get_dirac_input()
        dinput.write_input(output_dir, make_dir_if_not_present=make_dir_if_not_present)


def _load_yaml_config(fname):
    config = loadfn(str(MODULE_DIR / ("%s.yaml" % fname)))
    if "PARENT" in config:
        parent_config = _load_yaml_config(config["PARENT"])
        for k, v in parent_config.items():
            if k not in config:
                config[k] = v
            elif isinstance(v, dict):
                v_new = config.get(k, {})
                v_new.update(v)
                config[k] = v_new
    return config


class DictSet(DiracInputSet):

    def __init__(
            self,
            molecule,
            config_dict,
            files_to_transfer=None,
            user_inp_settings=None,
            user_mol_settings=None,
            use_structure_charge=False
    ):
        self._molecule = molecule
        self._config_dict = deepcopy(config_dict)
        self.files_to_transfer = files_to_transfer or {}
        self.user_inp_settings = user_inp_settings or {}
        self.user_mol_settings = user_mol_settings or {}
        self.use_molecule_charge = use_structure_charge

    @property
    def inp(self) -> Inp:
        """
        Returns: Inp object
        """
        settings = dict(self._config_dict["inp"])
        inp = Inp(settings)
        return inp

    @property
    def molecule(self) -> Molecule:
        return self._molecule

    @property
    def mol(self) -> Mol:
        """
        Returns: Mol object
        """
        mol = Mol(self.molecule)
        mol_settings = dict(self._config_dict["mol"])
        basis_type = mol_settings.get('basis_type', 'dyall.v2z')
        mol.basis_type = basis_type

        if 'basis_type' in self.user_mol_settings:
            mol.basis_type = self.user_mol_settings['basis_type']
        return mol

    def __str__(self):
        return self.__class__.__name__

    def write_input(self, output_dir, make_dir_if_not_present=True):
        super().write_input(
            output_dir=output_dir,
            make_dir_if_not_present=make_dir_if_not_present
        )
        for k,v in self.files_to_transfer.items():
            shutil.copy2(v, str(Path(output_dir) / k))


class AtomicDHFSet(DictSet):
    CONFIG = _load_yaml_config("AtomicDHFSet")

    def __init__(self, molecule, **kwargs):
        super().__init__(molecule, AtomicDHFSet.CONFIG, **kwargs)
        self.kwargs = kwargs


class AtomicCCSet(AtomicDHFSet):
    CONIFG = _load_yaml_config("AtomicCCSet")

    def __init__(
            self,
            molecule,
            prev_inp=None,
            prev_mol=None,
            **kwargs
    ):
        super(AtomicCCSet, self).__init__(molecule, **kwargs)
        # self._config_dict['inp'].update()
        if isinstance(prev_inp, str):
            prev_inp = Inp.from_file(prev_inp)
        if isinstance(prev_mol, str):
            prev_mol = Mol.from_file(prev_mol)

        self.prev_inp = prev_inp
        self.prev_mol = prev_mol
        self.kwargs = kwargs

    def inp(self):
        parent_inp = super().inp
        inp = (
            Inp(self.prev_inp)
            if self.prev_inp is not None
            else Inp(parent_inp)
        )

        inp.update(
            {

            }
        )

    def mol(self):
        pass


if __name__ == '__main__':
    pass
    # import sys

    # argv = sys.argv[1:]

    # relcc = ('**RELCC', ['.NEL_F1', '3 0 0 0 2 0 0 0', '.ENERGY',
    #                      '.PRINT', '1', '*CCENER',
    #                    '.MAXIT', '60', '.NTOL', '10'])
    # # inactivate = (closed_elec - elec_in_gas)//2
    # # gas1 = '10 12 / 6' # for elements from third row
    # # gas2 = '{0} {1} / 3'.format(open_elec+10, open_elec+12) # for p open-shell orbit
    # # gas3 = '{0} {1} / {2}'.format(open_elec+12, open_elec+12, nb_virtual)
    # # krci = ('*KRCICALC',['.CI PROGRAM', 'LUCIAREL','.INACTIVE',inactivate,
    # #                      '.GAS SHELLS',3, gas1,gas2,gas3,'.CIROOTS',ciroot1,'.CIROOTS',ciroot2,'.MAX CI','120', '.MXCIVE','60','.ANALYZ','.RSTRCI','rstr','.CHECKP'])
    # for f in argv:
    #     inp = Inpobj.from_file(f)
    #     inp.add_keywords(relcc)
    #     # inp.add_keywords(krci)
    #     inp.write_to_file('tmp_B.inp')

    def main():
        import os
        molecule = Molecule(['H'], [[0.0,0., 0.]])
        dhf = AtomicDHFSet(molecule=molecule, user_mol_settings={"basis_type":'dyall.acv4z'})
        print(dhf.inp)

        module_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.abspath(os.path.join(module_dir,'../unit_tests/', 'data'))
        dhf.write_input(output_dir=os.path.join(data_dir, 'dirac_input'), make_dir_if_not_present=True)

    main()
