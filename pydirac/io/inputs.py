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

import os
import re
import warnings
import subprocess
from monty.json import MSONable, jsanitize, MontyDecoder
from monty.os import cd

from pydirac import SETTINGS, __version__
from pydirac.core.settings import Settings
from pydirac.core.molecule import Molecule
from pydirac.utility.config import get_mol_by_custom_basis, \
    get_mol_by_default_basis
from pydirac.core.periodic_table import Element


class DiracInput(dict, MSONable):
    """
    Class to contain a set of dirac input objects corresponding to a run.
    """

    def __init__(self, inp, mol, optional_files=None, **kwargs):
        """
        Args:
            inp: Inp object
            mol: Mol object
            optional_files: Other input files supplied as a dict of {
            filename: object}. The object should follow standard conventions in
            implementing a as_dict() adn from_dict method.
            **kwargs:
        """
        super().__init__(**kwargs)
        self.update(
            {"inp":inp, "mol":mol}
        )
        if optional_files is not None:
            self.update(optional_files)

    def __str__(self):
        output = []
        for k,v in self.items():
            output.append(k)
            output.append(str(v))
            output.append("")
        return "\n".join(output)

    def as_dict(self) -> dict:
        """
        Returns:
            MSONable dict
        """
        d = {k: v.as_dict() for k, v in self.items()}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d: Dict representation.

        Returns:
            VaspInput
        """
        dec = MontyDecoder()
        sub_d = {"optional_files": {}}
        for k, v in d.items():
            #TODO: here, inp and mol is not filename but suffix
            if k in ["inp", "mol"]:
                sub_d[k.lower()] = dec.process_decoded(v)
            elif k not in ["@module", "@class"]:
                sub_d["optional_files"][k] = dec.process_decoded(v)
        return cls(**sub_d)

    def write_input(self, output_dir=".", make_dir_if_not_present=True):
        """
        Write DIRAC input to a directory.
        Args:
            output_dir: Directory to write to. Defaults to current directory
                (".").
            make_dir_if_not_present: Create the directory if not present.
                Defaults to True.
        Returns:
        """
        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for k, v in self.items():
            if v is not None:
                with open(os.path.join(output_dir, k), 'wt') as f:
                    f.write(v.__str__())

    @staticmethod
    def from_directory(input_dir, optional_files=None):
        """
        Read in a set of DIRAC input from a directory. Note that only the standard
        inp, mol files are read unless otional_filenames is specified.

        Args:
            input_dir: Directory to read DIRAC input from.
            opotional_files: Optional files to read in as well as a dict of {
            filename: Object type}. Object type must have a static method from_file.

        Returns:

        """
        sub_d = {}
        for fname, ftype in [
            ('PYDIRAC.inp', Inp),
            ('PYDIRAC.mol', Mol)
        ]:
            try:
                fullpath = os.path.join(input_dir, fname)
                sub_d[fname.lower()] = ftype.from_file(fullpath)
            except FileNotFoundError:
                sub_d[fname.lower()] = None

        sub_d["optional_files"] = {}
        if optional_files is not None:
            for fname, ftype in optional_files.items():
                sub_d["optional_files"][fname] = ftype.from_file(
                    os.path.join(input_dir, fname)
                )
        return DiracInput(**sub_d)

    def run_dirac(self, run_dir =".",
                  dirac_cmd: list = None,
                  output_file = "dirac.out",
                  err_file = "dirac.err"
                  ):
        """
        Write input files and run DIRAC.
        Args:
            run_dir: Where to write input files and do the run.
            dirac_cmd: Args to be supplied to run DIRAC.
            output_file: File to write output
            err_file: File to write err

        Returns:
        """
        self.write_input(output_dir=run_dir)
        dirac_cmd = dirac_cmd or SETTINGS.get("DIRAC_EXE")
        dirac_cmd = [os.path.expanduser(os.path.expandvars(t)) for t in dirac_cmd]
        if not dirac_cmd:
            raise RuntimeError(
                "You need to supply dirac_cmd or set the DIRAC_EXE in .pydirac.yaml to run DIRAC."
            )
        with cd(run_dir):
            with open(output_file, "w") as f_std, open(
                    err_file, "w", buffering=1
            ) as f_err:
                subprocess.check_call(dirac_cmd, stdout=f_std, stderr=f_err)


class Inp(dict, MSONable):
    """Class to represent DIRAC input file

    """
    _allowed_duplicated_list = ['CIROOTS']
    _top = ['dirac', 'analyze', 'hamiltonian', 'integrals', 'general']

    def __init__(self, params=None):
        super(Inp, self).__init__()
        if params:
            pass
        self.update(params)
        self._parse_input()

    def as_dict(self) -> dict:
        d = dict(self)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        return Inp({k: v for k, v in d.items() if k not in ("@module", "@class")})

    @staticmethod
    def add_sub_node(settinging_obj, dir_node, subdir_node,
                     keyword_node, value_list):
        """Add sub node to parent node
        """
        assert (type(dir_node) is not None)
        assert (type(value_list) == list)

        # whether add value to parent node
        # if yes: value_list has been added to parent node
        # if no: nothing has been done
        is_set = False

        if subdir_node is not None and keyword_node is not None:
            # -------------------------------------------------------
            # Case 1, for example:
            # -------------------------------------------------------
            # **ANALYZE
            # .MULPOP
            # *MULPOP
            # .VECPOP
            # 1..oo
            if len(value_list) > 0:
                settinging_obj[dir_node][subdir_node][keyword_node] = value_list
            else:
                settinging_obj[dir_node][subdir_node][keyword_node] = True
            is_set = True
        elif subdir_node is None and keyword_node is not None:
            # -------------------------------------------------------
            # Case 2, for example:
            # -------------------------------------------------------
            # **DIRAC
            # .TITLE
            # B, DOSSSS, KRCI
            if len(value_list) > 0:
                settinging_obj[dir_node][keyword_node] = value_list
            else:
                settinging_obj[dir_node][keyword_node] = True
            is_set = True
        elif subdir_node is not None and keyword_node is None:
            # if subdir_node is not None and keyword_node is None
            # -------------------------------------------------------
            # Case 3, for example:
            # -------------------------------------------------------
            # **INTEGRALS
            # *READINP
            # .UNCONTRACT
            #
            # when it's processing .UNCONTRACT, it's the case 3
            pass
        else:
            # if subdir_node is None and keyword_node is None
            # -------------------------------------------------------
            # Case 4, for example:
            # -------------------------------------------------------
            # **DIRAC
            #
            # when processing the part in which there is only directory
            pass
        return is_set

    @classmethod
    def from_string(cls, lines):
        setting = {}

        old_dir = None
        old_subdir = None
        old_dotkey = None
        curr_value_list = []

        for line in lines:
            if line.startswith('#'):
                # this is a comment line
                continue

            if line.startswith('**'):
                if old_dir is not None:
                    is_set = Inp.add_sub_node(setting, old_dir, old_subdir,
                                          old_dotkey, curr_value_list)
                    if is_set: curr_value_list = []
                if len(curr_value_list) > 0:
                    Inp.add_sub_node(setting, old_dir, old_subdir, old_dotkey,
                                 curr_value_list)

                # this is a directory
                dir_name = line.lstrip('**').rstrip()
                old_dir = dir_name
                # setting[old_dir] = Settings()
                setting[old_dir] = {}

                ## clear current value for subdir and keyword
                old_subdir = None
                old_dotkey = None
                curr_value_list.clear()

            elif line.startswith('*'):
                # this is a sub-directory
                if old_dir is not None:
                    is_set = Inp.add_sub_node(setting, old_dir, old_subdir,
                                          old_dotkey, curr_value_list)
                    if is_set: curr_value_list = []

                if line.startswith('*END'):
                    continue

                sub_dir_name = line.lstrip('*').rstrip()
                old_subdir = sub_dir_name

                assert (old_dir is not None)
                if sub_dir_name in setting[old_dir]:
                    ## if sub_dir_name in old_dir word list
                    setting[old_dir][sub_dir_name] = {}
                    setting[old_dir][sub_dir_name]['_en'] = True
                else:
                    ## this is a real second level directory
                    setting[old_dir][old_subdir] = {}

                ## clear current value for keyword
                old_dotkey = None
                curr_value_list = []

            elif line.startswith('.'):
                # this is dot keyword
                curr_dotkey = line.lstrip('.').rstrip()
                if old_dir is not None:
                    ## first we probabily have two same dotkey, e.g., CIROOTS
                    is_set = Inp.add_sub_node(setting, old_dir, old_subdir,
                                              old_dotkey, curr_value_list)
                    if is_set:
                        curr_value_list = []

                    if old_subdir is not None and curr_dotkey in setting[old_dir][
                        old_subdir].keys() and \
                            curr_dotkey.upper() in Inp._allowed_duplicated_list:
                        # first we need to clear curr_value_list and add this to the
                        # previous one
                        if is_set:
                            # curr_value_list = setting[old_dir][old_subdir][
                            #     curr_dotkey]
                            # curr_value_list.append('.' + curr_dotkey)

                            ## get id:
                            _tmp_k = []
                            for _k in setting[old_dir][old_subdir].keys():
                                if curr_dotkey == _k:
                                    continue
                                if curr_dotkey in _k:
                                    _tmp_k.append(int(_k.split('_id_')[-1]))
                            if len(_tmp_k):
                                new_id = max(_tmp_k) + 1
                            else:
                                new_id = 1
                            curr_dotkey = curr_dotkey + '_id_' + str(new_id)

                        old_dotkey = curr_dotkey
                        continue


                    # if old_subdir is None and curr_dotkey in setting[old_dir]:
                    #     print('{}, {}, {}'.format(old_dir, old_subdir,
                    #                               curr_dotkey))
                    #     curr_value_list.append('.' + curr_dotkey)
                    #     continue

                old_dotkey = curr_dotkey
                assert (not (old_dir is None and old_subdir is None))

                ## clear value for value list
                curr_value_list = []
            else:
                # this is value of keyword
                curr_value = line.rstrip()
                if len(curr_value):
                    curr_value_list.append(curr_value)
        return Inp.from_dict(setting)

    @classmethod
    def from_file(cls, filename):
        """
        Parse DIRAC input file to restore python Settings object
        """
        with open(filename, 'r') as f:
            lines = f.readlines()
        return Inp.from_string(lines)

    def get_string(self):
        """Transform all contents of ``input`` branch of ``settings``
        into string with blocks, subblocks, keys and values.

        On the highest level alphabetic order of iteration is modified:
        keys occuring in class attribute ``_top`` are printed first.
        See :ref:`dirac-input` for details.
        """
        is_empty = lambda x: isinstance(x, dict) and len(x) == 0

        def parse_key(key, value):
            id_ptn = re.compile('^(.*)_id_\d+$')
            m = id_ptn.match(key)
            if m:
                key = m.group(1)
                if key.upper() not in self._allowed_duplicated_list:
                    raise RuntimeError('Duplicated key does not support!')
            ret = '.' + key.upper() + '\n'
            if not (value is True or is_empty(value)):
                if isinstance(value, list):
                    for i in value:
                        ret += str(i) + '\n'
                else:
                    ret += str(value) + '\n'
            return ret

        def parse_block(block):
            enabler = '_en'
            ret = '**' + block.upper() + '\n'
            s = self[block]
            for k,v in s.items():
                if not isinstance(v, dict) or is_empty(v):
                    ret += parse_key(k, v)
            for k,v in s.items():
                if isinstance(v, dict) and enabler in v:
                    ret += parse_key(k, v[enabler])
            for k,v in s.items():
                if isinstance(v, dict) and len(v) > 0:
                    ret += '*' + k.upper() + '\n'
                    for kk,vv in v.items():
                        if kk != enabler:
                            ret += parse_key(kk, vv)
            return ret

        inp = ''
        for block in self._top:
            if block in self:
                inp += parse_block(block)
        for block in self:
            if block.lower() not in self._top:
                inp += parse_block(block)
        inp += '*END OF INPUT\n'
        return inp

    def __str__(self):
        return self.get_string()

    @property
    def use_wavefunc(self):
        return self._use_wavefunc

    @property
    def has_hamiltonian(self):
        return self._has_hamiltonian

    @property
    def hamiltonian(self):
        return self._hamiltonian

    @property
    def calc_type(self):
        return self._calc_type

    @property
    def calc_method(self):
        return self._calc_method

    @property
    def electric_field(self):
        return self._electric_field

    @property
    def ci_calc_orbit(self):
        if hasattr(self, '_ci_calc_orbit'):
            return self._ci_calc_orbit
        return {'occ':0, 'vir':0}

    def _find_real_key(self, k, keys):
        for _k in keys:
            if k in _k:
                return _k
        else:
            return None

    @property
    def wf_tag(self):
        return self._wf_tag

    def _parse_input(self):
        """Parse input from output file

        Returns:
            input settings
        """

        inp_settings = Settings(self)

        self._wf_tag = None
        for dk in inp_settings.dirac:
            if 'WAVE F' in dk:
                self._use_wavefunc = True
                self._wf_tag = self._find_real_key('WAVE F', self)
                break
        else:
            self._use_wavefunc = False
            self._has_hamiltonian = False
            self._hamiltonian = None
            self._calc_type = None
            self._calc_method = None
            self._electric_field = None
            return

        # check Hamiltonian
        if 'HAMILTONIAN' in inp_settings:
            self._has_hamiltonian =  True
            if 'DOSSSS' in inp_settings.hamiltonian:
                # this is 4-component Hamiltonian
                self._hamiltonian = '4C'
            elif 'X2C' in inp_settings.hamiltonian:
                # this is 2-component Hamiltonian
                self._hamiltonian = '2C'
            elif 'X2Cmmf' in inp_settings.hamiltonian:
                # need to do
                self._hamiltonian = '2C'
            else:
                warnings.warn('We do not know how many components, '
                              'please check inp file')
                self._hamiltonian = '?C'

            if 'NOSPIN' in inp_settings.hamiltonian:
                self._hamiltonian += '-SR'
            else:
                self._hamiltonian += '-SO'
        else:
            self._has_hamiltonian = False

        # check if CC or CI calculation
        if self._has_hamiltonian and 'RELCC' in inp_settings and \
                'RELCCSD' in inp_settings[self._wf_tag]:
            self._calc_method = 'CC'
        elif self._has_hamiltonian and 'KRCICALC' in inp_settings[self._wf_tag] \
                and 'KR CI'in inp_settings[self._wf_tag]:
            self._calc_method = 'CI'
            # here, we can obtain calc_orbit info from input directly
            gas_tag = None
            for k in inp_settings[self._wf_tag]['KRCICALC'].keys():
                if 'GAS' in k:
                    gas_tag = k
                    break
            gas_res = inp_settings[self._wf_tag]['KRCICALC'][gas_tag]
            # the first line is the number of GAS shell
            gas_nb = int(gas_res[0])
            assert gas_nb == len(gas_res[1:])

            nb_e = 0
            t_orb = []
            gas_ptn = re.compile(r'\s*\d+\s+(?P<nb_e>\d+)\s*/\s*(?P<nb_orb>\d+)')
            for l in gas_res[1:]:
                m = gas_ptn.match(l)
                if m:
                    t_orb.append(int(m.group('nb_orb')))
            m = gas_ptn.match(gas_res[-1])
            if m:
                nb_e = int(m.group('nb_e'))
                tot_orb = sum(t_orb)
                occ_orb = nb_e
                vir_orb = tot_orb * 2 - occ_orb
                self._ci_calc_orbit = {'occ': occ_orb, 'vir': vir_orb}
            else:
                self._ci_calc_orbit = {'occ':0, 'vir': 0}

        elif self._has_hamiltonian and 'SCF' in inp_settings[self._wf_tag]:
            if inp_settings[self._wf_tag]['SCF']['_en']:
                self._calc_method = 'SCF'
            else:
                self._calc_method = None
        # TODO: visual module calculation
        else:
            self._calc_method = None

        # check dipole or quadrupole
        if self._has_hamiltonian and 'OPERATOR' in inp_settings.hamiltonian:
            if ' ZZTHETA' in inp_settings.hamiltonian.operator:
                self._calc_type = 'Q'
            elif ' ZDIPLEN' in inp_settings.hamiltonian.operator:
                self._calc_type = 'D'
            else:
                warnings.warn('this is not finite-field calculations')
                self._calc_type = 'non_finite_field'
        else:
            if self._has_hamiltonian:
                warnings.warn('Hamiltonian is invalid, dir')
            else:
                warnings.warn('".OPERATOR" is not in Hamiltonian')
            self._calc_type = 'non_finite_field'

        # check electric field
        if self._calc_type in 'QD':
            # this is finit-field
            pos = inp_settings.hamiltonian.operator.index(' COMFACTOR')
            self._electric_field = inp_settings.hamiltonian.operator[pos+1]
        else:
            warnings.warn('Did not find electric field, '
                          'dir:')
            self._electric_field = 'null'

    def find_case(self, key):
        """Check if this instance contains a key consisting of the same
        letters as *key*, but possibly with different case. If found, return
        such a key. If not, return *key*.

        Args:
            key:

        Returns:

        """
        if not isinstance(key, str):
            return key
        lowkey = key.lower()
        for k in self:
            if k.lower() == lowkey:
                return k
        return key

    def __contains__(self, name):
        """Like regular ``__contains`__``, but ignore the case."""
        return dict.__contains__(self, self.find_case(name))

    def __getitem__(self, name):
        """Like regular ``__getitem__``, but ignore the case."""
        return dict.__getitem__(self, self.find_case(name))

    def __setitem__(self, name, value):
        """Like regular ``__setitem__``, but ignore the case and if the
        value is a dict, convert it to |Settings|.
        """
        dict.__setitem__(self, self.find_case(name), value)

    def __delitem__(self, name):
        """Like regular ``__detitem__``, but ignore the case.
        """
        return dict.__delitem__(self, self.find_case(name))

    def write_file(self, filename):
        with open(filename, "wt") as f:
            f.write(self.__str__())


class Mol(MSONable):
    """Class to represent DIRAC Mol file

    """

    def __init__(
            self,
            molecule:Molecule,
            comment: str = None,
            basis_type: str = None,
            basis_lib: str = 'BASIS',
    ):
        self.molecule = molecule
        self.comment = comment
        self.basis_type = basis_type
        self.basis_lib = basis_lib

    def as_dict(self) -> dict:
        d = {
            'molecule': self.molecule.as_dict(),
        }
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__

        return d

    @classmethod
    def from_file(cls, filename):
        # TODO, more info from mol file including basis sets.
        with open(filename, 'r') as f:
            lines = f.readlines()
        return Mol.from_string(lines)

    @classmethod
    def from_string_TODO(cls, lines):
        """Parse mol from mol file

         The first 3 lines are arbitrary comment lines. You can leave them
         blank or write some note to yourself but they have to be there.
         Line 4 specifies a cartesian GTO basis set (C) and 2 atom types
         (basis set types).

        Returns:
            Mol settings
        """

        #   # start_line = re.compile(r"^Contents of the molecule file\s*")
        #   #dash_line = re.compile(r"^-+")
        #   #empty_line = re.compile(r"^$")
        #   # inp_start_line = re.compile(r"^\*\*DIRAC\s*")
        #   nuclei_nb_line = re.compile(
        #       r"\s*(?P<nuclei>[-+]?(\d+(\.\d*)?|\d*\.\d+))\s+(?P<nb_atoms>\d+)")
        #   basis_line = re.compile(r"^\b(?:LARGE|EXPLICIT)\b\s+BASIS\s+(\S+)\s*")
        #   endline = re.compile(r"^FINISH")
        #   basis_summary_ptn = re.compile('^\s*(?P<gto_type>C)\s+(?P<nb_bs>\d+)\s*')

        #   # nuclei_nb_line =  r"\s*(?P<nuclei>[-+]?(\d+(\.\d*)?|\d*\.\d+))\s+(?P<nb_atoms>\d+)\s*"
        #   # coord_line = r"\w\d+\s+(?P<coord>(?:[-+]?(?:\d+(?:\.\d*)?|\d*\.\d+)\s+)+)"


        #   nuclei_id = 0
        #   nb_atoms = 0
        #   basis_type = None
        #   gto_type = 'C'
        #   nb_basis_set = 0
        #   coordinates = []
        #   atoms = []


        #   # the first three lines are comments
        #   comments = lines[0:3]

        #   m = basis_summary_ptn.match(lines[4])
        #   if m:
        #       nb_basis_set = m.group('nb_bs')
        #       gto_type = m.group('gto_type')
        #       nb_atoms = nb_basis_set

        #   lines = lines[5:]
        #   if nb_atoms > 1:
        #       for i in range(nb_atoms):
        #           pass


        #   for i, line in enumerate(lines):
        #       if endline.match(line):
        #           break
        #       elif nuclei_nb_line.match(line):
        #           m = nuclei_nb_line.match(line)
        #           nuclei_id = int(round(float(m.group('nuclei'))))
        #           nb_atoms = int(m.group('nb_atoms'))
        #           # print(cls.nuclei_id, cls.nb_atoms)
        #       elif basis_line.match(line):
        #           m = basis_line.match(line)
        #           if m:
        #               basis_type = m.group(1)
        #       else:
        #           # we do not interest with these lines
        #           continue

        #   if basis_type is None:
        #       basis_type = lines[4].strip()

        #   mol_settings = {'nuclei_id': nuclei_id, 'nb_atoms': nb_atoms,
        #                   'basis_type': basis_type}

    @classmethod
    def from_string(cls, context):
        """Parse mol from output file

        Returns:
            Mol settings
        """

        start_line = re.compile(r"^Contents of the molecule file\s*")
        dash_line = re.compile(r"^-+")
        empty_line = re.compile(r"^$")
        # inp_start_line = re.compile(r"^\*\*DIRAC\s*")
        nuclei_nb_line = re.compile(
            r"\s*(?P<nuclei>[-+]?(\d+(\.\d*)?|\d*\.\d+))\s+(?P<nb_atoms>\d+)")
        basis_line = re.compile(r"^\b(?:LARGE|EXPLICIT)\b\s+BASIS\s+(\S+)\s*")
        endline = re.compile(r"^FINISH")

        nuclei_id = 0
        nb_atoms = 0
        basis_type = None

        # with open(self.filename, 'r') as f:
        #     context = f.readlines()

        # for i, line in enumerate(context):
        #     if start_line.match(line):
        #         context = context[i + 1:]

        for i, line in enumerate(context):
            if endline.match(line):
                break
            elif dash_line.match(line) or empty_line.match(line):
                continue
            elif nuclei_nb_line.match(line):
                m = nuclei_nb_line.match(line)
                nuclei_id = int(round(float(m.group('nuclei'))))
                nb_atoms = int(m.group('nb_atoms'))
                # print(self.nuclei_id, self.nb_atoms)
            elif basis_line.match(line):
                m = basis_line.match(line)
                if m:
                    basis_type = m.group(1)
            else:
                # we do not interest with these lines
                continue

        if basis_type is None:
            basis_type = context[4].strip()
        basis_type = basis_type
        # TODO: we are doing atomic calculation
        assert nb_atoms == 1
        molecule = Molecule(atoms=[nuclei_id], corrdinates=[[0., 0., 0.,]])
        return cls(basis_type=basis_type, basis_lib='BASIS', molecule=molecule)

    @staticmethod
    def get_string(atom_info: Element, basis_type: str,
                      basis_lib: str = 'BASIS', ) -> None:
        atom_index = atom_info.Z
        atom_type = atom_info.symbol

        if basis_lib not in ['EXPLICIT', 'BASIS']:
            raise TypeError('Basis type should be "BASIS" or "EXPLICIT" '
                            'for builtin basis or custom basis.')

        if basis_lib == 'EXPLICIT':
            with open('basis/{0}.dat'.format(atom_type), 'r') as f:
                basis_info = f.read()
            template = get_mol_by_custom_basis(atom_type, atom_index,
                                               basis_lib, basis_info)
        elif basis_lib == 'BASIS':
            template = get_mol_by_default_basis(atom_type, atom_index,
                                                basis_type)
        else:
            raise TypeError('Basis type should be "BASIS" or "EXPLICIT" '
                            'for builtin basis or custom basis.')
        return template

    def __str__(self):
        basis_type = self.basis_type or 'null'
        return Mol.get_string(self.molecule.atomic_info, basis_type)

    def write_file(self, filename=None):
        if self.molecule.is_atom:
            fname = filename or self.molecule.atomic_info.symbol + '_' +\
                    self.basis_type + '.mol'
            with open(fname, 'w') as f:
                basis_type = self.basis_type or 'null'
                f.write(Mol.get_string(self.molecule.atomic_info.Z, basis_type))
        else:
            raise NotImplementedError('For many-atoms molecule, '
                                      'it has not been implemented yet! ')


if __name__ == '__main__':
    # setting = parse_dirac_input()
    # job = DiracJob(settings=setting)

    # with open('tmp2.inp', 'w') as f:
    #     f.write(job.get_input())

    def main():
        import os

        module_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.abspath(os.path.join(module_dir,'../unit_tests/', 'data'))
        dirac_inp = os.path.join(data_dir, 'tmp.inp')


        setting = Inp.from_file(dirac_inp)
        # setting['hello']['xxkk']= {'dfd':'s'}
        print(setting['WAVE FUNCTIONS']['KRCICALC'].keys())
        print(setting)
        print(dir(setting))
        #print(setting)
        #print(Settings(setting))
        # print(setting.as_dict())
        # job = DiracJob(settings=setting)
        # print(job.get_input())

    def main2():
        molecule = Molecule(['H'], [[0.0,0., 0.]])
        mol = Mol(molecule, basis_type='dyall.acv4z')

        print('*'*80)
        print(mol)
        print('*'*80)

    # main2()
    def main3():
        import os

        module_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.abspath(os.path.join(module_dir,'../unit_tests/', 'data'))
        dirac_inp = os.path.join(data_dir, 'tmp.inp')

        inp = Inp.from_file(dirac_inp)
        molecule = Molecule(['H'], [[0.0,0., 0.]])
        mol = Mol(molecule, basis_type='dyall.acv4z')

        dirac_input = DiracInput(inp=inp, mol=mol)
        dirac_input.write_input(os.path.join(data_dir, 'dirac_input'), make_dir_if_not_present=True)

    main3()

