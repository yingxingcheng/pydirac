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

import re
import warnings
from monty.json import MSONable, jsanitize
from pydirac.io.inputs import Inp, Mol, parse_dirac_input
from pydirac.core.settings import Settings
from pydirac.core.molecule import Molecule

class Output(MSONable):
    """Class to parse DIRAC output file

    Attributes:
        filename (str): output filename of DIRAC. Normally, a good output file
            contains `inp` file and `mol` file, thus one can parse them here
            together.

        input (Settings): an `input` settings
        mol (Settings): an `mol` settings
        inp_settings (Settings): an input settings
        mol_settings (Settings): a mol settings

    """

    def __init__(self, filename):
        self.filename = filename
        self.input = None
        self.mol = None
        self.basis_type = None
        self.calc_method = None
        self.calc_type = None
        self.electric_field = None
        self.energy_settings = None
        self.hamiltonian = None
        self.has_hamiltonian = False
        self.inp_settings = None
        self.input = None
        self.mol = None
        self.mol_settings = None
        self.nb_atoms = 0
        self.nuclei_id = None
        self.atom_info = None

        self.parse()

    def parse(self):
        """Parse the whole output file

        Parse the whole output file including `inp` file and `mol` file and
        `orbit` information and calculation `results` in the `output` file.

        Note:
            Input:
                1) SO or SR,
                2) zff,
                3) dipole or quadrupole,
                4) CC or CI,
                5) X2C or DOSSSS
            Mol:
                1) molecule
                2) basis sets: default (e.g., dyall.acv4z) or explicit

        Returns:
            None
        """

        self.parse_input()
        self.parse_mol()
        self.parse_orbit()

        self.parse_results()

        if not self.inp_settings:
            self.input = Inp(self.inp_settings)
        if not self.mol_settings:
            self.mol = Mol(self.mol_settings)


    def parse_input(self):
        """Parse input from output file

        Returns:
            input settings
        """

        start_line = re.compile(r"^Contents of the input file\s*")
        dash_line = re.compile(r"^-+")
        empty_line = re.compile(r"^$")
        inp_start_line = re.compile(r"^\*\*DIRAC\s*")
        endline = re.compile(r"^\s*\*END OF")

        with open(self.filename, 'r') as f:
            context = f.readlines()

        for i, line in enumerate(context):
            if start_line.match(line):
                context = context[i + 1:]

        inp_str = []
        is_start = False
        for i, line in enumerate(context):
            if is_start and not endline.match(line):
                inp_str.append(line)
            elif inp_start_line.match(line):
                is_start = True
                inp_str.append(line)
            elif dash_line.match(line) or empty_line.match(line):
                continue
            else:
                inp_str.append(line)
                break

        self.inp_settings = parse_dirac_input(''.join(inp_str))

        # check Hamiltonian
        if 'HAMILTONIAN' in self.inp_settings.input:
            self.has_hamiltonian =  True
            if 'DOSSSS' in self.inp_settings.input.hamiltonian:
                # this is 4-component Hamiltonian
                self.hamiltonian = '4C'
            elif 'X2C' in self.inp_settings.input.hamiltonian:
                # this is 2-component Hamiltonian
                self.hamiltonian = '2C'
            elif 'X2Cmmf' in self.inp_settings.input.hamiltonian:
                # need to do
                self.hamiltonian = '2C'
            else:
                warnings.warn('We do not know how many components, '
                              'please check inp file')
                self.hamiltonian = '?C'

            if 'NOSPIN' in self.inp_settings.input.hamiltonian:
                self.hamiltonian += ':SR'
            else:
                self.hamiltonian += ':SO'
        else:
            self.has_hamiltonian = False

        # check if CC or CI calculation
        if self.has_hamiltonian and 'RELCC' in self.inp_settings.input and \
                'RELCCSD' in self.inp_settings.input['WAVE FUNCTIONS']:
            self.calc_method = 'CC'
        elif self.has_hamiltonian and 'KRCICALC' in self.inp_settings.input['WAVE FUNCTIONS'] \
                and 'KR CI'in self.inp_settings.input['WAVE FUNCTIONS']:
            self.calc_method = 'CI'
        else:
            self.calc_method = None
            pass

        # check dipole or quadrupole
        if self.has_hamiltonian and 'OPERATOR' in self.inp_settings.input.hamiltonian:
            if ' ZZTHETA' in self.inp_settings.input.hamiltonian.operator:
                self.calc_type = 'Q'
            elif ' ZDIPLEN' in self.inp_settings.input.hamiltonian.operator:
                self.calc_type = 'D'
            else:
                warnings.warn('this is not finite-field calculations')
                self.calc_type = 'non_finite_field'
        else:
            warnings.warn('this is not finite-field calculations')
            self.calc_type = 'non_finite_field'

        # check electric field
        if self.calc_type in 'QD':
            # this is finit-field
            pos = self.inp_settings.input.hamiltonian.operator.index(' COMFACTOR')
            self.electric_field = self.inp_settings.input.hamiltonian.operator[pos+1]
        else:
            warnings.warn('Did not find electric field!')
            self.electric_field = 'NULL'

    def parse_mol(self):
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

        self.nuclei_id = 0
        self.nb_atoms = 0
        self.basis_type = None

        with open(self.filename, 'r') as f:
            context = f.readlines()

        for i, line in enumerate(context):
            if start_line.match(line):
                context = context[i + 1:]

        for i, line in enumerate(context):
            if endline.match(line):
                break
            elif dash_line.match(line) or empty_line.match(line):
                continue
            elif nuclei_nb_line.match(line):
                m = nuclei_nb_line.match(line)
                self.nuclei_id = int(round(float(m.group('nuclei'))))
                self.nb_atoms = int(m.group('nb_atoms'))
                # print(self.nuclei_id, self.nb_atoms)
            elif basis_line.match(line):
                m = basis_line.match(line)
                self.basis_type = m.group(1)
            else:
                # we do not interest with these lines
                continue

        self.mol_settings = Settings({'nuclei_id': self.nuclei_id,
                                      'nb_atoms': self.nb_atoms,
                                      'basis_type': self.basis_type})

    def parse_orbit(self, e_min=-20.0, e_max=25.0):
        """Parse orbits from output file

        Args:
            e_min (float): energy minimum
            e_max (float): energy maximum

        Returns:
            None
        """

        self.atom_info = Molecule.from_file(self.filename)
        if self.atom_info:
            print('The number of electrons at range ({0}, {1}): {2} '.format(
                e_min, e_max, self.atom_info.electron_count(-20, 25)))

        # TODO: get atom information and restore to a Settings object


    def parse_results(self):
        """Parse calculation results from output file

        Returns:
            None
        """

        if not self.inp_settings:
            self.parse_input()

        if not self.mol_settings:
            self.parse_mol()


        if self.calc_method in ['CI', 'CC']:
            if self.calc_method == 'CC':
                self.energy_settings = self._parse_coupled_cluster()
            else:
                self.energy_settings = self._parse_configuration_interaction()
        else:
            self.energy_settings = Settings()

        self.task_type = '-'.join([self.calc_type,self.hamiltonian,
                                   self.calc_method]) + '@' + self.basis_type.strip() + \
                         '@' + self.electric_field.strip()

    def get_task_type(self):
        """Get a task type from an output file

        Returns:
        """
        return self.task_type


    def _parse_coupled_cluster(self, has_T=True):
        """Deal with RELCCSD calculations

        Notes:
            An example with RELCCSD calculations

             Overview of calculated energies
            @ SCF energy :                             -5880.437843209574567
            @ MP2 correlation energy :                    -0.486037237589513
            @ CCSD correlation energy :                   -0.451464950662248
            @ 4th order triples correction :              -0.016640269726526
            @ 5th order triples (T) correction :           0.000318657059527
            @ 5th order triples -T  correction :           0.000362317692980
            @ Total MP2 energy :                       -5880.923880447164265
            @ Total CCSD energy :                      -5880.889308160236396
            @ Total CCSD+T  energy :                   -5880.905948429963246
            @ Total CCSD(T) energy :                   -5880.905629772903922
            @ Total CCSD-T  energy :                   -5880.905586112270612

        Returns:
            None
        """
        start_pattern = re.compile(r'^\s*Overview of calculated energies')
        scf_pattern = re.compile(r'^@ SCF energy :\s+(?P<energy>[-+]?(\d+(\.\d*)?|\d*\.\d+))')
        mp2_pattern = re.compile(r'^@ Total MP2 energy :\s+(?P<energy>[-+]?(\d+(\.\d*)?|\d*\.\d+))')
        ccsd_pattern = re.compile(r'^@ Total CCSD energy :\s+(?P<energy>[-+]?(\d+(\.\d*)?|\d*\.\d+))')
        ccsd_p_T_pattern = re.compile(r'^@ Total CCSD\(T\) energy :\s+(?P<energy>[-+]?(\d+(\.\d*)?|\d*\.\d+))')
        end_line = re.compile(r"^$")

        with open(self.filename, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if start_pattern.match(line):
                    lines = lines[i:]
                    break

        energies = {}
        if has_T:
            patterns = [scf_pattern, mp2_pattern, ccsd_pattern, ccsd_p_T_pattern]
            keys = ['scf_e', 'mp2_e', 'ccsd_e', 'ccsd_p_T_e']
        else:
            patterns = [scf_pattern, mp2_pattern, ccsd_pattern]
            keys = ['scf_e', 'mp2_e', 'ccsd_e']

        for i, line in enumerate(lines):
            if end_line.match(line):
                break

            for pattern, k in zip(patterns, keys):
                m = pattern.match(line)
                if m:
                    energies[k] = float(m.group('energy'))
                    break

        e_settings = Settings(energies)
        return e_settings


    def _parse_configuration_interaction(self):
        """Deal with KRCI calculations

        Notes:
            CI example:
                Final CI energies  =  -5880.8481422675641  -5880.8380283826846  -5880.8380141759062

                root    1 ...... converged!

                root    2 ...... converged!

                root    3 ...... converged!

        Returns:
            None
        """
        #start_line = re.compile(r'Final CI energies\s+=\s+(?P<energy>([-+]?(\d+(\.\d*)?|\d*\.\d+)\s+)+)\s+')
        energy_line = re.compile(r'Final CI energies\s+=\s+(?P<energy>(?:[-+]?(?:\d+(?:\.\d*)?|\d*\.\d+)\s+)+)\s+')
        conv_line = re.compile(r'(?:\s*root\s+\d+\s+\.+\s+(?:converged|unconverged)!\s*)+')

        with open(self.filename, 'r') as fin:
            context = fin.read()

        res_list = energy_line.findall(context, re.MULTILINE|re.DOTALL)
        conv_list = conv_line.findall(context, re.MULTILINE|re.DOTALL)

        energies = []
        for e_str in res_list:
            e_roots = [float(e) for e in e_str.split()]
            energies.append(e_roots)

        convergeds = []
        for syms in conv_list:
            converged = []
            for sym_res in syms.split('\n'):
                if not len(sym_res.strip()):
                    continue
                if 'unconverged' in sym_res:
                    converged.append(False)
                elif 'converged' in sym_res:
                    converged.append(True)
                else:
                    print(sym_res)
                    raise RuntimeError('There is no converged info')
            convergeds.append(converged)

        energy_settings = Settings({'ci_e': energies, 'ci_converged':convergeds})
        return energy_settings


    def _parse_hartree_fock(self):
        """Deal with HF calculations

        Returns:
            None
        """

        pass

    def _parse_visual_module(self):
        """Deal with visual module calculations

        Returns:
            None
        """

        pass

    def as_dict(self) -> dict:
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "basis_type": self.basis_type,
             "calc_method": self.calc_method,
             "electric_field": self.electric_field,
             "energy_settings": self.energy_settings.as_dict(),
             "filename": self.filename,
             'hamilton': self.hamiltonian,
             #'inp_settings': self.input.as_dict(),
             #'mol_settings': self.mol.as_dict(),
             'inp_settings': self.inp_settings.as_dict(),
             'mol_settings': self.mol_settings.as_dict(),
             'nb_atoms': self.nb_atoms,
             'nuclei_id': self.nuclei_id,
             'task_type': self.task_type,
             'atom_info': self.atom_info.as_dict()
             }
        return jsanitize(d, strict=True)

