#!/usr/bin/env python

import os
import sys
import subprocess
import re
import warnings
from pydirac.utility.read import parse_dirac_input
from pydirac.core.settings import Settings
from pydirac.core.atom import Atom


def get_symbol_and_charge(fname='atom.mol'):
    with open(fname, 'r') as f:
        context = f.read()

    pattern = r'^\s+(\d+)\.\s+(\d+\.?\d?)\s+'
    re_obj = re.compile(pattern)
    atoms = re.findall(pattern, context, re.MULTILINE)
    print(atoms)
    return atoms


# @ Total CCSD(T) energy :                     -15.138231245354447
def get_energy(fname, method='CCSD(T)'):
    """
    method: CCSD(T), MP2, SCF
    """
    with open(fname, 'r') as f:
        context = f.read()

    if '(' and ')' in method:
        # add \ for (
        method = method.replace('(', '\(').replace(')', '\)')
    pattern = r'^@.*{method} energy[\s:]*(-?[\d\.]+)'.format(
        **{'method': method})

    # re_obj = re.compile(pattern)
    energy = re.findall(pattern, context, re.MULTILINE)

    if energy:
        return energy[-1]
    else:
        raise RuntimeError(
            'Did not find energy for {0} in {1}'.format(method, fname))


class Outupt():
    """
    Class to parse DIRAC output file
    """

    def __init__(self, filename):
        self.filename = filename
        self.input = None
        self.mol = None

        if filename:
            with open(filename, 'r') as f:
                self.context = f.readlines
        else:
            raise ValueError("filename is None!")

    def parse(self):
        """
        Parse output file to create several objects or restore all info to
        a dict
        Returns
        -------

        # generate Input object and Mol object?
        # extract info:
            Input:
                1) SO or SR,
                2) zff,
                3) dipole or quadrupole,
                4) CC or CI,
                5) X2C or DOSSSS
            Mol:
                1) molecule
                2) basis sets: default (e.g., dyall.acv4z) or explicit

        """
        # check input file
        inp_settings = self.parse_input()
        self.input = Input(inp_settings)
        mol_settings = self.parse_input()
        self.mol = Mol(mol_settings)
        orbit_settings = self.parse_orbit()
        results_settings = self.parse_results()

    def parse_input(self):
        """
        Parse input from output file
        Parameters
        ----------
        filename: output filename

        Returns
        -------
        str

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

    def parse_mol(self):
        """
        Parse mol from output file
        Parameters
        ----------

        Returns
        -------

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
        """
        Parse orbits from output file
        Parameters
        ----------
        filename

        Returns
        -------

        """
        atom = Atom.from_file(self.filename)
        print('The number of electrons at range ({0}, {1}): {2} '.format(
            e_min, e_max, atom.electron_count(-20, 25)))

    def parse_results(self):
        """
        Parse calculation results from output file
        Parameters
        ----------
        filename

        Returns
        -------

        """
        if not hasattr(self, 'inp_settings'):
            self.parse_input()

        if not hasattr(self, 'mol_settings'):
            self.parse_mol()

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
            self.electric_field = float(self.inp_settings.input.hamiltonian.operator[pos+1])
        else:
            warnings.warn('Did not find electric field!')

        if self.calc_method in ['CI', 'CC']:
            if self.calc_method == 'CC':
                self.energy_settings = self._parse_coupled_cluster()
            else:
                self.energy_settings = self._parse_configuration_interaction()

    def get_task_type(self):
        """
        Get a task type from an output file
        Returns
        -------

        """
        if not (hasattr(self, 'hamiltonian') and hasattr(self, 'calc_type') and
                hasattr(self, 'calc_method') and hasattr(self, 'basis_type')):
            self.parse_results()

        self.task_type = '-'.join([self.calc_type,self.hamiltonian,
                       self.calc_method]) + '@' + self.basis_type.strip()

        return self.task_type


    def _parse_coupled_cluster(self, has_T=True):
        """
        Deal with RELCCSD calculations

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

        Returns
        -------

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
        """
        Deal with KRCI calculations
        Returns
        -------

        """
        pass

    def _parse_hartree_fock(self):
        """
        Deal with HF calculations
        Returns
        -------

        """
        pass

    def _parse_visual_module(self):
        """
        Deal with visual module calculations
        Returns
        -------

        """
        pass



class Input():
    """
    Class to represent DIRAC input file
    """

    def __init__(self, settings):
        self.settings = settings


class Mol():
    """
    Class to represent DIRAC Mol file
    """

    def __init__(self, settings):
        self.setting = settings


if __name__ == '__main__':
    pass
