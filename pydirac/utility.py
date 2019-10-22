#!/usr/bin/env python

import os
import sys
import subprocess
import re


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
        method = method.replace('(', '\(').replace(')', '\)')
    pattern = r'^@.*{method} energy[\s:]*(-?[\d\.]+)'.format(
        **{'method': method})

    # re_obj = re.compile(pattern)
    energy = re.findall(pattern, context, re.MULTILINE)
    print(energy)
    return energy


def execute(command, stdout_file_name='', accepted_errors=[]):
    """
    Runs the command.

    Raises:
        - AcceptedError
        - SubprocessError
    """

    if sys.version_info[0] == 2:
        process = subprocess.Popen(command,
                                   Stdin=Subprocess.Pipe,
                                   Stdout=Subprocess.Pipe,
                                   Stderr=Subprocess.Pipe)
    Else:
        Process = Subprocess.Popen(Command,
                                   Stdin=Subprocess.Pipe,
                                   Stdout=Subprocess.Pipe,
                                   Stderr=Subprocess.PIPE,
                                   encoding='utf8')
    stdout, stderr = process.communicate()
    return stdout


def run(binary_dir='', inp_file, mol_file,
        args='', accepted_errors=[], print_args=False):

    # launch_script = os.path.normpath(os.path.join(binary_dir, 'pam-dirac'))
    # if not os.path.exists(launch_script):
    #     sys.stderr.write('ERROR: launch script %s not found\n' % launch_script)
    #     sys.stderr.write(
    #         '       have you set the correct --binary-dir (or -b)?\n')
    #     sys.stderr.write('       try also --help\n')
    #     sys.exit(-1)
    #
    if sys.platform == "win32":
        dirac_exe = 'dirac.x.exe'
    else:
        dirac_exe = 'dirac.x'

    # launcher = 'python "%s" --dirac=%s --noarch --nobackup %s' % (
    #    launch_script, os.path.join(self.binary_dir, dirac_exe), args)

    launcher = 'pam-dirac --noarch {args}'.format(**{'args': args})
    command = launcher + ' --inp=%s --mol=%s' % (inp_file, mol_file)
    execute(command=command, accepted_errors=accepted_errors)

    # for inp in inp_files:
    #     inp_no_suffix = os.path.splitext(inp)[0]
    #     for mol in mol_files:
    #         mol_no_suffix = os.path.splitext(mol)[0]
    #         out = '%s_%s.out' % (inp_no_suffix, mol_no_suffix)
    #


if __name__ == '__main__':
    # get_symbol_and_charge(context)
    # get_energy(output)

    print(execute('ls'))
