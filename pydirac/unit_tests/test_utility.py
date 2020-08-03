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

from pydirac.analysis.utility import *
from pydirac.io.outputs import Output
import os
import glob

module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))
# dirac_inp = os.path.join(data_dir, 'tmp.inp')
In_res_dir = os.path.abspath(os.path.join(data_dir, 'In_so_res'))
In_q_so_dir = os.path.abspath(os.path.join(data_dir, 'In_q_so'))
In_q_mrci_dir = os.path.abspath(os.path.join(data_dir, 'In_q_mrci_res'))
In_d_mrci_dir = os.path.abspath(os.path.join(data_dir, 'In_mrci'))


def test_get_energy():
    current = os.getcwd()
    os.chdir(In_q_so_dir)
    for calc_dir in glob.glob('*'):
        if os.path.isdir(calc_dir):
            os.chdir(os.path.join(In_q_so_dir, calc_dir))
            if 'JOB_DONE' in glob.glob('*'):
                outfiles = glob.glob('*.out')
                if len(outfiles) > 1 :
                    raise RuntimeError('There are two output file in current '
                                       'directory: {0}'.format(calc_dir))
                else:
                    outfile = outfiles[0]
                    energy = get_energy(outfile, method = 'CCSD(T)')
                    print('{0} CCSD(T): {1}'.format(calc_dir, energy))
            os.chdir('..')
    os.chdir(current)


def test_output_object():

    current = os.getcwd()
    os.chdir(In_q_so_dir)
    for calc_dir in glob.glob('*'):
        if os.path.isdir(calc_dir):
            os.chdir(os.path.join(In_q_so_dir, calc_dir))
            if 'JOB_DONE' in glob.glob('*'):
                outfiles = glob.glob('*.out')
                if len(outfiles) > 1 :
                    raise RuntimeError('There are two output file in current '
                                       'directory: {0}'.format(calc_dir))
                else:
                    outfile = outfiles[0]
                    # energy = get_energy(outfile, method = 'CCSD(T)')
                    output_obj = Output(outfile)
                    settings = output_obj.parse_input()
                    #print(output_obj.inp_settings)
                    #settings = output_obj.parse_mol()
                    #settings = output_obj.parse_orbit()
                    #output_obj.parse_results()
                    #print(output_obj.energy_settings)
                    print(output_obj.get_task_type())
            os.chdir('..')
    os.chdir(current)

def test_output_object_CI():

    current = os.getcwd()
    os.chdir(In_d_mrci_dir)
    for calc_dir in glob.glob('dyall*'):
        if os.path.isdir(calc_dir):
            os.chdir(os.path.join(In_d_mrci_dir, calc_dir))
            if 'JOB_DONE' in glob.glob('*'):
                outfiles = glob.glob('*.out')
                if len(outfiles) > 1 :
                    raise RuntimeError('There are two output file in current '
                                       'directory: {0}'.format(calc_dir))
                else:
                    outfile = outfiles[0]
                    # energy = get_energy(outfile, method = 'CCSD(T)')
                    output_obj = Output(outfile)
                    #settings = output_obj.parse_input()
                    #settings = output_obj.parse_mol()
                    #settings = output_obj.parse_orbit()
                    #output_obj.parse_results()
                    #print(output_obj.energy_settings)
                    #print(output_obj.get_task_type())

                    res = output_obj.parse_results()
                    print(output_obj.energy_settings)

                    res_dict = output_obj.as_dict()
                    print(res_dict)


            os.chdir('..')
    os.chdir(current)
