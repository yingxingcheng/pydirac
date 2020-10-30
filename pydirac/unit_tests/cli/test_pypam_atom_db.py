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

from pydirac.cli.pypam_atom_db import get_polarizability


import os
import numpy as np
import pathlib

module_dir = pathlib.Path(__file__).resolve().parent.parent
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))
In_res_dir = os.path.abspath(os.path.join(data_dir, 'In_so_res'))
In_q_so_dir = os.path.abspath(os.path.join(data_dir, 'In_q_so'))
In_q_mrci_dir = os.path.abspath(os.path.join(data_dir, 'In_q_mrci_res'))
Br_q_mrci_dir = os.path.abspath(os.path.join(data_dir, 'Br_q_mrci'))
acv4z_new_dir = os.path.abspath(os.path.join(data_dir, 'acv4z_new'))
wrong_dir = os.path.abspath(os.path.join(data_dir, 'wrong'))
Rn_q_so_dir = os.path.abspath(os.path.join(data_dir, 'Rn_q_so'))
K_mrci_dir = os.path.abspath(os.path.join(data_dir, 'K_mrci'))
B_mrci_dir = os.path.abspath(os.path.join(data_dir, 'B_mrci_wrong'))


def test_get_polarizability():
    print(module_dir)
    print(get_polarizability(In_q_so_dir, calc_dir_patters=['dyall'], deepth=1))
    # get_polarizability(In_res_dir)
    # get_polarizability(Br_q_mrci_dir, deepth=1)
    # get_polarizability(acv4z_new_dir, deepth=1)
    # get_polarizability(wrong_dir, deepth=1)
    print(get_polarizability(Rn_q_so_dir, deepth=1))
    # get_polarizability(K_mrci_dir, deepth=1)
    print(get_polarizability(B_mrci_dir, deepth=0))

# test_get_polarizability()

def calc_ci_average(data, word, nb=2, verbos=False):
    tmp_res = [v[0] for k, v in data.items() if
               word in k and nb <= int(k.split('_')[-1])]
    if len(tmp_res):
        ave = sum(tmp_res)/ len(tmp_res)
        if verbos:
            if nb == 2:
                print('Average of root 1 and root 2 with {0} is: {1:.3f}'.format(word, ave))
            elif nb > 2:
                print('Average of root 1-{2} and root 2 with {0} is: {1:.3f}'.format(word, ave, nb))
            else:
                print('Ground-state value with {0} is: {1:.3f}'.format(word, ave))
        return ave

# def ond_directory():

def test_extract_info():
    print(module_dir)
    # B_res = get_polarizability(B_mrci_dir, deepth=0)
    res = get_polarizability(In_q_so_dir, calc_dir_patters=['dyall'], deepth=1)
    for k,v in res.items():
        # curr_dir, calc_dir, sub_dir
        if len(v) < 1:
            continue
        print(k)
        if k in ['curr_dir', 'calc_dir']:
            for kk, vv in v.items():
                calc_type, data = kk, vv
                print('{0:15s} {1:>50s}'.format('Parameter:', calc_type))

                for w in range(1, 5):
                    ave = calc_ci_average(data,'sym_' + str(w), nb=-1)
                    if ave is not None:
                        print('{0:15s} {1:50.3f}'.format('Final:',ave))

#test_extract_info()
test_get_polarizability()
