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

from pydirac.core.periodic_table import Element
from pydirac.io.inputs import Mol

MODULE_PATH = os.path.dirname(os.path.abspath(__file__))


def basis_helper(filename = 'ANO-RCC'):
    """
    Custom basis set generate script.
    Parameters
    ----------
    filename: file name

    Returns
    -------

    """
    with open(filename, 'r') as f:
        textlines = f.readlines()

    textlines = textlines[43:-1]
    res_list = []
    element_mess = []
    for line in textlines:
        first_word = line.strip().split()[0]
        if first_word == '!' and 'functions' not in line:
            if len(element_mess)>3:
                res_list.append(''.join(element_mess))
                element_mess = []
            element_mess.append(line)
            continue

        elif first_word == 'a':
            element_mess.append('!' + line)
            continue

        elif first_word == 'H':
            nb_exp = int(line.strip().split()[1])
            tmp_str = 'f{0:4d}{1:5d}\n'.format(nb_exp,0)
            element_mess.append(tmp_str)
            continue
        else:
            element_mess.append(line)

    res_list.append(''.join(element_mess))

    if not os.path.exists('basis'):
        os.makedirs('basis')

    assert(len(res_list) == 96)
    for i in range(1,97):
        symbol = Element(i).symbol
        with open('basis/{0}.dat'.format(symbol),'w') as f:
            nb_block = res_list[i-1].count('functions')
            f.write('LARGE EXPLICIT {0} {1}\n'.format(nb_block, '1 '*nb_block))
            f.write(res_list[i-1])
            f.write('FINISH\n')


def get_explicit_basis_default_basis(input_strings, fname_out = None):
    """
    Get explicit basis set from DIRAC default basis library, which means that we
    copy basis info from library to a new file which use 'EXPLICIT' keyword to
    specify basis set.

    In this way, we can add more diffuse functions based on the default library
    by using even-term methods.

    At present, this function supports basis sets as followings:
        1) dyall.acv3z
        2) dyall.acv4z
        3) dyall.cv3z
        4) dyall.cv4z

    Args:
        filename: normal file in which basis sets specified by default, e.g.,
                  dyall.acv4z and dyall.acv3z

    Returns:

    """
    bs_helper = DyallBasisHelper()

    # with open(fname_in, 'r') as f:
    #     lines = f.read()

    lines = input_strings.strip().split('\n')
    Z = int(float(lines[4].split()[0]))
    e_symbol = Element(Z).symbol
    basis_type = lines[6].split()[-1].split('.')[-1]

    bs = bs_helper.get_basis_set(basis_type, e_symbol)

    new_strings = []
    # lines[2] = 'c-' + lines[2]
    lines[2] = 'c-dyall.' + basis_type
    new_strings += lines[:6]
    nb_sym = len(bs.sym_dict)
    new_strings.append('LARGE EXPLICIT {} {}'.format(nb_sym, ' '.join(['1']*nb_sym)))
    new_strings.append(bs.to_string())
    new_strings.append('FINISH')

    if fname_out is None:
        fname_out = e_symbol + '_c-dyall.' + basis_type + '.mol'

    with open(fname_out, 'w') as f:
        f.write('\n'.join(new_strings))


def get_custom_basis_from_ele(ele_type, basis_type, fname_out=None):
    e = Element(ele_type)
    input_strings = Mol.get_mol_by_default_basis(e.symbol, e.Z, basis_type)
    get_explicit_basis_default_basis(input_strings, fname_out=fname_out)



class DyallBasisHelper(object):

    def __init__(self):
        pass

    def get_basis_set(self, basis_type, ele_type):
        if basis_type not in ['acv3z', 'acv4z', 'cv3z', 'cv4z']:
            raise RuntimeError('basis_type must be one of "acv3z", "acv4z", '
                               '"cv3z", and "cv4z"')

        fname = os.path.join(MODULE_PATH, 'dyall.{}'.format(basis_type))
        with open(fname, 'r') as f:
            lines = f.readlines()

        comment = []
        i_start = 0
        for i, line in enumerate(lines):
            match = re.match(r'^\$', line)
            if match:
                comment.append(line.rstrip())
                continue
            else:
                i_start = i
                break

        if i_start == 0:
            raise RuntimeError('Something wrong in the basis set file')

        i_start -= 1
        e = Element(ele_type)

        lines = lines[i_start:]
        for j, line in enumerate(lines):
            match = re.match(r'^\$ {}\s*\n?$'.format(e.symbol), line)
            # if line.startswith('$ {}'.format(e.symbol)):
            if match:
                j_start = j
                break
        else:
            raise RuntimeError("This basis set does not support element "
                               "type {}".format(e.symbol))

        lines = lines[j_start:]
        tmp_strings = []
        for line in lines:
            match = re.match(r'^\s*\n?$', line)
            # if len(line.strip()) == 0:
            if match:
                break
            else:
                tmp_strings.append(line.rstrip())

        dyall_basis = DyallBasis.from_string('\n'.join(tmp_strings))
        return dyall_basis


class DyallBasis(object):

    def __init__(self, ele_Z, sym_dict):
        self.e = Element(ele_Z)
        self.symbol = self.e.symbol
        self.Z = self.e.Z
        self.sym_dict = sym_dict

    @classmethod
    def from_string(cls, bs_text):
        bs_info = bs_text.strip().split('\n')
        # the first line is a comment that include element symbol info
        symbol = bs_info[0].split()[-1]
        Z = int(bs_info[1].split()[-1])
        assert Element(symbol).Z == Z
        sym_dict = {}
        curr_sym = None
        for line in bs_info[2:]:
            if line.startswith('$'):
                sym = line.split()[1]
                if sym not in sym_dict:
                    sym_dict[sym] = []
                curr_sym = sym
            elif len(line.strip().split()) > 1:
                continue
            else:
                sym_dict[curr_sym].append(line.rstrip())
        return DyallBasis(Z, sym_dict)


    def to_string(self):
        tmp_strings = []
        # tmp_strings.append('$ {}'.format(self.symbol))
        # tmp_strings.append('a {}'.format(self.Z))
        for sym in 'spdfghi':
            if sym in self.sym_dict:
                tmp_strings.append('$ {} functions'.format(sym))
                nb_exp = len(self.sym_dict[sym])
                tmp_strings.append('f{0:>4}    0    0'.format(nb_exp))
                for exp in self.sym_dict[sym]:
                    tmp_strings.append('    {0:10.8E}'.format(float(exp)))
        return '\n'.join(tmp_strings)

    def __str__(self):
        return self.to_string()


if __name__ == '__main__':
    # basis_helper(filename='dirac_basis')
    bs_helper = DyallBasisHelper()
    bs = bs_helper.get_basis_set('acv4z', 'Kr')
    print(bs)
    pass
