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
from pydirac.utility.config import get_element_by_id


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
        name = get_element_by_id(i)
        with open('basis/{0}.dat'.format(name),'w') as f:
            nb_block = res_list[i-1].count('functions')
            f.write('LARGE EXPLICIT {0} {1}\n'.format(nb_block, '1 '*nb_block))
            f.write(res_list[i-1])
            f.write('FINISH\n')


if __name__ == '__main__':
    # basis_helper(filename='dirac_basis')
    pass
