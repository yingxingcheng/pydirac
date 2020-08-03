#!/usr/bin/env python

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


def get_symbol_and_charge(filename = 'atom.mol'):
    """Get Info from `mol` file

    Args:
        filename (str): filename of `mol` file

    Returns:
        A list of atoms
    """

    with open(filename, 'r') as f:
        context = f.read()

    pattern = r'^\s+(\d+)\.\s+(\d+\.?\d?)\s+'
    re_obj = re.compile(pattern)
    atoms = re.findall(pattern, context, re.MULTILINE)
    print(atoms)
    return atoms


def get_energy(filename, method='CCSD(T)'):
    """Get energy from a CC calculation.

    Give a CC calculation output file specified by `fname` and which type energy
    do you need given by `method`, the options are 'CCSD(T)', 'CCSD', 'MP2' or
    'SCF'

    Args:
        filename (str): the filename of RELCCSD calculation output
        method (str): energy type

    Returns:
        A list of energies
    """

    with open(filename, 'r') as f:
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
            'Did not find energy for {0} in {1}'.format(method, filename))


if __name__ == '__main__':
    pass
