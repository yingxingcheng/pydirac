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

import os

from pydirac.core.molecule import Molecule
from pydirac.io.outputs import Output
from pydirac.io.sets import AtomicCISet

module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))


def test_write_file():
    """
    test write_file
    """
    fname = os.path.join(data_dir, 'Cu_DHF.out')
    out = Output(fname)
    atom = out.parse_orbit()
    print(atom.as_dict())
    print(atom.nb_closed_elec(), atom.nb_open_elec())
    print(atom.as_dict())

    #fout = os.path.join(data_dir, 'PYDIRAC.inp')
    ci_set = AtomicCISet.from_prev_dhf_calc(fname)
    print(ci_set.inp)
    print(ci_set.mol)

test_write_file()
