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

from pydirac.core.molecule import Molecule
from pydirac.core.molecular_orbitals import AtomicOrbital
import os
from pathlib import Path
from monty.os import cd


module_dir = Path(__file__).parent.parent.resolve()
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))


def test_get_orbital_info():
    """
    test get orbital info
    """
    # fname = os.path.join(data_dir, 'Kr_DHF.out')
    with cd(os.path.join(data_dir,'Li')):
        for f in ['Li_D-CC-SR.out', 'Li_D-CC-SO.out']:
            mol = Molecule.from_file(filename=f)
            for ao in mol.aos:
                print(ao)
            print(mol.nao(-10, 10))
            print(mol.nb_closed_ao(-10, 10))
            print(mol.nb_open_ao(-10, 10))
            print(mol.nb_virtual_ao(-10, 10))

