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
from pydirac.io.basic import *

module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))

dirac_inp = os.path.join(data_dir, 'He.inp')


def test_dossss_scf():
    get_dossss_scf_inp('He', filename_out=dirac_inp, is_dipole=True)
    get_dossss_scf_inp('He', filename_out=dirac_inp, is_dipole=False)


def test_dossss_relccsd():
    get_dossss_relccsd_inp('He', filename_out=dirac_inp, is_dipole=True)
    get_dossss_relccsd_inp('He', filename_out=dirac_inp, is_dipole=False)


def test_x2c_scf():
    get_x2c_scf_inp('He', filename_out=dirac_inp, is_dipole=True, is_nospin=True)
    get_x2c_scf_inp('He', filename_out=dirac_inp, is_dipole=True, is_nospin=False)
    get_x2c_scf_inp('He', filename_out=dirac_inp, is_dipole=False, is_nospin=True)
    get_x2c_scf_inp('He', filename_out=dirac_inp, is_dipole=False, is_nospin=False)


def test_x2c_relccsd():
    #get_x2c_relccsd_inp('He', filename_out=dirac_inp, is_dipole=False, is_nospin=False)
    #get_x2c_relccsd_inp('He', filename_out=dirac_inp, is_dipole=True, is_nospin=False)
    get_x2c_relccsd_inp('He', filename_out=dirac_inp, is_dipole=True, is_nospin=True)
    get_x2c_relccsd_inp('He', filename_out=dirac_inp, is_dipole=False, is_nospin=True)
