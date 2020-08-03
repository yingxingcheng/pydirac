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

from mendeleev import get_table
from mendeleev.plotting import periodic_plot
from pydirac.data.data_api import get_data, DataType


def get_plot(data_type: int, cmap='viridis', *args, **kwargs):
    if data_type == DataType.dp_SO:
        attr = 'dipole_polarizability_SO'
        title = 'Dipole polarizability with spin-orbital relativistic effects'
    elif data_type == DataType.dp_SR:
        attr = 'dipole_polarizability_SR'
        title = 'Dipole polarizability with scalar relativistic effects'
    elif data_type == DataType.qp_SR:
        attr = 'quadrupole_polarizability_SR'
        title = 'Quadrupole polarizability with scalar relativistic effects'
    elif data_type == DataType.qp_SO:
        attr = 'quadrupole_polarizability_SO'
        title = 'Quadrupole polarizability with spin-orbital relativistic effects'
    elif data_type == DataType.qp_SO_acv3z:
        attr = 'quadrupole_polarizability_SO_acv3z'
        title = 'Quadrupole polarizability with spin-orbital ' \
                'relativistic effects using d-aug-dyall-acv3z basis'
    else:
        raise TypeError('No data type {0}'.format(data_type))
    data = get_data(data_type)
    ptable = get_table('elements')
    ptable[attr] =  data
    periodic_plot(ptable, attribute=attr, colorby='attribute',
                  title=title, cmap=cmap, *args, **kwargs)

