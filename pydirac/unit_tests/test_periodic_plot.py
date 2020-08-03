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

from pydirac.utility.periodic_plot import get_plot, DataType


def test_table():
    # get_plot(DataType.qp_SO_acv3z)
    # get_plot(DataType.qp_SR)
    # get_plot(DataType.dp_SR)
    # get_plot(DataType.dp_SO)
     get_plot(DataType.qp_SO)
    # get_plot(DataType.qp_SO_acv3z, cmap='spring')
    # get_plot(DataType.qp_SO_acv3z, cmap='plasma')
    # get_plot(DataType.qp_SO_acv3z, cmap='magma')
    # get_plot(DataType.qp_SO_acv3z, cmap='inferno')
