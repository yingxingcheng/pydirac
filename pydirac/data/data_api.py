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

# from mendeleev import element
import numpy as np
import csv
from pydirac.core.periodic import Element


module_dir = os.path.dirname(os.path.abspath(__file__))


class DataType:
    dp_SR = 0
    dp_SO = 1
    qp_SR = 2
    qp_SO = 3
    qp_SO_acv3z = 4


def get_data(data_type: int) -> list:
    if data_type == DataType.dp_SO:
        fname = "dp_SO.csv"
    elif data_type == DataType.dp_SR:
        fname = "dp_SR.csv"
    elif data_type == DataType.qp_SR:
        fname = "qp_SR.csv"
    elif data_type == DataType.qp_SO:
        fname = "qp_SO.csv"
    elif data_type == DataType.qp_SO_acv3z:
        fname = "qp_SO_acv3z.csv"
    else:
        raise TypeError("No data type {0}".format(data_type))

    # data = np.loadtxt(fname=fname, delimiter=',', dtype=str)
    atomic_info = [np.NAN] * 118
    fname = os.path.join(module_dir, fname)
    with open(fname, "r") as f:
        reader = csv.reader(f)
        data = [(Element(row[0]).atomic_number - 1, row[1]) for row in reader]

    for atomic_nb, v in data:
        atomic_info[atomic_nb] = float(v)
    return atomic_info


if __name__ == "__main__":
    res = get_data(DataType.dp_SR)
    print(res)
