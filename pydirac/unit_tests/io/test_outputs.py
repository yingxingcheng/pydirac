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
import pathlib
from pydirac.io.outputs import Output

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

output_wrong_file = os.path.abspath(os.path.join(data_dir, 'parse_output_wrong.out'))


def test_output():
       o = Output(output_wrong_file)
       print(o.as_dict())
