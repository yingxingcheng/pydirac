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

from pydirac.io.inputs import Inp
# from pydirac.io.sets import DiracJob
import os
from pydirac.core.settings import Settings

module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))
dirac_inp = os.path.join(data_dir, 'tmp.inp')


def test_read_func():
    input = Inp.from_file(dirac_inp)
    # input['ANALYZE']['MULPOP']['VECPOP'] = 22
    # input['ANALYZE']['MULPOP']['vecpop_id_01'] = 22
    print(input['WAVE FUNCTIONS']['KRCICALC']['CIROOTS'])
    print(input)
    # print(Settings(input))
    #print(DiracJob(input).get_input())
    #print(Settings(input))
    #input.write_file('tmp2.inp')
    # print(Settings(input))
    # print(setting.as_dict())
    # job = DiracJob(settings=setting)
    # print(job)

test_read_func()