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
from pydirac.core.settings import Settings
from pydirac.io.sets import DiracJob

module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))


def test_write():
    #set single calculation parameters (single point, TZ2P/PW91)
    sett = Settings()
    sett.input.dirac.title = 'B, DOSSSS, KRCI'
    sett.input.dirac.analyze = True
    sett.input.dirac['WAVE FUNCTIONS'] = True
    sett.input.analyze.mulpop._en = True
    sett.input.analyze.mulpop.vecpop = '1..oo'
    sett.input.hamiltonian.dossss = True
    sett.input.hamiltonian.operator._en = True
    sett.input.hamiltonian.operator = [' ZDIPLEN', ' COMFACTOR', ' 0.01']
    sett.input.integrals.readinp.uncontract = True
    sett.input.general.pcmout = True
    sett.input['WAVE FUNCTIONS'].scf._en = True
    sett.input['WAVE FUNCTIONS']['KR CI'] = True
    sett.input['WAVE FUNCTIONS'].resolve = True
    sett.input['WAVE FUNCTIONS'].scf['CLOSED SHELL'] = 4
    sett.input['WAVE FUNCTIONS'].scf['OPEN SHELL'] = [1, '1/6']
    sett.input['WAVE FUNCTIONS'].scf['EVCCNV'] = '1.0D-9 5.0D-8'
    sett.input['WAVE FUNCTIONS'].scf['MAXITR'] = 90
    sett.input['WAVE FUNCTIONS'].krcicalc['CI PROGRAM'] = 'LUCIAREL'
    sett.input['WAVE FUNCTIONS'].krcicalc.inactive = 1
    sett.input['WAVE FUNCTIONS'].krcicalc['GAS SHELLS'] = \
    [3, '0 2 / 1', '1 3 / 3', '3 3 / 10']
    sett.input['WAVE FUNCTIONS'].krcicalc['MAX CI'] = 60
    sett.input['WAVE FUNCTIONS'].krcicalc['NOOCCN'] = True
    sett.input['WAVE FUNCTIONS'].krcicalc['DIPMOM'] = True
    sett.input['WAVE FUNCTIONS'].krcicalc['RSTRCI'] = 0
    sett.input['WAVE FUNCTIONS'].krcicalc['CIROOTS'] = '3  3'

    sett.runscript.pam.mol = '/Users/yxcheng/PhD/plams_tutorial/B.mol'

    job = DiracJob(settings=sett)

    fout = os.path.join(data_dir, 'tmp.inp')
    with open(fout, 'w') as f:
        f.write(job.get_input())
