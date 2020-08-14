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
from pydirac.io.inputs import Inp

module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))


def test_write():
    #set single calculation parameters (single point, TZ2P/PW91)
    sett = Settings()
    sett.dirac.title = 'B, DOSSSS, KRCI'
    sett.dirac.analyze = True
    wf_tag = 'WAVE FUNCTIONS'
    sett.set_nested(('DIRAC', wf_tag), True)
    sett.analyze.mulpop._en = True
    sett.analyze.mulpop.vecpop = '1..oo'
    sett.hamiltonian.dossss = False
    sett.hamiltonian.x2c = True
    sett.hamiltonian.nospin = True
    sett.hamiltonian.operator = [' ZDIPLEN', ' COMFACTOR', ' 0.01']
    sett.integrals.readinp.uncontract = True
    sett.general.pcmout = True

    wave_func = Settings()
    wave_func['KR CI'] = True
    wave_func.resolve = True

    scf = Settings()
    scf._en = True
    scf['CLOSED SHELL'] = 4
    scf['OPEN SHELL'] = [1, '1/6']
    scf.evccnv = '1.0D-9 5.0D-8'
    scf.maxitr = 90

    wave_func.scf = scf
    sett[wf_tag] = wave_func

    krci = Settings()
    krci['CI PROGRAM'] = 'LUCIAREL'
    krci['INACTIVE'] = 1
    krci['GAS SHELLS'] = [3, '0 2 / 1', '1 3 / 3', '3 3 / 10']
    krci['MAX CI'] = 60
    krci['NOOCCN'] = True
    krci['DIPMOM'] = True
    krci['RSTRCI'] = 0
    krci['CIROOTS_id_0'] = '3  3'
    krci['CIROOTS_id_1'] = '4  3'
    sett.set_nested(('WAVE FUNCTIONS','KRCICALC'), krci)

    job = Inp(sett.as_dict())
    print(job)

test_write()
