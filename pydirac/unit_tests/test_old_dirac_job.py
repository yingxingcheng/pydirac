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

from pydirac.io.inputs import Inp


module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))


def test_old_dirac_job():
    pass
    # scratch_dir = os.path.join(data_dir, 'scratch_dir')
    # OldDiracJob.get_all_input('He', 'dyall.v2z', scratch_dir = scratch_dir,
    #                           calc_type=JobType.d_DOSSSS_SCF, suffix='d_4s_scf')
    # OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
    #                           calc_type=JobType.q_DOSSSS_SCF, suffix = 'q_4s_scf')
    # OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
    #                           calc_type=JobType.d_DOSSSS_RELCCSD, suffix='d_4s_cc')
    # OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
    #                           calc_type=JobType.q_DOSSSS_RELCCSD, suffix='q_4s_cc')
    # OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
    #                           calc_type=JobType.d_X2C_NOSPIN_SCF,suffix='d_2c_nospin_scf')
    # OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
    #                           calc_type=JobType.q_X2C_NOSPIN_SCF, suffix='q_2c_nospin_scf')
    # OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
    #                           calc_type=JobType.d_X2C_NOSPIN_RELCCSD,suffix='d_2c_nospin_cc')
    # OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
    #                           calc_type=JobType.q_X2C_NOSPIN_RELCCSD, suffix='q_2c_nospin_cc')
