#!/usr/bin/env python
import os

from pydirac.io.jobs import OldDiracJob, JobType


module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))


def test_old_dirac_job():
    scratch_dir = os.path.join(data_dir, 'scratch_dir')
    OldDiracJob.get_all_input('He', 'dyall.v2z', scratch_dir = scratch_dir,
                              calc_type=JobType.d_DOSSSS_SCF, suffix='d_4s_scf')
    OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
                              calc_type=JobType.q_DOSSSS_SCF, suffix = 'q_4s_scf')
    OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
                              calc_type=JobType.d_DOSSSS_RELCCSD, suffix='d_4s_cc')
    OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
                              calc_type=JobType.q_DOSSSS_RELCCSD, suffix='q_4s_cc')
    OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
                              calc_type=JobType.d_X2C_NOSPIN_SCF,suffix='d_2c_nospin_scf')
    OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
                              calc_type=JobType.q_X2C_NOSPIN_SCF, suffix='q_2c_nospin_scf')
    OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
                              calc_type=JobType.d_X2C_NOSPIN_RELCCSD,suffix='d_2c_nospin_cc')
    OldDiracJob.get_all_input('He', 'dyall.acv2z', scratch_dir = scratch_dir,
                              calc_type=JobType.q_X2C_NOSPIN_RELCCSD, suffix='q_2c_nospin_cc')
