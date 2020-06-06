import os
from pydirac.input.basic import *

module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))

dirac_inp = os.path.join(data_dir, 'He.inp')


def test_dossss_scf():
    get_dossss_scf_inp('He', filename_out=dirac_inp, is_dipole=True)
    get_dossss_scf_inp('He', filename_out=dirac_inp, is_dipole=False)


def test_x2c_scf():
    get_x2c_scf_inp('He', filename_out=dirac_inp, is_dipole=True, is_nospin=True)
    get_x2c_scf_inp('He', filename_out=dirac_inp, is_dipole=True, is_nospin=False)
    get_x2c_scf_inp('He', filename_out=dirac_inp, is_dipole=False, is_nospin=True)
    get_x2c_scf_inp('He', filename_out=dirac_inp, is_dipole=False, is_nospin=False)


def test_x2c_relccsd():
    # get_x2c_relccsd_inp('He', filename_out=dirac_inp, is_dipole=False, is_nospin=False)
    # get_x2c_relccsd_inp('He', filename_out=dirac_inp, is_dipole=True, is_nospin=False)
    get_x2c_relccsd_inp('He', filename_out=dirac_inp, is_dipole=True, is_nospin=True)
    # get_x2c_relccsd_inp('He', filename_out=dirac_inp, is_dipole=True, is_nospin=False)
