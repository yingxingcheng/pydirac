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

import sys
# from mendeleev import element, Element
from pydirac.core.periodic_table import Element
from pydirac.utility.config import *

dir_path = os.path.dirname(os.path.abspath(__file__))
__all__ = ['get_dossss_scf_inp', 'get_x2c_relccsd_inp', 'get_x2c_scf_inp',
           'input_from_calctype', 'get_dossss_relccsd_inp']


def input_from_calctype(atom_info:Element, calc_type:str,
                        filename_out:str=None):
    """
    Get basis input for DIRAC
    Parameters
    ----------
    atom_info
    calc_type
    filename_out

    Returns
    -------

    """
    # create input file of DIRAC
    elec_conf = get_dirac_shell_str(atom_info.atomic_number)
    func_name = 'get_' + calc_type + '_str'
    context = getattr(sys.modules['pydirac.utility.config'], func_name)\
        (elec_conf, atom_info.symbol)
    inp_fname = atom_info.symbol + '_' + calc_type + '.inp'

    filename_out = filename_out or inp_fname
    with open(filename_out, 'w') as f:
        f.write(context)


# calc_funcs = ['d_DOSSS_SCF', 'd_X2C_SCF', 'd_X2C_NOSPIN_SCF',
#               'q_DOSSSS_SCF', 'q_X2C_NOSPIN_SCF',
#               'd_X2C_NOSPIN_RELCCSD']

def get_dossss_scf_inp(atom_type, filename_out=None, is_dipole=True):
    """
    Get DOSSSS SCF input file
    Parameters
    ----------
    atom_type (int or str): aotm type
    filename_out (str or None): output filename
    is_dipole (bool): True for dipole and False for quadrupole

    Returns
    -------
    None
    """
    calc_str = 'DOSSSS_SCF'
    if is_dipole:
        calc_str = 'd_' + calc_str
    else:
        calc_str = 'q_' + calc_str
    input_from_calctype(Element(atom_type), filename_out=filename_out,
                        calc_type=calc_str)


def get_dossss_relccsd_inp(atom_type, filename_out=None, is_dipole=True):
    """
    Get DOSSSS RELCCSD input file
    Parameters
    ----------
    atom_type (int or str): aotm type
    filename_out (str or None): output filename
    is_dipole (bool): True for dipole and False for quadrupole

    Returns
    -------
    None
    """
    calc_str = 'DOSSSS_RELCCSD'
    if is_dipole:
        calc_str = 'd_' + calc_str
    else:
        calc_str = 'q_' + calc_str
    input_from_calctype(Element(atom_type), filename_out=filename_out,
                        calc_type=calc_str)

def get_x2c_scf_inp(atom_type, filename_out=None, is_dipole=True,
                    is_nospin=False):
    """
    Get X2C SCF calculation input.
    Parameters
    ----------
    atom_type (int or str): aotm type
    filename_out (str or None): output filename
    is_dipole (bool): True for dipole and False for quadrupole

    Returns
    -------
    None

    """
    calc_str = 'X2C'
    if is_dipole:
        calc_str = 'd_' + calc_str
    else:
        calc_str = 'q_' + calc_str

    if is_nospin:
        calc_str += '_NOSPIN_SCF'
    else:
        calc_str += '_SCF'

    input_from_calctype(Element(atom_type), filename_out=filename_out,
                        calc_type=calc_str)


def get_x2c_relccsd_inp(atom_type, filename_out=None, is_dipole=True,
                        is_nospin=False):
    """
    Get X2C RELCCSD input
    Parameters
    ----------
    atom_type (int or str): aotm type
    filename_out (str or None): output filename
    is_dipole (bool): True for dipole and False for quadrupole

    Returns
    -------
    None
    """
    calc_str = 'X2C'
    if is_dipole:
        calc_str = 'd_' + calc_str
    else:
        calc_str = 'q_' + calc_str

    if is_nospin:
        calc_str += '_NOSPIN'

    calc_str += '_RELCCSD'

    input_from_calctype(Element(atom_type), filename_out=filename_out,
                        calc_type=calc_str)


if __name__ == '__main__':
    input_from_calctype(Element('He'), 'd_DOSSSS_SCF')
