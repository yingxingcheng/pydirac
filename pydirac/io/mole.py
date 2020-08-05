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

#from mendeleev import Element, element
from pydirac.utility.config import get_mol_by_custom_basis, \
    get_mol_by_default_basis
from pydirac.core.periodic_table import Element


def get_mole_file(atom_info:Element, basis_type:str, filename_out:str,
                  basis_choice : str = 'BASIS', ) -> None:
    atom_index = atom_info.atomic_number
    atom_type = atom_info.symbol

    if basis_choice not in ['EXPLICIT', 'BASIS']:
        raise TypeError('Basis type should be "BASIS" or "EXPLICIT" '
                        'for builtin basis or custom basis.')


    if basis_choice == 'EXPLICIT':
        with open('basis/{0}.dat'.format(atom_type), 'r') as f:
            basis_info = f.read()
        template = get_mol_by_custom_basis(atom_type, atom_index,
                                           basis_choice, basis_info)
    elif basis_choice == 'BASIS':
        template = get_mol_by_default_basis(atom_type, atom_index, basis_type)
    else:
        raise TypeError('Basis type should be "BASIS" or "EXPLICIT" '
                        'for builtin basis or custom basis.')


    fname = filename_out or atom_type + '_' + basis_type + '.mol'
    with open(fname, 'w') as f:
        f.write(template)

