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

from monty.json import MSONable
from pydirac.core.periodic_table import Element


class Molecule(MSONable):

    def __init__(self, atoms=None, corrdinates=None):

        self.atoms = atoms or []
        if len(self.atoms) == 1:
            self.is_atom = True
        else:
            self.is_atom = False

        self.coordinates = corrdinates or []
        assert len(self.atoms) == len(self.coordinates)

        # TODO: we need to check if this is a atomic calculation
        self.atomic_info = 'null'
        if self.is_atom:
            self.atomic_info = Element(self.atoms[0])


    def as_dict(self) -> dict:
        d = {
            'atoms': self.atoms,
            'coordinates' : self.coordinates,
            'is_atom': self.is_atom,
            'atomic_info': self.atomic_info.as_dict()
        }
        return d


if __name__ == '__main__':
    pass

