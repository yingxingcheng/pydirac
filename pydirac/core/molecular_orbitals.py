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

from monty.json import MSONable, jsanitize


class AtomicOrbital(MSONable):
    """Orbit object to restore basis information about given orbit.
    """

    def __init__(self, orbit_sym, orbit_type,orbit_energy, orbit_degenerate, orbit_frac):
        """
        :param orbit_sym: orbit symmetry
        :param orbit_type: orbital type, e.g., open-shell or closed-shell
        :param orbit_energy: orbital energy (a.u.)
        :param orbit_degenerate: orbital degenerate
        :param orbit_frac: orbital fraction
        """
        self.orbit_sym = orbit_sym
        self.orbit_type = orbit_type
        self.orbit_energy = float(orbit_energy)
        self.orbit_frac = float(orbit_frac)
        self.orbit_degenerate = int(orbit_degenerate)
        pass

    def get_symmetry(self):
        return self.orbit_sym

    def get_type(self):
        return self.orbit_type

    def get_energy(self):
        return self.orbit_energy

    def get_fraction(self):
        return self.orbit_frac

    def get_degeneracy(self):
        return self.orbit_degenerate

    def __str__(self):
        return str(self.orbit_sym) + ' ' + str(self.orbit_type) + \
               ' ' +str(self.orbit_energy) + ' ' + \
               str(self.orbit_degenerate) + \
               ' ' + str(self.orbit_frac)

    def __repr__(self):
        return str(self.orbit_sym) + ' ' + str(self.orbit_type) + \
               ' ' +str(self.orbit_energy) + ' ' + \
               str(self.orbit_degenerate) + ' ' + \
               str(self.orbit_frac)

    def __eq__(self, other):
        is_equal = True
        if self.orbit_sym != other.orbit_sym:
            is_equal = False
        if self.orbit_type != other.orbit_type:
            is_equal = False
        if abs(self.orbit_energy - other.orbit_energy) > 0.000001:
            is_equal = False
        if abs(self.orbit_frac - other.orbit_frac) > 0.01:
            is_equal = False
        if self.orbit_degenerate != other.orbit_degenerate:
            is_equal = False

        return is_equal

    def as_dict(self) -> dict:
        d = {'orbit_degenderate': self.orbit_degenerate,
             'orbit_energy': self.orbit_energy,
             'orbit_frac': self.orbit_frac,
             'orbit_sym': self.orbit_sym,
             'orbit_type': self.orbit_type}
        return jsanitize(d, strict=True)

