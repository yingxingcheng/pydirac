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

from math import inf
from monty.json import MSONable, jsanitize
from enum import Enum
import numpy as np


class OrbitalType(Enum):
    CLOSED_SHELL = 0
    OPEN_SHELL = 1
    VIRTUAL_SHELL = 2
    ALL_SHELL = 3


class AtomicOrbital(MSONable):
    """Orbit object to restore basis information about given orbit."""

    def __init__(self, sym: str, type: OrbitalType, energy: float, degen: int, frac: float):
        """Create a AtomicOrbital object

        Args:
            sym: orbit symmetry
            type: orbit type, e.g., open-shell or closed-shell
            energy:  orbit energy (a.u)
            degen: orbit degeneration
            frac:
        """
        self.sym = sym
        self.type = type
        self.energy = float(energy)
        self.frac = float(frac)
        self.degen = int(degen)

    def is_valid(self):
        if self.type == OrbitalType.CLOSED_SHELL:
            return np.allclose(self.frac, 1.0)
        elif self.type is OrbitalType.VIRTUAL_SHELL:
            return np.allclose(self.frac, 0.0)
        else:
            return True

    def __str__(self):
        info_list = [self.sym, self.type, self.energy, self.degen, self.frac]
        return " ".join([str(info) for info in info_list])

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        # TODO: why do we need to compare two orbit if they are equal?
        is_equal = True
        if self.sym != other.sym:
            is_equal = False
        if self.type != other.type:
            is_equal = False
        if abs(self.energy - other.energy) > 0.000001:
            is_equal = False
        if abs(self.frac - other.frac) > 0.01:
            is_equal = False
        if self.degen != other.degen:
            is_equal = False
        return is_equal

    def as_dict(self) -> dict:
        ao_type = "closed_shell"
        if self.type is OrbitalType.OPEN_SHELL:
            ao_type = "open_shell"
        elif self.type is OrbitalType.VIRTUAL_SHELL:
            ao_type = "virtual_shell"

        d = {
            "degen": self.degen,
            "energy": self.energy,
            "frac": self.frac,
            "sym": self.sym,
            "type": ao_type,
        }
        return jsanitize(d, strict=True)

    def from_dict(cls, d):
        return AtomicOrbital(**d)


class MoleculeOrbitals(MSONable):
    def __init__(self, aos=None):
        self.aos = aos
        self.closed_shells = []
        self.open_shells = []
        self.virtual_shells = []
        self.occ_open_shell = 0.0

        if self.aos:
            self.analyse()

    def analyse(self):
        # TODO: can we get if this is an atomic calculation
        for ao in self.aos:
            # TODO: do we need to check degeneration?
            if ao.type is OrbitalType.CLOSED_SHELL:
                self.closed_shells.append(ao)
            elif ao.type is OrbitalType.OPEN_SHELL:
                self.open_shells.append(ao)
            elif ao.type is OrbitalType.VIRTUAL_SHELL:
                self.virtual_shells.append(ao)
            else:
                raise RuntimeError("Unknown orbital type !")

    def type_found(self, test_type, ao_type=OrbitalType.ALL_SHELL):
        if ao_type is OrbitalType.ALL_SHELL:
            return True
        else:
            return test_type is ao_type

    def nelec(self, e_min=-inf, e_max=inf, ao_type=OrbitalType.ALL_SHELL) -> int:
        """Count the number of orbital within the energy range [e_min, e_max]

        Args:
            e_min: the minimum of energy
            e_max: the maximum of energy

        Returns:
            The number of electrons located in the energy range
        """
        nb_elec = sum(
            [
                ao.degen * ao.frac
                for ao in self.aos
                if self.type_found(ao.type, ao_type) and e_min <= ao.energy <= e_max
            ]
        )
        return round(nb_elec)

    def nb_closed_elec(self, e_min=-inf, e_max=inf) -> int:
        """The number of electrons in closed-shell AOs"""
        return self.nelec(e_min, e_max, ao_type=OrbitalType.CLOSED_SHELL)

    def nb_open_elec(self, e_min=-inf, e_max=inf) -> int:
        """The number of electrons in open-shell AOs"""
        return self.nelec(e_min, e_max, ao_type=OrbitalType.OPEN_SHELL)

    def nb_virtual_elec(self) -> int:
        """The number of electrons in virtual AOs"""
        return 0

    def nao(self, e_min=-inf, e_max=inf, ao_type=OrbitalType.ALL_SHELL) -> int:
        """The number of AOs by specified AO type and energy range.

        Args:
            e_min: the minimum of energy
            e_max: the maximum of energy
            ao_type: atomic orbital type

        Returns:
            The number of AOs
        """
        nb_ao = sum(
            [
                ao.degen
                for ao in self.aos
                if self.type_found(ao.type, ao_type) and e_min <= ao.energy <= e_max
            ]
        )
        return round(nb_ao)

    def nb_open_ao(self, e_min=-inf, e_max=inf) -> int:
        return self.nao(e_min, e_max, ao_type=OrbitalType.OPEN_SHELL)

    def nb_virtual_ao(self, e_min=-inf, e_max=inf) -> int:
        """Count the number of virtual orbital within the energy
        range [e_min, e_max]

        Args:
            e_min: the minimum of energy
            e_max: the maximum of energy

        Returns:
            number of virtual orbitals
        """
        return self.nao(e_min, e_max, ao_type=OrbitalType.VIRTUAL_SHELL)

    def nb_closed_ao(self, e_min=-inf, e_max=inf) -> int:
        """Count the number of closed orbitals within the energy
        range [e_min, e_max]

        Args:
            e_min: the minimum of energy
            e_max: the maximum of energy

        Returns:
            number of virtual orbitals
        """
        return self.nao(e_min, e_max, ao_type=OrbitalType.CLOSED_SHELL)

    def get_ao_and_elec(self, e_min, e_max):
        # obtain the number of electrons
        nb_closed_elec = self.nb_closed_elec(e_min, e_max)
        nb_open_elec = self.nb_open_elec(e_min, e_max)
        nb_total_elec = self.nelec(e_min, e_max)

        # obtain orbit info
        ## res_orbits = [nb_closed_shells, nb_open_shells, nb_virtual_shells]
        nb_closed_shell = self.nb_closed_ao(e_min, e_max)
        nb_open_shell = self.nb_open_ao(e_min, e_max)
        nb_vir_shell = self.nb_virtual_ao(e_min, e_max)
        return (
            nb_closed_elec,
            nb_open_elec,
            nb_total_elec,
            nb_closed_shell,
            nb_open_shell,
            nb_vir_shell,
        )

    # @classmethod
    # def from_file(cls, filename):
    #     """Get eigenvalues of different symmetry based on 'Eigenvalues' section.

    #     """
    #     if filename.endswith('.out'):
    #         return Molecule.from_dirac_output(filename=filename)
    #     else:
    #         raise NotImplementedError('Not support other format but "*.out"')

    def as_dict(self) -> dict:
        d = {
            "closed_shells": self.closed_shells,
            "occ_open_shell": self.occ_open_shell,
            "open_shells": self.open_shells,
            "virtual_shells": self.virtual_shells,
        }
        return jsanitize(d, strict=True)
