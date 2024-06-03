from enum import Enum
from math import inf

import numpy as np
from monty.json import MSONable, jsanitize

__all__ = ["OrbitalType", "AtomicOrbital", "MoleculeOrbitals"]


class OrbitalType(Enum):
    """
    An enumeration class representing the type of an orbital.

    Attributes
    ----------
    CLOSED_SHELL : int
        A closed shell orbital.
    OPEN_SHELL : int
        An open shell orbital.
    VIRTUAL_SHELL : int
        A virtual orbital.
    ALL_SHELL : int
        A combination of closed, open, and virtual orbitals.

    """

    CLOSED_SHELL = 0
    OPEN_SHELL = 1
    VIRTUAL_SHELL = 2
    ALL_SHELL = 3


class AtomicOrbital(MSONable):
    """
    A class representing an atomic orbital.

    Parameters
    ----------
    sym : str
        The symmetry of the orbital.
    type : OrbitalType
        The type of the orbital.
    energy : float
        The energy of the orbital in atomic units.
    degen : int
        The degeneracy of the orbital.
    frac : float
        The fractional occupancy of the orbital.

    Attributes
    ----------
    sym : str
        The symmetry of the orbital.
    type : OrbitalType
        The type of the orbital.
    energy : float
        The energy of the orbital in atomic units.
    degen : int
        The degeneracy of the orbital.
    frac : float
        The fractional occupancy of the orbital.

    Methods
    -------
    is_valid() -> bool
        Checks if the atomic orbital is valid.
    __str__() -> str
        Returns a string representation of the atomic orbital.
    __repr__() -> str
        Returns a string representation of the atomic orbital.
    __eq__(other: AtomicOrbital) -> bool
        Checks if two atomic orbitals are equal.
    as_dict() -> dict
        Returns a dictionary representation of the atomic orbital.
    from_dict(cls, d: dict) -> AtomicOrbital
        Returns an AtomicOrbital object created from a dictionary.

    """

    def __init__(self, sym: str, type: OrbitalType, energy: float, degen: int, frac: float):
        """
        Initialize an AtomicOrbital object.

        Parameters
        ----------
        sym : str
            The symmetry of the orbital.
        type : OrbitalType
            The type of the orbital.
        energy : float
            The energy of the orbital in atomic units.
        degen : int
            The degeneracy of the orbital.
        frac : float
            The fractional occupancy of the orbital.
        """
        self.sym = sym
        self.type = type
        self.energy = float(energy)
        self.frac = float(frac)
        self.degen = int(degen)

    def is_valid(self) -> bool:
        """
        Check if the atomic orbital is valid.

        Returns
        -------
        bool
            True if the atomic orbital is valid, False otherwise.
        """

        if self.type == OrbitalType.CLOSED_SHELL:
            return np.allclose(self.frac, 1.0)
        elif self.type is OrbitalType.VIRTUAL_SHELL:
            return np.allclose(self.frac, 0.0)
        else:
            return True

    def __str__(self) -> str:
        """
        Return a string representation of the atomic orbital.

        Returns
        -------
        str
            A string representation of the atomic orbital.
        """

        info_list = [self.sym, self.type, self.energy, self.degen, self.frac]
        return " ".join([str(info) for info in info_list])

    def __repr__(self) -> str:
        """
        Return a string representation of the atomic orbital.

        Returns
        -------
        str
            A string representation of the atomic orbital.
        """

        return self.__str__()

    def __eq__(self, other) -> bool:
        """
        Check if two atomic orbitals are equal.

        Parameters
        ----------
        other : AtomicOrbital
            The other atomic orbital.

        Returns
        -------
        bool
            True if the two atomic orbitals are equal, False
        """
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
        """
        Return a dictionary representation of the atomic orbital.

        Returns
        -------
        dict
            A dictionary with the following keys:
            - sym: str, the symmetry of the orbital
            - type: str, the type of the orbital
            - energy: float, the energy of the orbital in atomic units
            - degen: int, the degeneracy of the orbital
            - frac: float, the fractional occupancy of the orbital
        """

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

    @classmethod
    def from_dict(cls, d):
        """
        Create an `AtomicOrbital` object from a dictionary.

        Args:
            d (dict): The dictionary containing the `AtomicOrbital` object data.

        Returns:
            An `AtomicOrbital` object instance.
        """
        return AtomicOrbital(**d)


class MoleculeOrbitals(MSONable):
    """
    A class representing the orbitals of a molecule.

    Parameters
    ----------
    aos : list of AtomicOrbital objects, optional
        The atomic orbitals of the molecule.

    Attributes
    ----------
    aos : list of AtomicOrbital objects
        The atomic orbitals of the molecule.
    closed_shells : list of AtomicOrbital objects
        The closed-shell atomic orbitals of the molecule.
    open_shells : list of AtomicOrbital objects
        The open-shell atomic orbitals of the molecule.
    virtual_shells : list of AtomicOrbital objects
        The virtual atomic orbitals of the molecule.
    occ_open_shell : float
        The fractional occupancy of the open-shell orbitals.

    Methods
    -------
    analyse()
        Analyzes the atomic orbitals and separates them into different categories.
    type_found(test_type: OrbitalType, ao_type: OrbitalType) -> bool
        Checks if the atomic orbital type matches the specified type.
    nelec(e_min=-inf, e_max=inf, ao_type=OrbitalType.ALL_SHELL) -> int
        Counts the number of electrons in the specified energy range and atomic orbital type.
    nb_closed_elec(e_min=-inf, e_max=inf) -> int
        Counts the number of electrons in the closed-shell atomic orbitals within the specified energy range.
    nb_open_elec(e_min=-inf, e_max=inf) -> int
        Counts the number of electrons in the open-shell atomic orbitals within the specified energy range.
    nb_virtual_elec() -> int
        Counts the number of electrons in the virtual atomic orbitals.
    nao(e_min=-inf, e_max=inf, ao_type=OrbitalType.ALL_SHELL) -> int
        Counts the number of atomic orbitals in the specified energy range and atomic orbital type.
    nb_open_ao(e_min=-inf, e_max=inf) -> int
        Counts the number of open-shell atomic orbitals within the specified energy range.
    nb_virtual_ao(e_min=-inf, e_max=inf) -> int
        Counts the number of virtual atomic orbitals within the specified energy range.
    nb_closed_ao(e_min=-inf, e_max=inf) -> int
        Counts the number of closed-shell atomic orbitals within the specified energy range.
    get_ao_and_elec(e_min, e_max)
        Returns the number of electrons and atomic orbitals in different categories within the specified energy range.
    as_dict() -> dict
        Returns a dictionary representation of the MoleculeOrbitals object.
    """

    def __init__(self, aos=None):
        """
        Initialize a MoleculeOrbitals object.

        Parameters
        ----------
        aos : list of AtomicOrbital objects, optional
            The atomic orbitals of the molecule.
        """
        self.aos = aos
        self.closed_shells = []
        self.open_shells = []
        self.virtual_shells = []
        self.occ_open_shell = 0.0

        if self.aos:
            self.analyse()

    def analyse(self):
        """
        Analyze the atomic orbitals and separates them into different categories.
        """
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
        """
        Check if the atomic orbital type matches the specified type.

        Parameters
        ----------
        test_type : OrbitalType
            The atomic orbital type to check against.
        ao_type : OrbitalType, optional
            The specified atomic orbital type to match against. Defaults to OrbitalType.ALL_SHELL.

        Returns
        -------
        bool
            True if the atomic orbital type matches the specified type, False otherwise.
        """
        if ao_type is OrbitalType.ALL_SHELL:
            return True
        else:
            return test_type is ao_type

    def nelec(self, e_min=-inf, e_max=inf, ao_type=OrbitalType.ALL_SHELL) -> int:
        """
        Count the number of electrons in the specified energy range and atomic orbital type.

        Parameters
        ----------
        e_min : float, optional
            The minimum energy of the range. Defaults to -inf.
        e_max : float, optional
            The maximum energy of the range. Defaults to inf.
        ao_type : OrbitalType, optional
            The specified atomic orbital type to count. Defaults to OrbitalType.ALL_SHELL.

        Returns
        -------
        int
            The number of electrons located in the energy range.
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
        """
        Count the number of electrons in the closed-shell atomic orbitals within the specified energy range.

        Parameters
        ----------
        e_min : float, optional
            The minimum energy of the range. Defaults to -inf.
        e_max : float, optional
            The maximum energy of the range. Defaults to inf.

        Returns
        -------
        int
            The number of electrons located in the energy range in the closed-shell atomic orbitals.
        """
        return self.nelec(e_min, e_max, ao_type=OrbitalType.CLOSED_SHELL)

    def nb_open_elec(self, e_min=-inf, e_max=inf) -> int:
        """
        Count the number of electrons in the open-shell atomic orbitals within the specified energy range.

        Parameters
        ----------
        e_min : float, optional
            The minimum energy of the range. Defaults to -inf.
        e_max : float, optional
            The maximum energy of the range. Defaults to inf.

        Returns
        -------
        int
            The number of electrons located in the energy range in the open-shell atomic orbitals.
        """
        return self.nelec(e_min, e_max, ao_type=OrbitalType.OPEN_SHELL)

    def nb_virtual_elec(self) -> int:
        """
        Count the number of electrons in the virtual atomic orbitals.

        Returns
        -------
        int
            The number of electrons in the virtual atomic orbitals.
        """
        return 0

    def nao(self, e_min=-inf, e_max=inf, ao_type=OrbitalType.ALL_SHELL) -> int:
        """
        Count the number of atomic orbitals in the specified energy range and atomic orbital type.

        Parameters
        ----------
        e_min : float, optional
            The minimum energy of the range. Defaults to -inf.
        e_max : float, optional
            The maximum energy of the range. Defaults to inf.
        ao_type : OrbitalType, optional
            The specified atomic orbital type to count. Defaults to OrbitalType.ALL_SHELL.

        Returns
        -------
        int
            The number of atomic orbitals located in the energy range.
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
        """
        Count the number of open-shell atomic orbitals within the specified energy range.

        Parameters
        ----------
        e_min : float, optional
            The minimum energy of the range. Defaults to -inf.
        e_max : float, optional
            The maximum energy of the range. Defaults to inf.

        Returns
        -------
        int
            The number of open-shell atomic orbitals located in the energy range.
        """
        return self.nao(e_min, e_max, ao_type=OrbitalType.OPEN_SHELL)

    def nb_virtual_ao(self, e_min=-inf, e_max=inf) -> int:
        """
        Count the number of virtual atomic orbitals within the specified energy range.

        Parameters
        ----------
        e_min : float, optional
            The minimum energy of the range. Defaults to -inf.
        e_max : float, optional
            The maximum energy of the range. Defaults to inf.

        Returns
        -------
        int
            The number of virtual atomic orbitals located in the energy range.
        """
        return self.nao(e_min, e_max, ao_type=OrbitalType.VIRTUAL_SHELL)

    def nb_closed_ao(self, e_min=-inf, e_max=inf) -> int:
        """
        Count the number of closed-shell atomic orbitals within the specified energy range.

        Parameters
        ----------
        e_min : float, optional
            The minimum energy of the range. Defaults to -inf.
        e_max : float, optional
            The maximum energy of the range. Defaults to inf.

        Returns
        -------
        int
            The number of closed-shell atomic orbitals located in the energy range.
        """
        return self.nao(e_min, e_max, ao_type=OrbitalType.CLOSED_SHELL)

    def get_ao_and_elec(self, e_min, e_max):
        """
        Obtain the number of electrons and atomic orbitals in the specified energy range.

        Parameters
        ----------
        e_min : float
            The minimum energy of the range.
        e_max : float
            The maximum energy of the range.

        Returns
        -------
        tuple
            A tuple containing the number of closed-shell electrons, open-shell electrons, total electrons,
            closed-shell atomic orbitals, open-shell atomic orbitals, and virtual atomic orbitals.
        """
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

    def as_dict(self) -> dict:
        """
        Return a dictionary representation of the MoleculeOrbitals object.

        Returns
        -------
        dict
            A dictionary containing the closed-shell, open-shell, and virtual atomic orbitals, as well as
            the occupation number of open-shell atomic orbitals.
        """
        d = {
            "closed_shells": self.closed_shells,
            "occ_open_shell": self.occ_open_shell,
            "open_shells": self.open_shells,
            "virtual_shells": self.virtual_shells,
        }
        return jsanitize(d, strict=True)
