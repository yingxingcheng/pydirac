from monty.json import MSONable

from pydirac.core.periodic import Element

__all__ = ["Molecule"]


class Molecule(MSONable):
    """
    A class representing a molecule.

    Parameters
    ----------
    atoms : list, optional
        A list of atomic symbols representing the atoms in the molecule.
    coordinates : list, optional
        A list of tuples representing the cartesian coordinates of the atoms in the molecule.

    Attributes
    ----------
    atoms : list
        A list of atomic symbols representing the atoms in the molecule.
    coordinates : list
        A list of tuples representing the cartesian coordinates of the atoms in the molecule.
    is_atom : bool
        A boolean indicating whether the molecule is a single atom.
    atomic_info : str or Element
        A string representing the atomic symbol or an Element object representing the atomic properties
        of the atom in the molecule. If the molecule is not a single atom, atomic_info is set to "null".

    Methods
    -------
    as_dict() -> dict
        Returns a dictionary representation of the molecule.

    """

    def __init__(self, atoms=None, coordinates=None):
        """
        Initialize a Molecule object.

        Raises
        ------
        AssertionError
            If the length of the atoms and coordinates lists do not match.
        """

        self.atoms = atoms or []
        if len(self.atoms) == 1:
            self.is_atom = True
        else:
            self.is_atom = False

        self.coordinates = coordinates or []
        assert len(self.atoms) == len(self.coordinates)

        self.atomic_info = "null"
        if self.is_atom:
            self.atomic_info = Element(self.atoms[0])

    def as_dict(self) -> dict:
        """
        Returns a dictionary representation of the molecule.

        Returns
        -------
        dict
            A dictionary with the following keys:
            - atoms: list of atomic symbols
            - coordinates: list of tuples representing the cartesian coordinates of the atoms in the molecule
            - is_atom: boolean indicating whether the molecule is a single atom
            - atomic_info: dictionary representing the atomic properties of the atom in the molecule, if the molecule is a single atom.
        """
        d = {
            "atoms": self.atoms,
            "coordinates": self.coordinates,
            "is_atom": self.is_atom,
            "atomic_info": self.atomic_info.as_dict() if self.is_atom else "null",
        }
        return d
