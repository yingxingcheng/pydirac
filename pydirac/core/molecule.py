from monty.json import MSONable
from pydirac.core.periodic import Element


class Molecule(MSONable):
    def __init__(self, atoms=None, corrdinates=None):

        self.atoms = atoms or []
        if len(self.atoms) == 1:
            self.is_atom = True
        else:
            self.is_atom = False

        self.coordinates = corrdinates or []
        assert len(self.atoms) == len(self.coordinates)

        self.atomic_info = "null"
        if self.is_atom:
            self.atomic_info = Element(self.atoms[0])

    def as_dict(self) -> dict:
        d = {
            "atoms": self.atoms,
            "coordinates": self.coordinates,
            "is_atom": self.is_atom,
            "atomic_info": self.atomic_info.as_dict(),
        }
        return d
