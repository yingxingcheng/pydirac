import six
from monty.json import MSONable, jsanitize

__all__ = ["Element"]
# fmt: off
_PERIODIC_TABLE = [
    'dummy', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S',
    'Cl', 'Ar','K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge',
    'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
    'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
    'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
    'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
]
# fmt: on

_SPECIAL_ELEMENTS = {
    24: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^1 3d^5",
    29: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^1 3d^10",
    41: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^4",
    42: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^5",
    44: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^7",
    45: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^8",
    46: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 4s^10",
    47: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10",
    57: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10 5p^6 6s^2 5d^1",
    58: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10 5p^6 6s^2 4f^1 5d^1",
    64: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10 5p^6 6s^2 4f^7 5d^1",
    78: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^1 4f^14 5d^9",
    79: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^1 4f^14 5d^10",
    89: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 6d^1",
    90: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 6d^2",
    91: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^2 6d^1",
    92: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^3 6d^1",
    93: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^4 6d^1",
    96: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^7 6d^1",
    103: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^14 7p^1",
    110: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^1 5f^14 6d^9",
}


class Element(MSONable):
    """
    Represents a chemical element in the periodic table.

    Attributes
    ----------
    atomic_number : int
        The atomic number of the element.
    symbol : str
        The symbol of the element.
    """

    ROW_SIZES = (2, 8, 8, 18, 18, 32, 32)

    def __init__(self, id):
        """
        Initialize the Element object.

        Parameters
        ----------
        id : str or int
            The symbol or atomic number of the element.
        """
        if isinstance(id, (six.string_types, int)):
            idx, symbol = self._get_element(id)
        else:
            raise ValueError("Expected a <str> or <int>, got: {0:s}".format(type(id)))

        self.atomic_number = idx
        self.symbol = symbol

    def _get_element(self, id):
        """
        Get the element by the given id (either symbol or atomic number).

        Parameters
        ----------
        id : str or int
            The symbol or atomic number of the element.

        Returns
        -------
        tuple
            A tuple containing the atomic number and symbol of the element.
        """

        def get_symbol(id):
            """
            Get element by given index.
            """
            if type(id) == int:
                return _PERIODIC_TABLE[id]
            else:
                raise RuntimeError("Type id should be int")

        def get_id(symbol):
            """
            Get element number by given element type.
            """

            idx = _PERIODIC_TABLE.index(symbol)
            if idx >= 0:
                return idx
            else:
                raise RuntimeError("Wrong element type!")

        if isinstance(id, six.string_types):
            return get_id(id), id
        else:
            return id, get_symbol(id)

    @property
    def block(self):
        """
        Get the block character (s, p, d, or f) of the element.

        Returns
        -------
        str
            The block character of the element.
        """
        if (self.is_actinoid or self.is_lanthanoid) and self.atomic_number not in [71, 103]:
            return "f"
        if self.is_actinoid or self.is_lanthanoid:
            return "d"
        if self.group in [1, 2]:
            return "s"
        if self.group in range(13, 19):
            return "p"
        if self.group in range(3, 13):
            return "d"
        raise ValueError("unable to determine block")

    @property
    def is_lanthanoid(self):
        """
        Check if the element is a lanthanoid.

        Returns
        -------
        bool
            True if the element is a lanthanoid, False otherwise.
        """
        return 56 < self.atomic_number < 72

    @property
    def is_actinoid(self):
        """
        Check if the element is an actinoid.

        Returns
        -------
        bool
            True if the element is an actinoid, False otherwise.
        """
        return 88 < self.atomic_number < 104

    @property
    def row(self):
        """
        Get the periodic table row of the element.

        Returns
        -------
        int
            The row of the element in the periodic table.
        """
        z = self.Z
        total = 0
        if 57 <= z <= 71:
            return 8
        if 89 <= z <= 103:
            return 9
        for i, size in enumerate(Element.ROW_SIZES):
            total += size
            if total >= z:
                return i + 1
        return 8

    @property
    def group(self):
        """
        Get the periodic table group of the element.

        Returns
        -------
        int
            The group of the element in the periodic table.
        """
        z = self.atomic_number
        if z == 1:
            return 1
        if z == 2:
            return 18
        if 3 <= z <= 18:
            if (z - 2) % 8 == 0:
                return 18
            if (z - 2) % 8 <= 2:
                return (z - 2) % 8
            return 10 + (z - 2) % 8

        if 19 <= z <= 54:
            if (z - 18) % 18 == 0:
                return 18
            return (z - 18) % 18

        if (z - 54) % 32 == 0:
            return 18
        if (z - 54) % 32 >= 18:
            return (z - 54) % 32 - 14
        return (z - 54) % 32

    @property
    def group_symbol(self):
        """
        Get the group symbol of the element.

        Returns
        -------
        str
            The group symbol of the element.
        """
        return "group-{0}".format(self.group)

    @property
    def Z(self):
        """
        Get the atomic number of the element.

        Returns
        -------
        int
            The atomic number of the element.
        """
        return self.atomic_number

    @property
    def period(self):
        """
        Get the period of the element in the periodic table.

        Returns
        -------
        int
            The period of the element in the periodic table.
        """
        return self.row

    def as_dict(self) -> dict:
        """
        Get a dictionary representation of the Element object.

        Returns
        -------
        dict
            The dictionary representation of the Element object.
        """
        d = {
            "Z": self.Z,
            "block": self.block,
            "row": self.row,
            "group": self.group,
            "group_symbol": self.group_symbol,
            "period": self.period,
        }
        return jsanitize(d)

    def get_elec_config(self):
        """
        Get the electron configuration of the element.

        Returns
        -------
        str
            The electron configuration of the element.
        """
        n = self.Z
        if n in _SPECIAL_ELEMENTS.keys():
            return _SPECIAL_ELEMENTS[n]

        rule = "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p"
        nb_dict = {"s": 2, "p": 6, "d": 10, "f": 14}
        orbitals = [(i, nb_dict[i[-1]]) for i in rule.split()]
        output = []

        for orbital, size in orbitals:
            k = min(size, n)
            output.append("%s^%d" % (orbital, k))
            n -= k
            if n <= 0:
                break

        orbital_info = " ".join(output)
        return orbital_info

    def get_dirac_shell_str(self, return_dict=False):
        """
        Returns the electron configurations in DIRAC according to the number of electrons.

        Parameters
        ----------
        return_dict : bool, optional
            Whether to return the results as a dictionary or string. Default is False.

        Returns
        -------
        str or dict
            The electron configurations in DIRAC.

        Notes
        -----
        For example:

        - for NSYM = 1
            ```text
            .CLOSED SHELL
            2
            .OPEN SHELL
            1
            1/2
            ```

        - or for NSYM != 1
            ```text
            .CLSED SHELL
            2 0
            .OPEN SHELL
            1
            1/2,0
            ```
        """

        def get_shell(elec_config, NSYM=1):
            """
            Determines the number of electrons in each shell and subshell for a given electronic configuration.

            Parameters
            ----------
            elec_config : str
                The electronic configuration of the atom. For example, "1s2 2s2 2p6" represents the configuration of Neon.
            NSYM : int, optional
                The number of symmetry operations, defaults to 1.

            Returns
            -------
            cs_str : str or None
                The string representation of the number of electrons in closed shells, or None if there are no closed shells.
            os_str : str or None
                The string representation of the number of electrons in open shells, or None if there are no open shells.

            Notes
            -----
            This function uses the Aufbau principle to determine the electron configuration. The order of orbitals is as follows:
            "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p"
            The maximum number of electrons that can be in each subshell is as follows:
            s: 2
            p: 6
            d: 10
            f: 14
            Any subshell with less than the maximum number of electrons is considered open shell, otherwise it is considered closed shell.

            Examples
            --------
            >>> get_shell("1s2 2s2 2p6")
            ('10', None)

            >>> get_shell("1s2 2s2 2p5")
            (None, '1/2,0\\n1/2,0\\n2/6,0')

            """
            rule = "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p"
            nb_dict = {"s": 2, "p": 6, "d": 10, "f": 14}
            cs_even, cs_odd = 0, 0
            os_even, os_odd = [], []
            for i in elec_config.split():
                # 1s^2
                t, nb = i.split("^")
                if int(nb) < nb_dict[t[-1]]:
                    # open_shell
                    if t[-1] in ["s", "d"]:
                        # even
                        os_even.append((t, int(nb)))
                    else:
                        os_odd.append((t, int(nb)))
                else:
                    if t[-1] in ["s", "d"]:
                        # even
                        cs_even += int(nb)
                    else:
                        cs_odd += int(nb)

            cs_total = cs_even + cs_odd
            os_total = os_even + os_odd
            os_total = sorted(os_total, key=lambda x: rule.index(x[0]))
            nb_os = len(os_total)

            cs_str, os_str = [], []
            if NSYM == 1:
                if cs_total == 0:
                    cs_str = None
                else:
                    cs_str.append("{0}".format(cs_total))
                if not len(os_total):
                    os_str = None
                else:
                    os_str.append("{0}".format(nb_os))
                    for i in os_total:
                        os_str.append("{0}/{1}".format(i[-1], nb_dict[i[0][1]]))
            else:
                # for closed shell: 10, 10
                # for open shell:
                # 2
                # 5/10,0
                # 1/2,0
                if cs_total:
                    cs_str.append("{0} {1}\n".format(cs_even, cs_odd))
                else:
                    cs_str = None

                if nb_os:
                    os_str.append("{0}\n".format(nb_os))
                    for i in os_total:
                        if i[0][-1] in ["s", "d"]:
                            os_str.append("{0}/{1},{2}\n".format(i[-1], nb_dict[i[0][1]], 0))
                        else:
                            os_str.append("{0}/{1},{2}\n".format(i[-1], 0, nb_dict[i[0][1]]))
                else:
                    os_str = None
            if cs_str:
                cs_str = "\n".join(cs_str)
            if os_str:
                os_str = "\n".join(os_str)
            return cs_str, os_str

        elec_config = self.get_elec_config()
        closed_shell, open_shell = get_shell(elec_config)

        res_str = []
        res_dict = {}
        res_str.append("# {0} atom ".format(self.symbol) + "and Econf: {0}".format(elec_config))
        if closed_shell:
            res_str.append(".CLOSED SHELL")
            res_str.append(closed_shell)
            res_dict["CLOSED SHELL"] = closed_shell
        if open_shell:
            res_str.append(".OPEN SHELL")
            res_str.append(open_shell)
            res_dict["OPEN SHELL"] = open_shell

        res_str = "\n".join(res_str)
        if return_dict:
            return res_dict
        else:
            return res_str
