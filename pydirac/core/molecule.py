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

import re
from monty.json import MSONable, jsanitize
from pydirac.core.periodic_table import Element
from pydirac.core.molecular_orbitals import AtomicOrbital


class Molecule(MSONable):

    def __init__(self, OrbitInfo_list=None):
        self.closed_shells = []
        self.open_shells = []
        self.virtual_shells = []
        self.occ_open_shell = 0.0

        if OrbitInfo_list:
            OrbitInfo_list = list(OrbitInfo_list)
            for o_info in OrbitInfo_list:
                self._extract_info(o_info)

        total_nb_elec = self.total_nb_elec()
        if total_nb_elec:
            self.info = Element(total_nb_elec)
        else:
            self.info = 'NULL'

    def _extract_info(self, o_info):
        if 'close' in o_info.get_type():
            self.closed_shells.append(o_info)
            # if o_info not in self.closed_shells:
            #     self.closed_shells.append(o_info)
            # else:
            #     raise warnings.warn("I've met a similar orbit "
            #                          "in closed_shells, SKIP it!")
        elif 'open' in o_info.get_type():
            self.open_shells.append(o_info)
            # if o_info not in self.open_shells:
            #     self.open_shells.append(o_info)
            # else:
            #     raise warnings.warn("I've met a similar orbit "
            #                          "in open_shells, SKIP it!")
        elif 'virtual' in o_info.get_type():
            self.virtual_shells.append(o_info)
            # if o_info not in self.virtual_shells:
            #     self.virtual_shells.append(o_info)
            # else:
            #     raise warnings.warn("I've met a similar orbit "
            #                          "in virtual_shells, SKIP it!")
        else:
            raise RuntimeError('There is a unknown type orbit!')

    def add_orbit(self, orbit):
        self._extract_info(orbit)

    def electron_count(self, e_min, e_max):
        """
        Count electron number if energy range given
        Parameters
        ----------
        e_min
        e_max

        Returns
        -------

        """

        nb_open_shell = int(round(sum([o.orbit_degenerate *
                                       o.orbit_frac for o in self.open_shells if
                                       e_min <= o.orbit_energy <= e_max])))

        nb_closed_shell = int(round(sum([o.orbit_degenerate *
                                         o.orbit_frac for o in
                                         self.closed_shells if
                                         e_min <= o.orbit_energy <= e_max])))

        return nb_open_shell + nb_closed_shell

    def orbit_count(self, res_all=False, min_e=-999999999,
                    max_e=999999999, res_closed=False,
                    c_min_e=0.0, c_max_e=0.0, res_open=False,
                    o_min_e=0.0, o_max_e=0.0,
                    res_virtual=False, v_min_e=0.0, v_max_e=0.0):
        """
        Get orbital info based on energy range specified by user.
        :param res_all: (bool) if true, all orbitals will be taken into consider
        :param min_e:
        :param max_e:
        :param res_closed:
        :param c_min_e:
        :param c_max_e:
        :param res_open:
        :param o_min_e:
        :param o_max_e:
        :param res_virtual:
        :param v_min_e:
        :param v_max_e:
        :return:
        """

        if res_closed:
            res_closed_shells = [o for o in self.closed_shells
                                 if c_min_e <= o.orbit_energy <= c_max_e]
        elif res_all:
            res_closed_shells = [o for o in self.closed_shells
                                 if min_e <= o.orbit_energy <= max_e]
        else:
            res_closed_shells = [o for o in self.closed_shells]

        if res_open:
            res_open_shells = [o for o in self.open_shells
                               if o_min_e <= o.orbit_energy <= o_max_e]
        elif res_all:
            res_open_shells = [o for o in self.open_shells
                               if min_e <= o.orbit_energy <= max_e]
        else:
            res_open_shells = [o for o in self.open_shells]

        if res_virtual:
            res_virtual_shells = [o for o in self.virtual_shells
                                  if v_min_e <= o.orbit_energy <= v_max_e]
        elif res_all:
            res_virtual_shells = [o for o in self.virtual_shells
                                  if min_e <= o.orbit_energy <= max_e]
        else:
            res_virtual_shells = [o for o in self.virtual_shells]

        c_closed = sum([o.orbit_degenerate for o in res_closed_shells])
        c_open = sum([o.orbit_degenerate for o in res_open_shells])
        c_virtual = sum([o.orbit_degenerate for o in res_virtual_shells])

        print('The number of closed-shell orbits is : {0}'.format(c_closed))
        print('The number of open-shell orbits is : {0}'.format(c_open))
        print('The number of virtual orbits is : {0}'.format(c_virtual))
        print('The total number of orbits is: {0}'.format(sum([c_closed,
                                                               c_open,
                                                               c_virtual])))
        return c_closed, c_open, c_virtual

    def closed_elec(self, verbos=True):
        nb_elec = int(round(sum([o.orbit_degenerate * o.orbit_frac
                                 for o in self.closed_shells])))
        if verbos:
            print(
                'The number of closed-shell electrons is : {0}'.format(nb_elec))
        return nb_elec

    def openshell_elec(self, verbos=True):
        nb_elec = int(round(sum([o.orbit_degenerate *
                                 o.orbit_frac for o in self.open_shells])))
        if verbos:
            print('The number of open-shell electrons is : {0}'.format(nb_elec))
        return nb_elec

    def total_nb_elec(self):
        """
        Get total number of electrons
        """
        return self.closed_elec(verbos=False) + self.openshell_elec(
            verbos=False)

    def virtualshell_elec(self):
        print('The number of virtual electrons is : 0')
        # return sum([ o.orbit_degenerate for o in self.virtual_shells])
        return 0

    @classmethod
    def from_file(cls, filename):
        """Get eigenvalues of different symmetry based on 'Eigenvalues' section.
        """

        start_line = re.compile(r"\s+Eigenvalues\s+")
        symline = re.compile(r"^\* \b(?:Boson|Fermion)\b symmetry (.*)")
        open_shell = re.compile(r"^\s+\*\s+Open shell #\d+, f = (\d+(\.\d*)?)")
        close_shell = re.compile(r"^\s+\*\s+Closed shell, f = (\d+(\.\d*)?)")
        virtual_shell = re.compile(
            r"\s+\*\s+Virtual eigenvalues, f = (\d+(\.\d*)?)")
        number = re.compile(
            r"(?P<num>[-+]?(\d+(\.\d*)?|\d*\.\d+))\s+\(\s*(?P<degenerate>\d+)\)")
        number_line = re.compile(r"\s+[-+]?(\d+(\.\d+)?|\d*\.\d+) .*")
        endline = re.compile(r"^\* HOMO - LUMO")

        with open(filename, 'r') as f:
            lines = f.readlines()

        for i, l in enumerate(lines):
            if start_line.match(l):
                #        print(l)
                break
        # print(i)
        lines = lines[i + 1:]

        info_list = {}
        cur_sym = ''
        o_type = ''
        occ_frac = 0.00
        res_list = []
        for l in lines:
            symmetry_match = symline.match(l)
            if symmetry_match:
                # add new branch
                cur_sym = symmetry_match.group(1)
                if cur_sym in info_list.keys():
                    pass
                else:
                    info_list[cur_sym] = {}

            openshell_match = open_shell.match(l)
            if openshell_match:
                occ_frac = float(openshell_match.group(1))
                o_type = 'open-shell'
                if o_type in info_list[cur_sym].keys():
                    pass
                else:
                    info_list[cur_sym][o_type] = {}

            closeshell_match = close_shell.match(l)
            if closeshell_match:
                occ_frac = float(closeshell_match.group(1))
                # occ_frac = 1.0
                o_type = 'close-shell'
                if o_type in info_list[cur_sym].keys():
                    pass
                else:
                    info_list[cur_sym][o_type] = {}

            virtualshell_match = virtual_shell.match(l)
            if virtualshell_match:
                occ_frac = float(virtualshell_match.group(1))
                # occ_frac = 0.0
                o_type = 'virtual-shell'
                if o_type in info_list[cur_sym].keys():
                    pass
                else:
                    info_list[cur_sym][o_type] = {}

            if number_line.match(l):
                line_nums = [float(m.group("num")) for m in number.finditer(l)]
                de_list = [float(m.group("degenerate")) for m in
                           number.finditer(l)]
                for e, d in zip(line_nums, de_list):
                    oinfo = AtomicOrbital(cur_sym, o_type, e, d, occ_frac)
                    res_list.append(oinfo)

            if endline.match(l):
                break

        return Molecule(res_list)

    def as_dict(self) -> dict:
        info_dict = {}
        if type(self.info) == Element:
            for attr in ['atomic_number', 'symbol', 'dipole_polarizability']:
                info_dict[attr] = getattr(self.info, attr)
            # TODO:
            # for k in [ key for key in dir(self.info) if not key.startswith('_')]:
            #     if not callable(getattr(self.info, k)):
            #         info_dict[k] = getattr(self.info, k)
        else:
            info_dict = 'null'

        d = {'closed_shells': self.closed_shells,
             'info': info_dict,
             'occ_open_shell': self.occ_open_shell,
             'open_shells': self.open_shells,
             'virtual_shells': self.virtual_shells}
        return jsanitize(d, strict=True)


if __name__ == '__main__':
    import sys

    argv = sys.argv[1:]

    for fname in argv:
        Atom1 = Molecule.from_file(fname)
        print('Original orbit information'.center(80, '-'))
        Atom1.orbit_count()

        Atom1.closed_elec()
        Atom1.openshell_elec()
        Atom1.virtualshell_elec()
        # for i in Atom1.virtual_shells:
        #     print(i)
        # print('-'*80)

        print('')
        print('Selected orbit information'.center(80, '-'))
        Atom1.orbit_count(res_virtual=True, v_min_e=-1, v_max_e=10.0,
                          res_closed=True, c_min_e=-3.0, c_max_e=0)


