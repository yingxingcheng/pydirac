#!/usr/bin/env python

"""
Function: Extract orbital information from SCF calculation output file.
Author: Yingxing Cheng
Date: 10/20/2019
"""

import re


class OrbitInfo(object):
    """
    Orbit object to restore basis information about given orbit.
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
        return str(self.orbit_sym) + ' ' + str(self.orbit_type) + ' ' +str(self.orbit_energy) + ' ' + str(self.orbit_degenerate) + ' ' + str(self.orbit_frac)

    def __repr__(self):
        return str(self.orbit_sym) + ' ' + str(self.orbit_type) + ' ' +str(self.orbit_energy) + ' ' + str(self.orbit_degenerate) + ' ' + str(self.orbit_frac)

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


class Atom(object):

    def __init__(self, OrbitInfo_list = None):
        self.closed_shells = []
        self.open_shells = []
        self.virtual_shells = []
        self.occ_open_shell = 0.0

        if OrbitInfo_list:
            OrbitInfo_list = list(OrbitInfo_list)
            for o_info in OrbitInfo_list:
                self._extract_info(o_info)

    def _extract_info(self, o_info):
        if 'close' in o_info.get_type():
            if o_info not in self.closed_shells:
                self.closed_shells.append(o_info)
            else:
                raise RuntimeWarning("I've met a similar orbit in closed_shells, SKIP it!")
        elif 'open' in o_info.get_type():
            if o_info not in self.open_shells:
                self.open_shells.append(o_info)
            else:
                raise RuntimeWarning("I've met a similar orbit in open_shells, SKIP it!")
        elif 'virtual' in o_info.get_type():
            if o_info not in self.virtual_shells:
                self.virtual_shells.append(o_info)
            else:
                raise RuntimeWarning("I've met a similar orbit in virtual_shells, SKIP it!")
        else:
            raise RuntimeError('There is a unknown type orbit!')


    def add_orbit(self, orbit):
        self._extract_info(orbit)


    def orbit_count(self,  res_all = False, min_e = -999999999, max_e= 999999999,res_closed=False,
                    c_min_e = 0.0, c_max_e = 0.0, res_open =False, o_min_e = 0.0, o_max_e = 0.0,
                    res_virtual=  False,v_min_e = 0.0, v_max_e = 0.0):

        if res_closed:
            res_closed_shells = [o for o in self.closed_shells if o.orbit_energy >= c_min_e
                                 and o.orbit_energy <=c_max_e]
        elif res_all:
            res_closed_shells = [o for o in self.closed_shells if o.orbit_energy >= min_e
                                 and o.orbit_energy <= max_e]
        else:
            res_closed_shells = [o for o in self.closed_shells]


        if res_open:
            res_open_shells = [o for o in self.open_shells if o.orbit_energy >=o_min_e
                               and o.orbit_energy <= o_max_e]
        elif res_all:
            res_open_shells = [o for o in self.open_shells if o.orbit_energy >= min_e
                               and o.orbit_energy <= max_e]
        else:
            res_open_shells = [o for o in self.open_shells]


        if res_virtual:
            res_virtual_shells= [o for o in self.virtual_shells if o.orbit_energy >= v_min_e
                                 and o.orbit_energy <= v_max_e]
        elif res_all:
            res_virtual_shells= [o for o in self.virtual_shells if o.orbit_energy >= min_e
                                 and o.orbit_energy <= max_e]
        else:
            res_virtual_shells= [o for o in self.virtual_shells]


        c_closed = sum([ o.orbit_degenerate for o in res_closed_shells])
        c_open = sum([ o.orbit_degenerate for o in res_open_shells])
        c_virtual = sum([ o.orbit_degenerate for o in res_virtual_shells])
        
        print('The number of closed-shell orbits is : {0}'.format(c_closed))
        print('The number of open-shell orbits is : {0}'.format(c_open))
        print('The number of virtual orbits is : {0}'.format(c_virtual))
        print('The total number of orbits is: {0}'.format(sum([c_closed, c_open,c_virtual])))
        return c_closed, c_open, c_virtual

    def closed_elec(self):
        nb_elec = int(round(sum([ o.orbit_degenerate * o.orbit_frac for o in self.closed_shells])))
        print('The number of closed-shell electrons is : {0}'.format(nb_elec))
        return nb_elec


    def openshell_elec(self):
        nb_elec= int(round(sum([ o.orbit_degenerate * o.orbit_frac for o in self.open_shells])))
        print('The number of open-shell electrons is : {0}'.format(nb_elec))
        return nb_elec
        
    def virtualshell_elec(self):
        print('The number of virtual electrons is : 0')
        #return sum([ o.orbit_degenerate for o in self.virtual_shells])
        return 0

    @classmethod
    def from_file(cls, filename):
        """Get eigenvalues of different symmetry based on 'Eigenvalues' section.
        """
    
        start_line = re.compile(r"\s+Eigenvalues\s+")
        symline = re.compile(r"^\* \b(?:Boson|Fermion)\b symmetry (.*)")
        open_shell = re.compile(r"^\s+\*\s+Open shell \#\d+, f = (\d+(\.\d*)?)")
        close_shell = re.compile(r"^\s+\*\s+Closed shell, f = (\d+(\.\d*)?)")
        virtual_shell = re.compile(r"\s+\*\s+Virtual eigenvalues, f = (\d+(\.\d*)?)")
        number = re.compile(r"(?P<num>[-+]?(\d+(\.\d*)?|\d*\.\d+))\s+\(\s*(?P<degenerate>\d+)\)")
        number_line = re.compile(r"\s+[-+]?(\d+(\.\d+)?|\d*\.\d+) .*")
        endline = re.compile(r"^\* HOMO - LUMO")
    
    
        with open(filename, 'r') as f:
            lines = f.readlines()
    
        for i, l in enumerate(lines):
            if start_line.match(l):
        #        print(l)
                break
        #print(i)
        lines = lines[i+1:]
    
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
                #occ_frac = 1.0
                o_type = 'close-shell'
                if o_type in info_list[cur_sym].keys():
                    pass
                else:
                    info_list[cur_sym][o_type] = {}
    
            virtualshell_match = virtual_shell.match(l)
            if virtualshell_match:
                occ_frac = float(virtualshell_match.group(1))
                #occ_frac = 0.0
                o_type = 'virtual-shell'
                if o_type in info_list[cur_sym].keys():
                    pass
                else:
                    info_list[cur_sym][o_type] = {}
    
    
            if number_line.match(l):
                line_nums = [float(m.group("num")) for m in number.finditer(l)]
                de_list = [float(m.group("degenerate")) for m in number.finditer(l)]
                for e, d in zip(line_nums, de_list):
                    oinfo = OrbitInfo(cur_sym, o_type,e, d, occ_frac )
                    res_list.append(oinfo)
    
            if endline.match(l):
                break
        # print(res_list)
        #fo
        #    print(i)
        return Atom(res_list)


if __name__ == '__main__':
    import sys

    argv = sys.argv[1:]

    for fname in argv:
        Atom1 = Atom.from_file(fname)
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
        Atom1.orbit_count(res_virtual = True, v_min_e = -1, v_max_e = 10.0, res_closed = True, c_min_e = -3.0 , c_max_e = 0)
        
    
