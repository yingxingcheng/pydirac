#!/usr/bin/env python

"""
Function: Create input file from the previous output file, e.g., SCF calculation.
Author: Yingxing Cheng
Date: 10/21/2019
"""
from pydirac.get_orbit_info import OrbitInfo, Atom
from pydirac.dirac_input import Inpobj


def create_mrci_inp(fin, fout='PYDIRAC.inp'):
    """
    Create MRCI input based on the previous SCF calculation.
    """

    atom = Atom.from_file(fin)
    # obtain the number of electrons
    closed_elec = atom.closed_elec()
    open_elec = atom.openshell_elec()

    # obtain orbit info
    res_orbits = atom.orbit_count(res_virtual=True, v_min_e=-1, v_max_e=10.0)
    #res_orbits = atom.orbit_count(res_virtual=True, v_min_e=-1, v_max_e=2.0)
    nb_openshell_orbit = res_orbits[1]
    nb_virtual = res_orbits[-1]

    nb_elec = closed_elec + open_elec
    # whether atom contain d orbits
    has_d = True
    if nb_elec <= 20:
        has_d = False

    if has_d:
        elec_in_gas = 10 + 2 + open_elec
        gas1 = '10 12 / 6' 
        gas2 = '{0} {1} / 3'.format(open_elec + 10, open_elec + 12)
        gas3 = '{0} {1} / {2}'.format(open_elec + 12, open_elec + 12, nb_virtual//2)
    else:
        if nb_elec <= 10:
            elec_in_gas = nb_elec
            gas1 = '2 4 / 2 !1s2s' 
            gas2 = '{0} {1} / 3 !2p'.format(nb_elec - 2, nb_elec)
            gas3 = '{0} {1} / {2}'.format(nb_elec, nb_elec, nb_virtual//2)
        else:
            elec_in_gas = 2 + open_elec
            gas1 = '0 2 / 1' 
            gas2 = '{0} {1} / 3'.format(open_elec, open_elec + 2)
            gas3 = '{0} {1} / {2}'.format(open_elec + 2, open_elec + 2, nb_virtual//2)

    # how many roots should we compute
    if open_elec == 1:
        if nb_openshell_orbit == 1:
            # s-block elements
            ciroot1 = '1 1'
            ciroot2 = '2 1'
        else:

            ciroot1 = '3 3'
            ciroot2 = '4 3'
    elif open_elec == 2:
        # # IV A elements
        # ciroot1 = '1 1'
        # ciroot2 = '2 1'
        # IV A elements
        ciroot1 = '1 9'
        ciroot2 = '2 6'
    elif open_elec == 3:
        # # V A elements
        # ciroot1 = '3 2'
        # ciroot2 = '4 2'
        # V A elements
        ciroot1 = '3 10'
        ciroot2 = '4 10'
    elif open_elec == 4:
        # # VI A elements
        # ciroot1 = '1 2'
        # ciroot2 = '2 2'
        # VI A elements
        ciroot1 = '1 9'
        ciroot2 = '2 6'
    elif open_elec == 5:
        # VII A elements
        ciroot1 = '3 3'
        ciroot2 = '4 3'
    elif open_elec == 0:
        # closed-shell elementss
        ciroot1 = '1 1'
        ciroot2 = '2 1'
    else:
        raise RuntimeError('This script only support main group elements')


    inactivate = (nb_elec - elec_in_gas) // 2
    krci = ('*KRCICALC', ['.CI PROGRAM', 'LUCIAREL', '.INACTIVE', str(inactivate),
                          '.GAS SHELLS', '3', gas1, gas2, gas3, '.CIROOTS', ciroot1, '.CIROOTS', ciroot2,
                          '.MAX CI', '120', '.MXCIVE', '60', '.ANALYZ', '.RSTRCI', 'rstr', '.CHECKP'])


    inp = Inpobj.from_file(fin)
    inp.add_keywords(krci)
    inp.write_to_file(fout)


if __name__ == '__main__':
    import sys

    argv = sys.argv[1:]
    for filename in argv:
        create_mrci_inp(fin=filename)
