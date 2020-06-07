#!/usr/bin/env python

"""
Function: Create input file from the previous output file, e.g., SCF calculation.
Author: Yingxing Cheng
Date: 10/21/2019
"""
from pydirac.core.atom import Atom
from pydirac.input.jobs import Inpobj


__all__ = ['get_mrci_inp']

def get_info_for_s_block(atom, v_min_e, v_max_e, calc_type='quadrupole'):
    """
    atom_info: information about specified atom
    Return:
         nb_active_elec: (int) the number of electrons in GAS
         gas_list: (list) a list for GAS setup
         root_list: (list) root for different symmetry
    """
    assert atom.info.block == 's'

    # obtain the number of electrons
    nb_closed_elec = atom.closed_elec()
    nb_open_elec = atom.openshell_elec()
    nb_total_elec = atom.total_nb_elec()

    # obtain orbit info
    ## res_orbits = [nb_closed_shells, nb_open_shells, nb_virtual_shells]
    orbitals_dirac_info = atom.orbit_count(res_virtual=True, v_min_e=v_min_e, v_max_e=v_max_e)
    nb_closed_shell = orbitals_dirac_info[0]
    nb_open_shell = orbitals_dirac_info[1]
    nb_vir_shell = orbitals_dirac_info[-1]

    assert (nb_closed_elec == nb_closed_shell)
    assert (nb_open_elec <= nb_open_shell)
    period = atom.info.period

    gas_list = []
    root_list = []

    if nb_open_shell == 0:
        # for nobel gas elements
        if period == 1:
            # He
            nb_active_elec = 2
            gas1 = '0 2 / 1'
            gas2 = '2 2 / {0}'.format(nb_vir_shell // 2)
            gas_list.extend([gas1, gas2])
        elif period == 2:
            # Be
            nb_active_elec = 4
            gas1 = '2 4 / 2'
            gas2 = '4 4 / {0}'.format(nb_vir_shell // 2)
            gas_list.extend([gas1, gas2])
        elif period == 3:
            # Mg
            nb_active_elec = 10
            gas1 = '0 2 / 1'
            gas2 = '8 10 / 4'
            gas3 = '10 10 / {0}'.format(nb_vir_shell // 2)
            gas_list.extend([gas1, gas2, gas3])
            pass
        else:
            # Ca Sr Ba Ra
            # including d electrons
            nb_active_elec = 20
            gas1 = '10 12 / 6 ! (n-1)s2 (n-2)d10'
            gas2 = '18 20 / 4 ! (n-1)p6 (n)s2'
            gas3 = '20 20 / {0} ! virtual orbitals'.format(nb_vir_shell // 2)
            gas_list.extend([gas1, gas2, gas3])
    else:
        # for s-block open-shell elements
        if period > 3:
            # K Rb Cs Fr
            # including d electrons
            nb_active_elec = 19
            gas1 = '10 12 / 6 ! (n-1)s2 (n-2)d10'
            gas2 = '17 19 / 4 ! (n-1)p6 (n)s1'
            gas3 = '19 19 / {0} ! virtual orbitals'.format(nb_vir_shell // 2)
            gas_list.extend([gas1, gas2, gas3])
        elif period == 3:
            # for Na
            nb_active_elec = 9
            gas1 = '0 2 / 1'
            gas2 = '7 9 / 4'
            gas3 = '9 9 / {0}'.format(nb_vir_shell // 2)
            gas_list.extend([gas1, gas2, gas3])
        elif period == 2:
            # for Li
            nb_active_elec = 3
            gas1 = '1 3 / 2 ! 1s'
            gas2 = '3 3 / {0}'.format(nb_vir_shell // 2)
            gas_list.extend([gas1, gas2])
        else:
            # H
            assert period == 1
            nb_active_elec = 1
            gas1 = '0 1 / 1 ! 1s'
            gas2 = '1 1 / {0}'.format(nb_vir_shell // 2)
            gas_list.extend([gas1, gas2])


    # specify many roots should we compute
    if nb_open_elec == 1:
        # IA: H Li Na K Rb Cs Fr
        ciroot1 = '3 1'
        root_list.extend([ciroot1])
    elif nb_open_elec == 0:
        # IIA: He Be Mg Ca Sr Ba Ra
        ciroot1 = '1 1'
        # ciroot2 = '2 1'
        root_list.extend([ciroot1])
    else:
        raise RuntimeError('The number of electrons in open shell '
                           'is wrong: {0}'.format(nb_open_elec))

    return nb_active_elec, gas_list, root_list


def get_info_for_p_block(atom, v_min_e, v_max_e, calc_type='quadrupole'):
    """
    atom_info: information about specified atom
    Return:
         nb_active_elec: (int) the number of electrons in GAS
         gas_list: (list) a list for GAS setup
         root_list: (list) root for different symmetry
    """
    assert atom.info.block == 'p'

    # obtain the number of electrons
    nb_closed_elec = atom.closed_elec()
    nb_open_elec = atom.openshell_elec()
    nb_total_elec = atom.total_nb_elec()

    # obtain orbit info
    ## res_orbits = [nb_closed_shells, nb_open_shells, nb_virtual_shells]
    orbitals_dirac_info = atom.orbit_count(res_virtual=True, v_min_e=v_min_e, v_max_e=v_max_e)
    nb_closed_shell = orbitals_dirac_info[0]
    nb_open_shell = orbitals_dirac_info[1]
    nb_vir_shell = orbitals_dirac_info[-1]

    assert (nb_closed_elec == nb_closed_shell)
    assert (nb_open_elec <= nb_open_shell)
    period = atom.info.period

    gas_list = []
    root_list = []

    if nb_open_shell == 0:
        # for nobel gas elements
        # if period == 1:
        #     # He
        #     nb_active_elec = 2
        #     gas1 = '0 2 / 1'
        #     gas2 = '2 2 / {0}'.format(nb_vir_shell // 2)
        #     gas_list.extend([gas1, gas2])
        if period == 2:
            # Ne
            nb_active_elec = 10
            gas1 = '2 4 / 2'
            gas2 = '8 10 / 3'
            gas3 = '10 10 / {0}'.format(nb_vir_shell // 2)
            gas_list.extend([gas1, gas2, gas3])
        else:
            # Ar, Kr, Xe, Rn, Og
            # including d electrons
            nb_active_elec = 18
            gas1 = '10 12 / 6 ! (n-1)s2 (n-2)d10'
            gas2 = '16 18 / 3 ! (n-1)p6'
            gas3 = '18 18 / {0} ! virtual orbitals'.format(nb_vir_shell // 2)
            gas_list.extend([gas1, gas2, gas3])
    else:
        # for p-block open-shell elements
        if period > 3:
            # including d electrons
            nb_active_elec = 10 + 2 + nb_open_elec
            gas1 = '10 12 / 6 ! (n-1)s2 (n-2)d10'
            gas2 = '{0} {1} / 3 ! (n-1)p6'.format(nb_open_elec + 10,
                                                  nb_open_elec + 12)
            gas3 = '{0} {1} / {2} ! virtual orbitals'.format(nb_open_elec + 12,
                                                             nb_open_elec + 12,
                                                             nb_vir_shell // 2)
            gas_list.extend([gas1, gas2, gas3])
        elif period == 3:
            # for Al Si P S Cl
            nb_active_elec = 2 + nb_open_elec
            gas1 = '0 2 / 1'
            gas2 = '{0} {1} / 3'.format(nb_open_elec,
                                        nb_open_elec + 2)
            gas3 = '{0} {1} / {2}'.format(nb_open_elec + 2,
                                          nb_open_elec + 2,
                                          nb_vir_shell // 2)
            gas_list.extend([gas1, gas2, gas3])
        else:
            # without d electrons
            assert (period == 2)
            # for B C N O F
            nb_active_elec = nb_total_elec
            gas1 = '2 4 / 2 ! 1s2s'
            gas2 = '{0} {1} / 3 ! 2p'.format(nb_total_elec - 2, nb_total_elec)
            gas3 = '{0} {1} / {2}'.format(nb_total_elec, nb_total_elec, nb_vir_shell // 2)
            gas_list.extend([gas1, gas2, gas3])


    # specify many roots should we compute
    if nb_open_elec == 1:
        # IIIA: B Al Ga In Tl Nh
        ciroot1 = '3 3'
        ciroot2 = '4 3'
        root_list.extend([ciroot1, ciroot2])
    elif nb_open_elec == 2:
        # IVA: C Si Ge Sn Pb Fl
        if calc_type == 'dipole':
            ciroot1 = '1 1'
            # ciroot2 = '2 1'
            root_list.extend([ciroot1])
        elif calc_type == 'quadrupole':
            ciroot1 = '1 9'
            ciroot2 = '2 6'
            root_list.extend([ciroot1, ciroot2])
        else:
            raise TypeError('Calc_type wrong! "dipole" or "quadrupole"')
    elif nb_open_elec == 3:
        # VA: N P As Sb Bi Mc
        if calc_type == 'dipole':
            ciroot1 = '3 2'
            ciroot2 = '4 2'
            root_list.extend([ciroot1, ciroot2])
        elif calc_type == 'quadrupole':
            ciroot1 = '3 10'
            # ciroot2 = '4 10'
            root_list.extend([ciroot1])
        else:
            raise TypeError('Calc_type wrong! "dipole" or "quadrupole"')
    elif nb_open_elec == 4:
        # # VIA: O S Se Te Po Lv
        if calc_type == 'dipole':
            ciroot1 = '1 3'
            ciroot2 = '2 2'
            root_list.extend([ciroot1, ciroot2])
        elif calc_type == 'quadrupole':
            ciroot1 = '1 9'
            ciroot2 = '2 6'
            root_list.extend([ciroot1, ciroot2])
        else:
            raise TypeError('Calc_type wrong! "dipole" or "quadrupole"')
    elif nb_open_elec == 5:
        # VIIA: F Cl Br I At Ts
        ciroot1 = '3 3'
        ciroot2 = '4 3'
        root_list.extend([ciroot1, ciroot2])
    elif nb_open_elec == 0:
        # closed-shell elementss
        ciroot1 = '1 1'
        # ciroot2 = '2 1'
        root_list.extend([ciroot1])
    else:
        raise RuntimeError('The number of electrons in open shell '
                           'is wrong: {0}'.format(nb_open_elec))

    return nb_active_elec, gas_list, root_list


def get_info_for_d_block(atom, v_min_e, v_max_e, calc_type='quadrupole'):
    """
    atom_info: information about specified atom
    Return:
         nb_active_elec: (int) the number of electrons in GAS
         gas_list: (list) a list for GAS setup
         root_list: (list) root for different symmetry
    """
    assert atom.info.block == 'd'

    # obtain the number of electrons
    nb_closed_elec = atom.closed_elec()
    nb_open_elec = atom.openshell_elec()
    nb_total_elec = atom.total_nb_elec()

    # obtain orbit info
    ## res_orbits = [nb_closed_shells, nb_open_shells, nb_virtual_shells]
    orbitals_dirac_info = atom.orbit_count(res_virtual=True, v_min_e=v_min_e, v_max_e=v_max_e)
    nb_closed_shell = orbitals_dirac_info[0]
    nb_open_shell = orbitals_dirac_info[1]
    nb_vir_shell = orbitals_dirac_info[-1]

    assert (nb_closed_elec == nb_closed_shell)
    assert (nb_open_elec <= nb_open_shell)
    period = atom.info.period

    gas_list = []
    root_list = []

    if nb_open_shell == 0:
        # IIB
        # Zn Cd Hg Cn
        # including d electrons
        nb_active_elec = 12
        gas1 = '8 10 / 5 ! (n-1)d10'
        gas2 = '10 12 / 1 ! (n)s2'
        gas3 = '12 12 / {0} ! virtual orbitals'.format(nb_vir_shell // 2)
        gas_list.extend([gas1, gas2, gas3])
    else:
        # for d-block open-shell elements
        if nb_open_shell == 2:
            # s is open
            # Cu Ag Au Rg
            assert (nb_open_elec == 1)
            nb_active_elec = 11
            gas1 = '8 10 / 5 ! (n-1)d10'
            gas2 = '9 11 / 1 ! (n)s1'
            gas3 = '11 11 / {0} ! virtual orbitals'.format(nb_vir_shell // 2)
            gas_list.extend([gas1, gas2, gas3])
        elif nb_open_shell == 10:
            # d is open
            nb_active_elec = 2 + nb_open_elec
            # normally there is no this case, since in DIRAC (n-1)d locates
            # before (n)s shells
            gas1 = '0 2 / 1 ! (n)s2'
            gas2 = '{0} {1} / 5 ! (n-1)d10'.format(nb_open_elec, nb_open_elec + 2)
            gas3 = '{0} {0} / {1} ! virtual orbitals'.format(nb_open_elec+2,
                                                             nb_vir_shell // 2)
            gas_list.extend([gas1, gas2, gas3])
        elif nb_open_shell == 12:
            # s + d are open
            nb_active_elec = nb_open_elec
            gas1 = '{0} {1} / 6 ! (n-1)d10 (n)s2'.format(nb_open_elec-2,
                                                         nb_open_elec)
            gas2 = '{0} {0} / {1} ! virtual orbitals'.format(nb_open_elec,
                                                             nb_vir_shell // 2)
            gas_list.extend([gas1, gas2])
        else:
            raise RuntimeError('Unknown number of open shell: {0}'.format(nb_open_shell))



    # specify many roots should we compute
    if nb_open_elec == 1:
        if nb_open_shell == 2:
            # only s is open
            ciroot1 = '1 1'
            # ciroot2 = '2 1'
            root_list.extend([ciroot1])
        elif nb_open_shell == 12:
            raise NotImplementedError('For two open-shell d elements, '
                                      'there is no implement!')
        else:
            raise RuntimeError('the number of open shell {0} '
                               'is error!'.format(nb_open_shell))
    elif nb_open_elec == 0:
        # closed-shell elementss
        ciroot1 = '1 1'
        # ciroot2 = '2 1'
        root_list.extend([ciroot1])
    else:
        raise NotImplementedError('For open-shell d elements and more than '
                                  'one open-shell element, it has not '
                                  'been implmented yet!')

    return nb_active_elec, gas_list, root_list


def get_info_for_f_block(atom, v_min_e, v_max_e, calc_type='quadrupole'):
    """
    atom_info: information about specified atom
    Return:
         nb_active_elec: (int) the number of electrons in GAS
         gas_list: (list) a list for GAS setup
         root_list: (list) root for different symmetry
    """
    assert atom.info.block == 'f'

    # obtain the number of electrons
    nb_closed_elec = atom.closed_elec()
    nb_open_elec = atom.openshell_elec()
    nb_total_elec = atom.total_nb_elec()

    # obtain orbit info
    ## res_orbits = [nb_closed_shells, nb_open_shells, nb_virtual_shells]
    orbitals_dirac_info = atom.orbit_count(res_virtual=True, v_min_e=v_min_e, v_max_e=v_max_e)
    nb_closed_shell = orbitals_dirac_info[0]
    nb_open_shell = orbitals_dirac_info[1]
    nb_vir_shell = orbitals_dirac_info[-1]

    assert (nb_closed_elec == nb_closed_shell)
    assert (nb_open_elec <= nb_open_shell)
    period = atom.info.period

    gas_list = []
    root_list = []

    raise NotImplementedError('For open-shell f elements, it has not '
                              'been implmented yet!')



def get_mrci_inp(filename_input, filename_out='PYDIRAC.inp', v_min_e=-1, v_max_e=10.0):
    """
    Create MRCI input based on the previous SCF calculation.
    """

    atom = Atom.from_file(filename_input)

    # this is for p block elements
    if atom.info.block == 'p':
        nb_active_elec, gas_list, root_list = \
            get_info_for_p_block(atom, v_min_e=v_min_e, v_max_e=v_max_e, calc_type='quadrupole')
    elif atom.info.block == 's':
        nb_active_elec, gas_list, root_list = \
            get_info_for_s_block(atom, v_min_e=v_min_e, v_max_e=v_max_e, calc_type='quadrupole')
    elif atom.info.block == 'd':
        nb_active_elec, gas_list, root_list = \
            get_info_for_d_block(atom, v_min_e=v_min_e, v_max_e=v_max_e, calc_type='quadrupole')
    elif atom.info.block == 'f':
        nb_active_elec, gas_list, root_list = \
            get_info_for_f_block(atom, v_min_e=v_min_e, v_max_e=v_max_e, calc_type='quadrupole')
    else:
        raise RuntimeError('Unknown block element!')

    total_nb_elec = atom.total_nb_elec()
    print('total number electrons is {0}'.format(total_nb_elec))
    inactivate = (atom.total_nb_elec() - nb_active_elec) // 2
    nb_gas_shell = len(gas_list)
    nb_root_sym = len(root_list)

    krci_setup = ['.CI PROGRAM', 'LUCIAREL', '.INACTIVE', str(inactivate), '.GAS SHELLS']
    krci_setup.append(str(nb_gas_shell))
    krci_setup.extend(gas_list)
    for i in range(nb_root_sym):
        krci_setup.append('.CIROOTS')
        krci_setup.append(root_list[i])
    krci_setup.extend(['.MAX CI', '120', '.MXCIVE', '60', '.ANALYZ',
                       '.RSTRCI', 'rstr', '.CHECKP', '.NOOCCN'])

    krci = ('*KRCICALC', krci_setup)

    inp = Inpobj.from_file(filename_input)
    inp.add_keywords(krci)
    inp.write_to_file(filename_out)


if __name__ == '__main__':
    import sys

    argv = sys.argv[1:]
    for filename in argv:
        get_mrci_inp(filename_input=filename)
