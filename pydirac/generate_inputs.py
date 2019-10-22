#!/usr/bin/env python

import os
from config import *

dir_path = os.path.dirname(os.path.abspath(__file__))


def get_elec_config(n):
    """Get electron configuration by specifying the number of electrons.
    """
    if n in special_elements.keys():
        return special_elements[n]

    rule = '1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p'
    nb_dict = {'s': 2, 'p': 6, 'd': 10, 'f': 14}
    orbitals = [(i, nb_dict[i[-1]]) for i in rule.split()]
    output = []

    for orbital, size in orbitals:
        k = min(size, n)
        output.append("%s^%d" % (orbital, k))
        n -= k
        if n <= 0:
            break

    orbital_info = ' '.join(output)
    return orbital_info


def get_shell(orbital_info, NSYM=1):
    """
    """
    rule = '1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p'
    nb_dict = {'s': 2, 'p': 6, 'd': 10, 'f': 14}
    cs_even, cs_odd = 0, 0
    os_even, os_odd = [], []
    for i in elec_config.split():
        # 1s^2
        t, nb = i.split('^')
        if int(nb) < nb_dict[t[-1]]:
            # open_shell
            if t[-1] in ['s', 'd']:
                # even
                os_even.append((t, int(nb)))
            else:
                os_odd.append((t, int(nb)))
        else:
            if t[-1] in ['s', 'd']:
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
            cs_str.append('{0}'.format(cs_total))
        if not len(os_total):
            os_str = None
        else:
            os_str.append('{0}'.format(nb_os))
            for i in os_total:
                os_str.append('{0}/{1}'.format(i[-1], nb_dict[i[0][1]]))
    else:
        # for closed shell: 10, 10
        # for open shell:
        # 2
        # 5/10,0
        # 1/2,0
        if cs_total:
            cs_str.append('{0} {1}\n'.format(cs_even, cs_odd))
        else:
            cs_str = None

        if nb_os:
            os_str.append('{0}\n'.format(nb_os))
            for i in os_total:
                if i[0][-1] in ['s', 'd']:
                    os_str.append(
                        '{0}/{1},{2}\n'.format(i[-1], nb_dict[i[0][1]], 0))
                else:
                    os_str.append(
                        '{0}/{1},{2}\n'.format(i[-1], 0, nb_dict[i[0][1]]))
        else:
            os_str = None

    if cs_str:
        cs_str = '\n'.join(cs_str)
    if os_str:
        os_str = '\n'.join(os_str)
    return cs_str, os_str


def get_dirac_shell_str(nb_elec, f=None):
    """Get electron configurations in DIRAC according to the number of electrons.

    For example:
    for NSYM = 1
    .CLOSED SHELL
    2
    .OPEN SHELL
    1
    1/2

    or for NSYM != 1
    .CLSED SHELL
    2 0
    .OPEN SHELL
    1
    1/2,0
    """

    def get_shell(orbital_info, NSYM=1):
        rule = '1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p'
        nb_dict = {'s': 2, 'p': 6, 'd': 10, 'f': 14}
        cs_even, cs_odd = 0, 0
        os_even, os_odd = [], []
        for i in elec_config.split():
            # 1s^2
            t, nb = i.split('^')
            if int(nb) < nb_dict[t[-1]]:
                # open_shell
                if t[-1] in ['s', 'd']:
                    # even
                    os_even.append((t, int(nb)))
                else:
                    os_odd.append((t, int(nb)))
            else:
                if t[-1] in ['s', 'd']:
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
                cs_str.append('{0}'.format(cs_total))
            if not len(os_total):
                os_str = None
            else:
                os_str.append('{0}'.format(nb_os))
                for i in os_total:
                    os_str.append('{0}/{1}'.format(i[-1], nb_dict[i[0][1]]))
        else:
            # for closed shell: 10, 10
            # for open shell:
            # 2
            # 5/10,0
            # 1/2,0
            if cs_total:
                cs_str.append('{0} {1}\n'.format(cs_even, cs_odd))
            else:
                cs_str = None

            if nb_os:
                os_str.append('{0}\n'.format(nb_os))
                for i in os_total:
                    if i[0][-1] in ['s', 'd']:
                        os_str.append(
                            '{0}/{1},{2}\n'.format(i[-1], nb_dict[i[0][1]], 0))
                    else:
                        os_str.append(
                            '{0}/{1},{2}\n'.format(i[-1], 0, nb_dict[i[0][1]]))
            else:
                os_str = None
        if cs_str:
            cs_str = '\n'.join(cs_str)
        if os_str:
            os_str = '\n'.join(os_str)
        return cs_str, os_str

    elec_config = get_elec_config(nb_elec)
    closed_shell, open_shell = get_shell(elec_config)
    res_str = []
    res_str.append('# {0} atom '.format(
        get_element_by_id(nb_elec)) + 'and Econf: {0}'.format(elec_config))
    if closed_shell:
        res_str.append(".CLOSED SHELL")
        res_str.append(closed_shell)
    if open_shell:
        res_str.append(".OPEN SHELL")
        res_str.append(open_shell)

    res_str = '\n'.join(res_str)
    if f:
        f.write(res_str)
    return res_str


def create_input_files(atom_type='He', basis_type='dyall.v2z', basis_choice='BASIS',
                       field_value=None, scratch_dir=None,
                       create_script=True, suffix=None):
    """Create input file according atom type and basis specified by users.
    """
    if scratch_dir is None:
        if suffix:
            scratch_dir = atom_type + suffix
        else:
            scratch_dir = atom_type
    if not os.path.exists(scratch_dir):
        os.makedirs(scratch_dir)

    # write molecuar file of DIRAC
    if atom_type is int:
        atom_index = atom_type
        atom_type = get_element_by_id(atom_index)
    else:
        atom_index = get_element_id_by_type(atom_type)

    if basis_choice == 'EXPLICIT':
        with open('basis/{0}.dat'.format(atom_type), 'r') as f:
            basis_info = f.read()
        template = get_mol_by_custom_basis(atom_type, atom_index,
                                           basis_choice, basis_info)
    else:
        template = get_mol_by_default_basis(atom_type, atom_index, basis_type)


    fname = atom_type + '_' + basis_type + '.mol'
    fname = os.path.join(scratch_dir, fname)
    with open(fname, 'w') as f:
        f.write(template)

    # create input file of DIRAC
    elec_conf = get_dirac_shell_str(atom_index)
    context = get_d_DOSSSS_SCF(elec_conf, atom_type)

    inp_fname = atom_type + '_' + basis_type + '.inp'
    #inp_fname = atom_type + '_' + 'Ml0.inp'
    inp_fname = os.path.join(scratch_dir, inp_fname)
    with open(inp_fname, 'w') as f:
        f.write(context)

    # create script to execute with different field values.
    if create_script:
        field_str = []
        if field_value is None:
            field_str.append('+0.0000')
        else:
            for f in field_value:
                if f > 0:
                    field_str.append('+{0:.4f}'.format(f))
                else:
                    field_str.append('{0:.4f}'.format(f))

        field_str = ' '.join(field_str)
        # TODO: support multiple methods (post-HF)
        reference_str = '' if 0.00 in field_value else '+0.0000'

        calc_method = 'CCSD(T)'
        res_fname = '_'.join(['res', atom_type,
                              calc_method.replace('(', '\(').
                              replace(')', '\)'),
                              basis_type]) + '.dat'
        script_tp = get_script(atom_type, field_str, res_fname, calc_method,
                               field_value, dir_path, reference_str, basis_type)

        spt_fname = 'run_{0}_{1}.sh'.format(atom_type, basis_type)
        spt_fname = os.path.join(scratch_dir, spt_fname)
        with open(spt_fname, 'w') as f:
            f.write(script_tp)

        #os.chmod(scratch_dir, 777)
        os.chmod(spt_fname, 0o777)


def basis_helper(filename = 'ANO-RCC'):
    with open(filename, 'r') as f:
        textlines = f.readlines()

    textlines = textlines[43:-1]
    res_list = []
    element_mess = []
    for line in textlines:
        first_word = line.strip().split()[0]
        if first_word == '!' and 'functions' not in line:
            if len(element_mess)>3:
                res_list.append(''.join(element_mess))
                element_mess = []
            element_mess.append(line)
            continue

        elif first_word == 'a':
            element_mess.append('!' + line)
            continue

        elif first_word == 'H':
            nb_exp = int(line.strip().split()[1])
            tmp_str = 'f{0:4d}{1:5d}\n'.format(nb_exp,0)
            element_mess.append(tmp_str)
            continue
        else:
            element_mess.append(line)

    res_list.append(''.join(element_mess))

    if not os.path.exists('basis'):
        os.makedirs('basis')

    assert(len(res_list) == 96)
    for i in range(1,97):
        name = get_element_by_id(i)
        with open('basis/{0}.dat'.format(name),'w') as f:
            nb_block = res_list[i-1].count('functions')
            f.write('LARGE EXPLICIT {0} {1}\n'.format(nb_block, '1 '*nb_block))
            f.write(res_list[i-1])
            f.write('FINISH\n')


    # with open('elec_config.dat', 'w') as f:
    #    for i in range(1, 119):
    #        get_dirac_shell_str(i, f)
            # print('"""', file=f)
    # get_dirac_shell_str(24)
    import numpy as np
    # for atom in ['Be', 'B', 'C', 'N']:
    #     create_input_files(atom_type=atom, field_value=np.linspace(
    #         0.0120, 0.002, 15), basis_type='cc-PVDZ')
    #

    #for atom in ['He']:
    #    create_input_files(atom_type=atom, field_value=np.linspace(
    #        0.0190, 0.000, 20), basis_type='cc-PVDZ')

    basis_helper(filename='dirac_basis')
