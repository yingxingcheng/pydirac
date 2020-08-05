#!/usr/bin/env python

"""
Function: this file is default configuration for generate_inputs.py
Authour: Yingxing Cheng
Email: Yingxing.cheng@ugent.be
Date: 2019-10-18
"""
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

import os
from pydirac.core.periodic_table import special_elements, Element



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

    def get_shell(elec_config, NSYM=1):
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
        Element(nb_elec).symbol) + 'and Econf: {0}'.format(elec_config))
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


def get_d_DOSSSS_SCF_str(elec_conf, atom_type, *arg, **kwarg):
    context = """**DIRAC
.TITLE
 {atom_type}, DOSSSS, SCF
.WAVE F
.ANALYZE
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
**HAMILTONIAN
.DOSSSS
.OPERATOR
 ZDIPLEN
 COMFACTOR
 zff
**INTEGRALS
*READINP
.UNCONTRACT
**WAVE FUNCTIONS
.SCF
.RESOLVE
*SCF
{elec_conf}
.EVCCNV
1.0D-9  5.0D-8
.MAXITR
60
*END OF""".format(**{'elec_conf': elec_conf, 'atom_type': atom_type})
    return context


def get_q_DOSSSS_SCF_str(elec_conf, atom_type, *arg, **kwarg):
    context = """**DIRAC
.TITLE
 {atom_type}, DOSSSS, SCF
.WAVE F
.ANALYZE
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
**HAMILTONIAN
.DOSSSS
.OPERATOR
 'Theta quadru-field'
 DIAGONAL
 ZZTHETA
 COMFACTOR
 zff
**INTEGRALS
*READINP
.UNCONTRACT
**WAVE FUNCTIONS
.SCF
.RESOLVE
*SCF
{elec_conf}
.EVCCNV
1.0D-9  5.0D-8
.MAXITR
60
*END OF""".format(**{'elec_conf': elec_conf, 'atom_type': atom_type})
    return context


def get_d_DOSSSS_RELCCSD_str(elec_conf, atom_type, *arg, **kwarg):
    context = """**DIRAC
.TITLE
 {atom_type}, DC, RELCC
.WAVE F
.4INDEX
.ANALYZE
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
**HAMILTONIAN
.DOSSSS
.OPERATOR
 'Theta quadru-field'
 DIAGONAL
 ZZTHETA
 COMFACTOR
 zff
**INTEGRALS
*READINP
.UNCONTRACT
**WAVE FUNCTIONS
.SCF
.RESOLVE
.RELCCSD
*SCF
{elec_conf}
.EVCCNV
1.0D-9  5.0D-8
.MAXITR
60
**MOLTRA
.ACTIVE
-20 2 0.1
**RELCC
.ENERGY
.PRINT
1
*CCENER
.MAXIT
60
.NTOL
14
.NOSDT
*END OF""".format(**{'elec_conf': elec_conf, 'atom_type': atom_type})
    return context


def get_q_DOSSSS_RELCCSD_str(elec_conf, atom_type, *arg, **kwarg):
    context = """**DIRAC
.TITLE
 {atom_type}, DC, RELCC
.WAVE F
.4INDEX
.ANALYZE
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
**HAMILTONIAN
.DOSSSS
.OPERATOR
 'Theta quadru-field'
 DIAGONAL
 ZZTHETA
 COMFACTOR
 zff
**INTEGRALS
*READINP
.UNCONTRACT
**WAVE FUNCTIONS
.SCF
.RESOLVE
.RELCCSD
*SCF
{elec_conf}
.EVCCNV
1.0D-9  5.0D-8
.MAXITR
60
**MOLTRA
.ACTIVE
-20 2 0.1
**RELCC
.ENERGY
.PRINT
1
*CCENER
.MAXIT
60
.NTOL
14
.NOSDT
*END OF""".format(**{'elec_conf': elec_conf, 'atom_type': atom_type})
    return context


def get_q_DOSSSS_SCF_property_str(elec_conf, atom_type, *arg, **kwarg):
    context = """**DIRAC
.TITLE
 {atom_type}, DOSSSS, SCF
.WAVE F
.ANALYZE
.PROPERTIES
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
**HAMILTONIAN
.DOSSSS
.OPERATOR
 'Theta quadru-field'
 DIAGONAL
 ZZTHETA
 COMFACTOR
 zff
**INTEGRALS
*READINP
.UNCONTRACT
**WAVE FUNCTIONS
.SCF
.RESOLVE
*SCF
{elec_conf}
.EVCCNV
1.0D-9  5.0D-8
.MAXITR
90
**PROPERTIES
.QUADRUPOLE
*END OF""".format(**{'elec_conf': elec_conf, 'atom_type': atom_type})
    return context


def get_d_X2C_SCF_str(elec_conf, atom_type, *arg, **kwarg):
    context = """**DIRAC
.TITLE
 {atom_type}, X2C, SCF
.WAVE F
.ANALYZE
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
**HAMILTONIAN
.X2C
.OPERATOR
 ZDIPLEN
 COMFACTOR
 zff
**INTEGRALS
*READINP
.UNCONTRACT
**WAVE FUNCTIONS
.SCF
.RESOLVE
*SCF
{elec_conf}
.EVCCNV
1.0D-9  5.0D-8
.MAXITR
90
*END OF""".format(**{'elec_conf': elec_conf, 'atom_type': atom_type})
    return context


def get_q_X2C_SCF_str(elec_conf, atom_type, *arg, **kwarg):
    context = """**DIRAC
.TITLE
 {atom_type}, X2C, SCF
.WAVE F
.ANALYZE
.PROPERTIES
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
**HAMILTONIAN
.X2C
.OPERATOR
 'Theta quadru-field'
 DIAGONAL
 ZZTHETA
 COMFACTOR
 zff
**INTEGRALS
*READINP
.UNCONTRACT
**WAVE FUNCTIONS
.SCF
.RESOLVE
*SCF
{elec_conf}
.EVCCNV
1.0D-9  5.0D-8
.MAXITR
60
**PROPERTIES
.QUADRUPOLE
*END OF""".format(**{'elec_conf': elec_conf, 'atom_type': atom_type})
    return context


def get_d_X2C_NOSPIN_SCF_str(elec_conf, atom_type, *arg, **kwarg):
    context = """**DIRAC
.TITLE
 {atom_type}, scalar relativistic, SCF
.WAVE F
.ANALYZE
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
**HAMILTONIAN
.X2C
.NOSPIN
.OPERATOR
 ZDIPLEN
 COMFACTOR
 zff
**INTEGRALS
*READINP
.UNCONTRACT
**WAVE FUNCTIONS
.SCF
.RESOLVE
*SCF
{elec_conf}
.EVCCNV
1.0D-9  5.0D-8
.MAXITR
90
*END OF""".format(**{'elec_conf': elec_conf, 'atom_type': atom_type})
    return context


def get_q_X2C_NOSPIN_SCF_str(elec_conf, atom_type, *arg, **kwarg):
    context = """**DIRAC
.TITLE
 {atom_type}, X2C_NOSPIN, SCF
.WAVE F
.ANALYZE
.PROPERTIES
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
**HAMILTONIAN
.X2C
.NOSPIN
.OPERATOR
 'Theta quadru-field'
 DIAGONAL
 ZZTHETA
 COMFACTOR
 zff
**INTEGRALS
*READINP
.UNCONTRACT
**WAVE FUNCTIONS
.SCF
.RESOLVE
*SCF
{elec_conf}
.EVCCNV
1.0D-9  5.0D-8
.MAXITR
60
**PROPERTIES
.QUADRUPOLE
*END OF""".format(**{'elec_conf': elec_conf, 'atom_type': atom_type})
    return context


def get_d_X2C_NOSPIN_RELCCSD_str(elec_conf, atom_type, *arg, **kwarg):
    context = """**DIRAC
.TITLE
 {atom_type}, scalar relativistic, RELCC
.WAVE F
.4INDEX
.ANALYZE
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
**HAMILTONIAN
.X2C
.NOSPIN
.OPERATOR
 ZDIPLEN
 COMFACTOR
 zff
**INTEGRALS
*READINP
.UNCONTRACT
**WAVE FUNCTIONS
.SCF
.RESOLVE
.RELCCSD
*SCF
{elec_conf}
.EVCCNV
1.0D-9  5.0D-8
.MAXITR
60
**RELCC
.ENERGY
.PRINT
1
*CCENER
.MAXIT
60
.NTOL
14
.NOSDT
*END OF""".format(**{'elec_conf': elec_conf, 'atom_type': atom_type})
    return context


def get_q_X2C_NOSPIN_RELCCSD_str(elec_conf, atom_type, *arg, **kwarg):
    context = """**DIRAC
.TITLE
 {atom_type}, scalar relativistic, RELCC
.WAVE F
.4INDEX
.ANALYZE
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
**HAMILTONIAN
.X2C
.NOSPIN
.OPERATOR
 'Theta quadru-field'
 DIAGONAL
 ZZTHETA
 COMFACTOR
 zff
**INTEGRALS
*READINP
.UNCONTRACT
**WAVE FUNCTIONS
.SCF
.RESOLVE
.RELCCSD
*SCF
{elec_conf}
.EVCCNV
1.0D-9  5.0D-8
.MAXITR
60
**RELCC
.ENERGY
.PRINT
1
*CCENER
.MAXIT
60
.NTOL
14
.NOSDT
*END OF""".format(**{'elec_conf': elec_conf, 'atom_type': atom_type})
    return context


def get_mol_by_custom_basis(atom_type, atom_index, basis_choice, basis_info):
    template = """INTGRL
{atom} atom
{basis}
C   1    2  X  Y
    {elec}.    1
{atom}    0.0000000000        0.0000000000        0.0000000000
{basis_info}
""".format(**{'atom': atom_type, 'elec': atom_index, 'basis': basis_choice,
              'basis_info': basis_info})
    return template


def get_mol_by_default_basis(atom_type, atom_index, basis_type):
    template = """INTGRL
{atom} atom
{basis}
C   1    2  X  Y
    {elec}.    1
{atom}    0.0000000000        0.0000000000        0.0000000000
LARGE BASIS {basis}
FINISH

""".format(**{'atom': atom_type, 'elec': atom_index, 'basis': basis_type})
    return template


def get_script(atom_type, field_str, res_fname, calc_method,
               field_value, dir_path, reference_str, basis_type):
    script_tp = """#!/bin/bash

if [ -f {res_fname} ]; then
    #mv {res_fname} {res_fname}.$$
    {calc_fname} {res_fname} {f_begin} {f_end} {f_nb}
    exit
fi

for i in {field_value}; do
    pam-dirac --noarch  --replace zff=$i --mol={atom_type}_{basis_type}.mol --inp={atom_type}_{basis_type}.inp
    grep "@ Total {calc_method} energy" {atom_type}_{basis_type}_{atom_type}_{basis_type}_zff\\=$i.out >> {res_fname}
done

if [ '{reference}' == '' ]; then
    grep "@ Total {calc_method} energy" {atom_type}_{basis_type}_{atom_type}_{basis_type}_zff\\=0.0000.out >> {res_fname}
else
    pam-dirac --noarch  --replace zff={reference} --mol={atom_type}.mol --inp={atom_type}.inp
    grep "@ Total {calc_method} energy" {atom_type}_{basis_type}_{atom_type}_{basis_type}_zff\\={reference}.out >> {res_fname}
fi

{calc_fname} {res_fname} {f_begin} {f_end} {f_nb}
""".format(**{'atom_type': atom_type, 'field_value': field_str,
              'res_fname': res_fname, 'calc_method': calc_method,
              'f_begin': field_value[0], 'f_end': field_value[-1],
              'f_nb': len(field_value),
              'calc_fname': os.path.join(dir_path, 'calc_polarizability.py'),
              'reference': reference_str, 'basis_type': basis_type})
    return script_tp
