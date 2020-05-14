#!/usr/bin/env python

"""
Function: this file is default configuration for generate_inputs.py
Authour: Yingxing Cheng
Email: Yingxing.cheng@ugent.be
Date: 2019-10-18
"""
import os

periodic_table=['dummy', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F',
           'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
           'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu','Zn',
           'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb','Sr','Y','Zr','Nb','Mo',
           'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In','Sn','Sb','Te','I','Xe',
           'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm','Sm','Eu','Gd','Tb',
           'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf','Ta','W','Re','Os','Ir',
           'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po','At','Rn','Fr','Ra',
           'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am','Cm','Bk','Cf','Es','Fm',
           'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh','Hs','Mt','Ds','Rg',
           'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

special_elements = {
        24: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^1 3d^5',
        29: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^1 3d^10',
        41: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^4',
        42: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^5',
        44: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^7',
        45: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^8',
        46: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 4s^10',
        47: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10',
        57: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10 5p^6 6s^2 5d^1',
        58: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10 5p^6 6s^2 4f^1 5d^1',
        64: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10 5p^6 6s^2 4f^7 5d^1',
        64: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10 5p^6 6s^2 4f^7 5d^1',
        78: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^1 4f^14 5d^9',
        79: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^1 4f^14 5d^10',
        79: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^1 4f^14 5d^10',
        89: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 6d^1',
        90: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 6d^2',
        91: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^2 6d^1',
        92: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^3 6d^1',
        93: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^4 6d^1',
        96: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^7 6d^1',
        103: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^14 7p^1',
        110: '1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^1 5f^14 6d^9',
    }


def get_element_by_id(id):
    """
    Get element by given index.
    """
    if type(id) == int:
        return periodic_table[id]
    else:
        raise RuntimeError('Type id should be int')


def get_element_id_by_type(e_type):
    """
    Get element number by given element type.
    """

    idx = periodic_table.index(e_type)
    if idx >= 0:
        return idx
    else:
        raise RuntimeError('Wrong element type!')


def get_d_DOSSSS_SCF(elec_conf, atom_type, *arg, **kwarg):
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
90
*END OF""".format(**{'elec_conf': elec_conf, 'atom_type': atom_type})
    return context


def get_d_X2C_SCF(elec_conf, atom_type, *arg, **kwarg):
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


def get_d_X2C_NOSPIN_SCF(elec_conf, atom_type, *arg, **kwarg):
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

def get_q_DOSSSS_SCF_p(elec_conf, atom_type, *arg, **kwarg):
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

def get_q_DOSSSS_SCF(elec_conf, atom_type, *arg, **kwarg):
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

def get_q_X2C_SCF(elec_conf, atom_type, *arg, **kwarg):
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

def get_q_X2C_NOSPIN_SCF(elec_conf, atom_type, *arg, **kwarg):
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


def get_d_X2C_NOSPIN_RELCCSD(elec_conf, atom_type, *arg, **kwarg):
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

def get_d_DOSSSS_RELCCSD(elec_conf, atom_type, *arg, **kwarg):
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
