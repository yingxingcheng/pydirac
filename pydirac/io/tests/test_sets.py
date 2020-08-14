from unittest import TestCase
import os

from pathlib import Path
from pydirac.io.sets import *
from pydirac.io.outputs import Output


MOUDLE_DIR = Path(__file__).resolve().parent
DATA_DIR = str(MOUDLE_DIR.parent.parent / 'unit_tests'/'data')

class TestAtomicCCSet(TestCase):
    def setUp(self) -> None:
        self.cc_set = AtomicCCSet(molecule=Molecule(['Mc'], [[0., 0., 0.]]))

    def test_inp(self):
        inp_from_code = self.cc_set.inp
        inp_from_dict = {
            'DIRAC': {
                'TITLE': 'DIRAC input file generated by PYDIRAC',
                'WAVE F': True,
                'ANALYZE': True,
                '4INDEX': True,
            },
            'ANALYZE': {
                'MULPOP': {
                    '_en': True,
                    'VECPOP': '1..oo',
                }
            },
            'HAMILTONIAN': {
                'DOSSSS': True,
                'X2C': False,
                'NOSPIN': False,
            },
            'INTEGRALS': {
                'READINP': 'UNCONTRACT',
            },
            'WAVE FUNCTIONS': {
                'RESOLVE': True,
                'SCF': {
                    '_en': True,
                    'EVCCNV': '1.0D-9 5.0D-8',
                    'MAXITR': 60,
                    'CLOSED SHELL': '112',
                    'OPEN SHELL': '1\n3/6',
                },
                'RELCCSD': True,
            },
            'MOLTRA': {
                'ACTIVE': '-20.0 25.0 0.01',
            },
            'RELCC': {
                'ENERGY': True,
                'PRINT': 1,
                'CCENER': True,
                'MAXIT': 60,
                'NTOL': 10,
                'NOSDT': False,
            }

        }
        self.assertDictEqual(inp_from_code, inp_from_dict)

    def test_mol(self):
        mol_from_code = {
            'basis_type': self.cc_set.mol.basis_type,
            'basis_lib': self.cc_set.mol.basis_lib
        }
        mol_from_dict = {
            'basis_type': 'dyall.v3z',
            'basis_lib': 'BASIS',
        }
        self.assertDictEqual(mol_from_code, mol_from_dict)

    def test_from_prev_calc(self):
        output = str(Path(
            MOUDLE_DIR) / '../../unit_tests/data/Rn_q_so/d-aug-dyall.acv3z_+0.00001')
        #print(output)
        new_cc = AtomicCCSet.from_prev_calc(output, no_scf=True,
                                            user_mol_settings={
                                                'basis_type': 's-aug-ANO-RCC'})
        # new_cc = AtomicCCSet(molecule=Molecule([7], [[0.0, 0.0,0.0]]),
        #                      hamiltonian_mode='2C',is_spinfree=True ,is_ff=True, ff_mode='D',
        #                      user_mol_settings={'basis_type':'s-aug-ANO-RCC'})

        print(new_cc.inp)
        print(new_cc.mol)
        data_dir = os.path.abspath(
            os.path.join(MOUDLE_DIR, '../unit_tests/', 'data'))
        new_cc.write_input(output_dir=os.path.join(data_dir, 'dirac_input'),
                           make_dir_if_not_present=True)


class TestAtomicCISet(TestCase):

    def test_inp(self):
        self.fail()

    def test_mol(self):
        self.fail()

    def test_from_prev_dhf_calc(self):
        fname = os.path.join(DATA_DIR, 'Cu_DHF.out')
        out = Output(fname)
        atom = out.parse_orbit()
        #print(atom.as_dict())
        #print(atom.nb_closed_elec(), atom.nb_open_elec())
        #print(atom.as_dict())

        # fout = os.path.join(DATA_DIR, 'PYDIRAC.inp')
        ci_set = AtomicCISet.from_prev_dhf_calc(fname)
        print(ci_set.inp)
        print(ci_set.mol)
