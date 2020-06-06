#!/usr/bin/env python
import unittest
from pydirac.dirac_input import Inpobj
from pydirac.get_orbit_info import Atom
from pydirac.create_input import create_mrci_inp
import mendeleev
import os


module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))


def test_write_file():
    """
    test write_file
    """
    fname = os.path.join(data_dir, 'Cu_DHF.out')
    atom = Atom.from_file(fname)
    print(atom.closed_elec(), atom.openshell_elec())
    print(atom.info.symbol)
    print(atom.info.period)
    if atom.info.group:
        print(atom.info.group.symbol)
    print(atom.info.block)

    fout = os.path.join(data_dir, 'PYDIRAC.inp')
    create_mrci_inp(fname, fout)



