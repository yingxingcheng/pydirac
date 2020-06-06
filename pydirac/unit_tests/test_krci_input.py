#!/usr/bin/env python
import os

from pydirac.utility.get_orbit_info import Atom
from pydirac.input.krci import get_mrci_inp

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
    get_mrci_inp(fname, fout)
