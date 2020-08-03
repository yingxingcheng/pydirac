#!/usr/bin/env python
from pydirac.io.jobs import Inpobj
import os


module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))


def test_get_orbital_info():
    """
    test get orbital info
    """
    fname = os.path.join(data_dir, 'Kr_DHF.out')
    inp = Inpobj.from_file(filename=fname)
    fname_out = os.path.join(data_dir, 'tmp.inp')
    inp.write_to_file(filename=fname_out)



