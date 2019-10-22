#!/usr/bin/env python
import unittest
from pydirac.dirac_input import Inpobj
import os


module_dir = os.path.dirname(os.path.abspath(__file__))


class Test_Inpobj(unittest.TestCase):
    """
    Test class for Inpobj
    """

    def test_write_file(self, filename):
        """
        test write_file
        """
        inp = Inpobj.from_file(filename)
        inp.write_to_file(filename)

