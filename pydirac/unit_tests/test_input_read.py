from pydirac.utility.read import parser_dirac_input
from pydirac.input.helper import DiracJob
import os

module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))
dirac_inp = os.path.join(data_dir, 'tmp.inp')


def test_read_func():
    setting = parser_dirac_input(filename=dirac_inp)
    job = DiracJob(settings=setting)
    print(job)

