from pydirac.utility.read import parse_dirac_input
from pydirac.input.jobs import DiracJob
import os

module_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(module_dir, 'data'))
dirac_inp = os.path.join(data_dir, 'tmp.inp')


def test_read_func():
    setting = parse_dirac_input(file_obj=dirac_inp)
    job = DiracJob(settings=setting)
    print(job)

