from pydirac.read import parser_dirac_input
from pydirac.dirac import DiracJob
from pathlib import Path
import os

DATA_PATH = Path('data')
dirac_inp = os.path.join(DATA_PATH, 'tmp.inp')


def test_read_func():
    setting = parser_dirac_input(filename=dirac_inp)
    job = DiracJob(settings=setting)

    # with open('tmp2.inp', 'w') as f:
    #     f.write(job.get_input())
    print(job)

