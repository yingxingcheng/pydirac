import glob
import os

import importlib.resources as pkg_resources

from pydirac.analysis.utility import *
from pydirac.io.outputs import Output

import tests.data  # Ensure your tests/data is a Python package (i.e., it contains an __init__.py file)

# Use pkg_resources to get the paths to the data directories
He_q_so_dir = str(pkg_resources.files(tests.data) / "He_q_so")
He_q_mrci_dir = str(pkg_resources.files(tests.data) / "He_q_mrci")
He_so_dir = str(pkg_resources.files(tests.data) / "He_so")
He_nr_dir = str(pkg_resources.files(tests.data) / "He_nr")
He_mrci_dir = str(pkg_resources.files(tests.data) / "He_mrci")
He_q_theta_dir = str(pkg_resources.files(tests.data) / "He_q_theta")


def test_get_keyword():
    assert get_keyword("C", "energy", "(core 2)[vir 3]") == "C@energy@(core 2)[vir 3]"
    assert get_keyword("H", "gradient", "(occ 3)[vir 2]") == "H@gradient@(occ 3)[vir 2]"
    assert get_keyword("O", "hessian", "(occ 4)[vir 4]") == "O@hessian@(occ 4)[vir 4]"


def test_get_orbital_info():
    assert get_orbital_info(2, 3) == "(core 2)[vir 3]"
    assert get_orbital_info(3, 2) == "(core 3)[vir 2]"
    assert get_orbital_info(4, 4) == "(core 4)[vir 4]"


def test_get_energy():
    current = os.getcwd()
    os.chdir(He_q_so_dir)
    try:
        for calc_dir in glob.glob("*"):
            if os.path.isdir(calc_dir):
                os.chdir(os.path.join(He_q_so_dir, calc_dir))
                if "JOB_DONE" in glob.glob("*"):
                    outfiles = glob.glob("*.out")
                    if len(outfiles) > 1:
                        raise RuntimeError(
                            "There are two output files in the current directory: {}".format(
                                calc_dir
                            )
                        )
                    else:
                        outfile = outfiles[0]
                        energy = get_energy(outfile, method="CCSD(T)")
                        print(f"{calc_dir} CCSD(T): {energy}")
                os.chdir("..")
    finally:
        os.chdir(current)


def test_output_object():
    current = os.getcwd()
    os.chdir(He_q_so_dir)
    try:
        for calc_dir in glob.glob("*"):
            if os.path.isdir(calc_dir):
                os.chdir(os.path.join(He_q_so_dir, calc_dir))
                if "JOB_DONE" in glob.glob("*"):
                    outfiles = glob.glob("*.out")
                    if len(outfiles) > 1:
                        raise RuntimeError(
                            "There are two output files in the current directory: {}".format(
                                calc_dir
                            )
                        )
                    else:
                        outfile = outfiles[0]
                        output_obj = Output(outfile)
                        settings = output_obj.parse_input()
                os.chdir("..")
    finally:
        os.chdir(current)


def test_output_object_CI():
    current = os.getcwd()
    os.chdir(He_mrci_dir)
    try:
        for calc_dir in glob.glob("dyall*"):
            if os.path.isdir(calc_dir):
                os.chdir(os.path.join(He_mrci_dir, calc_dir))
                if "JOB_DONE" in glob.glob("*"):
                    outfiles = glob.glob("*.out")
                    if len(outfiles) > 1:
                        raise RuntimeError(
                            "There are two output files in the current directory: {}".format(
                                calc_dir
                            )
                        )
                    else:
                        outfile = outfiles[0]
                        output_obj = Output(outfile)
                        res = output_obj.parse_results()
                        print(output_obj.energy_settings)

                        res_dict = output_obj.as_dict()
                        print(res_dict)

                os.chdir("..")
    finally:
        os.chdir(current)
