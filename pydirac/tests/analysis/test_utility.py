import importlib_resources
import os
import glob

from pydirac.analysis.utility import *
from pydirac.io.outputs import Output


data_root = importlib_resources.files("pydirac.tests.data")
He_q_so_dir = str(data_root / "He_q_so")
He_q_mrci_dir = str(data_root / "He_q_mrci")
He_so_dir = str(data_root / "He_so")
He_nr_dir = str(data_root / "He_nr")
He_mrci_dir = str(data_root / "He_mrci")
He_q_theta_dir = str(data_root / "He_q_theta")


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
    for calc_dir in glob.glob("*"):
        if os.path.isdir(calc_dir):
            os.chdir(os.path.join(He_q_so_dir, calc_dir))
            if "JOB_DONE" in glob.glob("*"):
                outfiles = glob.glob("*.out")
                if len(outfiles) > 1:
                    raise RuntimeError(
                        "There are two output file in current "
                        "directory: {0}".format(calc_dir)
                    )
                else:
                    outfile = outfiles[0]
                    energy = get_energy(outfile, method="CCSD(T)")
                    print("{0} CCSD(T): {1}".format(calc_dir, energy))
            os.chdir("..")
    os.chdir(current)


def test_output_object():

    current = os.getcwd()
    os.chdir(He_q_so_dir)
    for calc_dir in glob.glob("*"):
        if os.path.isdir(calc_dir):
            os.chdir(os.path.join(He_q_so_dir, calc_dir))
            if "JOB_DONE" in glob.glob("*"):
                outfiles = glob.glob("*.out")
                if len(outfiles) > 1:
                    raise RuntimeError(
                        "There are two output file in current "
                        "directory: {0}".format(calc_dir)
                    )
                else:
                    outfile = outfiles[0]
                    output_obj = Output(outfile)
                    settings = output_obj.parse_input()
            os.chdir("..")
    os.chdir(current)


def test_output_object_CI():

    current = os.getcwd()
    os.chdir(He_mrci_dir)
    for calc_dir in glob.glob("dyall*"):
        if os.path.isdir(calc_dir):
            os.chdir(os.path.join(He_mrci_dir, calc_dir))
            if "JOB_DONE" in glob.glob("*"):
                outfiles = glob.glob("*.out")
                if len(outfiles) > 1:
                    raise RuntimeError(
                        "There are two output file in current "
                        "directory: {0}".format(calc_dir)
                    )
                else:
                    outfile = outfiles[0]
                    output_obj = Output(outfile)
                    res = output_obj.parse_results()
                    print(output_obj.energy_settings)

                    res_dict = output_obj.as_dict()
                    print(res_dict)

            os.chdir("..")
    os.chdir(current)
