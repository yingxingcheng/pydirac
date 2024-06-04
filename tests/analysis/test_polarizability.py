import importlib.resources as pkg_resources

import pytest

import tests.data  # Ensure your tests/data is a Python package (i.e., it contains an __init__.py file)
from pydirac import Settings
from pydirac.analysis.polarizability import get_polarizability

# Use pkg_resources to get the paths to the data directories
He_q_so_dir = str(pkg_resources.files(tests.data) / "He_q_so")
He_q_mrci_dir = str(pkg_resources.files(tests.data) / "He_q_mrci")
He_so_dir = str(pkg_resources.files(tests.data) / "He_so")
He_nr_dir = str(pkg_resources.files(tests.data) / "He_nr")
He_mrci_dir = str(pkg_resources.files(tests.data) / "He_mrci")
He_q_theta_dir = str(pkg_resources.files(tests.data) / "He_q_theta")


@pytest.mark.parametrize(
    "dirout_and_deepth",
    [
        (He_nr_dir, 1),
        (He_so_dir, 0),
        (He_mrci_dir, 0),
        (He_q_theta_dir, 1),
        (He_q_so_dir, 0),
        (He_q_mrci_dir, 0),
    ],
)
def test_calc_dipole_polarizability(dirout_and_deepth):
    res = get_polarizability(dirout_and_deepth[0], deepth=dirout_and_deepth[1])
    out = Settings(res)
    print(out)
