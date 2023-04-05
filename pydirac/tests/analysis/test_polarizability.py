import importlib_resources
import pytest

from pydirac.analysis.polarizability import get_polarizability
from pydirac import Settings


data_root = importlib_resources.files("pydirac.tests.data")
He_q_so_dir = str(data_root / "He_q_so")
He_q_mrci_dir = str(data_root / "He_q_mrci")
He_so_dir = str(data_root / "He_so")
He_nr_dir = str(data_root / "He_nr")
He_mrci_dir = str(data_root / "He_mrci")
He_q_theta_dir = str(data_root / "He_q_theta")


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
