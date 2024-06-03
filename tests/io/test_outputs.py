import importlib.resources as pkg_resources
import pytest

from pydirac.io.outputs import Output

import tests.data  # Ensure your tests/data is a Python package (i.e., it contains an __init__.py file)

data_root = pkg_resources.files(tests.data)
K_mrci_dir = data_root / "K_mrci"

# out_fn1 = str(data_root / "Li" / "Li_D-CC-SR.out")
# out_fn2 = str(data_root / "Li" / "Li_D-CC-SO.out")
# ref_dict = {out_fn1: [162, 2, 2, 158], out_fn2: [212, 2, 2, 208]}
#
# @pytest.mark.parametrize("out_fn", [out_fn1, out_fn2])
# @pytest.mark.parametrize("e_min", [-10.0])
# @pytest.mark.parametrize("e_max", [10.0])
# def test_get_orbital_info(out_fn, e_min, e_max):
#     out = Output(filename=out_fn)
#     out.parse_orbit()
#     nb_tot = out.mos.nao(e_min, e_max)
#     nb_occ = out.mos.nb_closed_ao(e_min, e_max)
#     nb_open = out.mos.nb_open_ao(e_min, e_max)
#     nb_vir = out.mos.nb_virtual_ao(e_min, e_max)
#     assert nb_tot == nb_occ + nb_open + nb_vir
#     assert ref_dict[out_fn] == [nb_tot, nb_occ, nb_open, nb_vir]


def test_relcc():
    out_fn = str(
        data_root
        / "He_q_so"
        / "d-aug-dyall.acv3z_+0.00001"
        / "He_d-aug-dyall.acv3z_He_d-aug-dyall.acv3z_zff=+0.00001.out"
    )
    o = Output(out_fn)
    res = o.as_dict()
    assert res["energy_settings"]["scf_e"] == pytest.approx(-2.861794767985585)
    assert res["energy_settings"]["mp2_e"] == pytest.approx(-2.895006044817857)
    assert res["energy_settings"]["ccsd_e"] == pytest.approx(-2.900862077849925)
    assert res["energy_settings"]["ccsd(t)_e"] == pytest.approx(-2.900862077849925)
    assert res["task_type"] == "Q-4C-DC-CC@d-aug-dyall.acv3z"


def test_mrci():
    out_fn = str(
        data_root
        / "He_mrci"
        / "dyall.acv4z_+0.001/He_dyall.acv4z_He_dyall.acv4z_zff=+0.001.out"
    )
    o = Output(out_fn)
    res = o.as_dict()
    assert res["energy_settings"]["scf_e"] == pytest.approx(-2.8618113380597565)
    assert res["energy_settings"]["ci_e"]["sym_1_root_1"] == pytest.approx(
        -2.8975548136776
    )
    assert res["task_type"] == "D-4C-DC-CI@dyall.acv4z"
