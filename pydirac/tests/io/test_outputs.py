import importlib_resources
import pytest

from pydirac.io.outputs import Output

data_root = importlib_resources.files("pydirac.tests.data")
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
        / "In_q_so"
        / "d-aug-dyall.acv3z_+0.00001"
        / "In_d-aug-dyall.acv3z_In_d-aug-dyall.acv3z_zff=+0.00001.out"
    )
    o = Output(out_fn)
    res = o.as_dict()
    assert pytest.approx(res["energy_settings"]["scf_e"], -5880.437843209575)
    assert pytest.approx(res["energy_settings"]["mp2_e"], -5880.923880447164)
    assert pytest.approx(res["energy_settings"]["ccsd_e"], -5880.889308160236)
    assert pytest.approx(res["energy_settings"]["ccsd(t)_e"], -5880.905629772904)
    assert res["task_type"] == "Q-4C-DC-CC@d-aug-dyall.acv3z"


def test_mrci():
    out_fn = str(
        data_root
        / "K_mrci"
        / "d-aug-dyall.cv3z_+0.001/K_d-aug-dyall.cv3z_K_d-aug-dyall.cv3z_zff=+0.001.out"
    )
    o = Output(out_fn)
    res = o.as_dict()
    assert pytest.approx(res["energy_settings"]["scf_e"], -601.5260549871577)
    assert pytest.approx(res["energy_settings"]["ci_e"]["sym_3_root_1"], -601.7982384883599)
    assert res["task_type"] == "D-4C-DC-CI@d-aug-dyall.cv3z"
