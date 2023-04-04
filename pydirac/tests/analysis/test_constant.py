from pydirac.analysis.constant import get_keyword, get_orbital_info


def test_get_keyword():
    assert get_keyword("C", "energy", "(core 2)[vir 3]") == "C@energy@(core 2)[vir 3]"
    assert get_keyword("H", "gradient", "(occ 3)[vir 2]") == "H@gradient@(occ 3)[vir 2]"
    assert get_keyword("O", "hessian", "(occ 4)[vir 4]") == "O@hessian@(occ 4)[vir 4]"


def test_get_orbital_info():
    assert get_orbital_info(2, 3) == "(core 2)[vir 3]"
    assert get_orbital_info(3, 2) == "(core 3)[vir 2]"
    assert get_orbital_info(4, 4) == "(core 4)[vir 4]"
