import pytest
import numpy as np
from pydirac import Element


def test_element_group():
    assert Element(1).group == 1

    group_1 = np.asarray([3, 11, 19, 37, 55, 87])
    for j in range(2):
        for i in group_1 + j:
            assert Element(i).group == 1 + j

    for i in list(range(57, 71)) + list(range(89, 103)):
        assert Element(i).group == None

    group_3 = np.array([21, 39, 71, 103], dtype=int)
    for j in range(9):
        for i in group_3 + j:
            assert Element(i).group == 3 + j

    group_13 = np.array([5, 13, 31, 49, 81, 113], dtype=int)
    for j in range(5):
        for i in group_13 + j:
            assert Element(i).group == 13 + j


def test_element_initialization():
    hydrogen = Element("H")
    assert hydrogen.atomic_number == 1
    assert hydrogen.symbol == "H"

    oxygen = Element(8)
    assert oxygen.atomic_number == 8
    assert oxygen.symbol == "O"

    with pytest.raises(ValueError):
        Element(999)  # Non-existent element

    with pytest.raises(ValueError):
        Element("Xx")  # Non-existent element


def test_element_properties():
    carbon = Element("C")
    assert carbon.atomic_number == 6
    assert carbon.symbol == "C"
    assert carbon.block == "p"
    assert carbon.is_lanthanoid is False
    assert carbon.is_actinoid is False
    assert carbon.row == 2
    assert carbon.group == 14
    assert carbon.group_symbol == "group-14"
    assert carbon.Z == 6
    assert carbon.period == 2

    lanthanum = Element(57)
    assert lanthanum.is_lanthanoid is True

    uranium = Element(92)
    assert uranium.is_actinoid is True


def test_element_as_dict():
    neon = Element("Ne")
    neon_dict = neon.as_dict()
    expected_dict = {
        "Z": 10,
        "block": "p",
        "row": 2,
        "group": 18,
        "group_symbol": "group-18",
        "period": 2,
    }
    assert neon_dict == expected_dict


def test_element_get_elec_config():
    helium = Element("He")
    assert helium.get_elec_config() == "1s^2"

    chromium = Element(24)
    assert chromium.get_elec_config() == "1s^2 2s^2 2p^6 3s^2 3p^6 4s^1 3d^5"


def test_element_get_dirac_shell_str():
    helium = Element("He")
    dirac_str = helium.get_dirac_shell_str()
    expected_str = "# He atom and Econf: 1s^2\n" ".CLOSED SHELL\n" "2"
    assert dirac_str == expected_str

    dirac_dict = helium.get_dirac_shell_str(return_dict=True)
    expected_dict = {
        "CLOSED SHELL": "2",
    }
    assert dirac_dict == expected_dict
