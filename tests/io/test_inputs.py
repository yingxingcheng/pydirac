import importlib.resources as pkg_resources

from pydirac.core.settings import Settings
from pydirac.io.inputs import Inp

import tests.data  # Ensure your tests/data is a Python package (i.e., it contains an __init__.py file)


def test_read_func():
    dirac_inp = str(
        pkg_resources.files(tests.data)
        / "He_mrci"
        / "dyall.acv4z_+0.001"
        / "He_dyall.acv4z.inp"
    )
    input = Inp.from_file(dirac_inp)
    print(input["DIRAC"])
    assert input["DIRAC"]["TITLE"] == [" He, DOSSSS, SCF"]
    assert input["DIRAC"]["ANALYZE"] == True
    assert input["DIRAC"]["WAVE F"] == True

    assert input["ANALYZE"]["MULPOP"]["_en"] == True
    assert input["ANALYZE"]["MULPOP"]["VECPOP"] == ["1..oo"]

    assert input["HAMILTONIAN"]["DOSSSS"] == True
    # not strip, so the space before keywords are kept
    assert input["HAMILTONIAN"]["OPERATOR"] == [" ZDIPLEN", " COMFACTOR", " zff"]
    assert input["INTEGRALS"]["READINP"]["UNCONTRACT"] == True

    assert input["WAVE FUNCTIONS"]["KR CI"] == True
    assert input["WAVE FUNCTIONS"]["RESOLVE"] == True
    assert input["WAVE FUNCTIONS"]["SCF"]["_en"] == True

    assert input["WAVE FUNCTIONS"]["SCF"]["CLOSED SHELL"] == ["2"]
    assert input["WAVE FUNCTIONS"]["SCF"]["EVCCNV"] == ["1.0D-9  5.0D-8"]
    assert input["WAVE FUNCTIONS"]["SCF"]["MAXITR"] == ["60"]

    assert input["WAVE FUNCTIONS"]["KRCICALC"]["CI PROGRAM"] == ["LUCIAREL"]
    assert input["WAVE FUNCTIONS"]["KRCICALC"]["INACTIVE"] == ["0"]
    assert input["WAVE FUNCTIONS"]["KRCICALC"]["GAS SHELLS"] == [
        "2",
        "0 2 / 1",
        "2 2 / 30",
    ]

    assert input["WAVE FUNCTIONS"]["KRCICALC"]["MAX CI"] == ["120"]
    assert input["WAVE FUNCTIONS"]["KRCICALC"]["MXCIVE"] == ["60"]
    assert input["WAVE FUNCTIONS"]["KRCICALC"]["NOOCCN"] == True
    assert input["WAVE FUNCTIONS"]["KRCICALC"]["RSTRCI"] == ["rstr"]


def test_write():
    # set single calculation parameters (single point, TZ2P/PW91)
    sett = Settings()
    sett.dirac.title = "B, DOSSSS, KRCI"
    sett.dirac.analyze = True
    wf_tag = "WAVE FUNCTIONS"
    sett.set_nested(("DIRAC", wf_tag), True)
    sett.analyze.mulpop._en = True
    sett.analyze.mulpop.vecpop = "1..oo"
    sett.hamiltonian.dossss = False
    sett.hamiltonian.x2c = True
    sett.hamiltonian.nospin = True
    sett.hamiltonian.operator = [" ZDIPLEN", " COMFACTOR", " 0.01"]
    sett.integrals.readinp.uncontract = True
    sett.general.pcmout = True

    wave_func = Settings()
    wave_func["KR CI"] = True
    wave_func.resolve = True

    scf = Settings()
    scf._en = True
    scf["CLOSED SHELL"] = 4
    scf["OPEN SHELL"] = [1, "1/6"]
    scf.evccnv = "1.0D-9 5.0D-8"
    scf.maxitr = 90

    wave_func.scf = scf
    sett[wf_tag] = wave_func

    krci = Settings()
    krci["CI PROGRAM"] = "LUCIAREL"
    krci["INACTIVE"] = 1
    krci["GAS SHELLS"] = [3, "0 2 / 1", "1 3 / 3", "3 3 / 10"]
    krci["MAX CI"] = 60
    krci["NOOCCN"] = True
    krci["DIPMOM"] = True
    krci["RSTRCI"] = 0
    krci["CIROOTS_id_0"] = "3  3"
    krci["CIROOTS_id_1"] = "4  3"
    sett.set_nested(("WAVE FUNCTIONS", "KRCICALC"), krci)

    job = Inp(sett.as_dict())

    ref_str = """**DIRAC
.TITLE
B, DOSSSS, KRCI
.ANALYZE
.WAVE FUNCTIONS
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
**HAMILTONIAN
.X2C
.NOSPIN
.OPERATOR
 ZDIPLEN
 COMFACTOR
 0.01
**INTEGRALS
*READINP
.UNCONTRACT
**GENERAL
.PCMOUT
**WAVE FUNCTIONS
.KR CI
.RESOLVE
.SCF
*SCF
.CLOSED SHELL
4
.OPEN SHELL
1
1/6
.EVCCNV
1.0D-9 5.0D-8
.MAXITR
90
*KRCICALC
.CI PROGRAM
LUCIAREL
.INACTIVE
1
.GAS SHELLS
3
0 2 / 1
1 3 / 3
3 3 / 10
.MAX CI
60
.NOOCCN
.DIPMOM
.RSTRCI
0
.CIROOTS
3  3
.CIROOTS
4  3
*END OF INPUT
"""
    assert str(job) == ref_str
