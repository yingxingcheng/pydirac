# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Pydirac: PYthon tool for DIRAC software.
#  Copyright (C) 2020-2020 The Pydirac Development Team
#
#  This file is part of Pydirac.
#
#  Pydirac is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 3
#  of the License, or (at your option) any later version.
#
#  Pydirac is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/>
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import re
import warnings
from monty.json import MSONable, jsanitize
from pydirac.io.inputs import Inp, Mol
from pydirac.core.orbitals import OrbitalType, AtomicOrbital, MoleculeOrbitals


class Output:
    """Class to parse DIRAC output file

    Attributes:
        filename (str): output filename of DIRAC. Normally, a good output file
            contains `inp` file and `mol` file, thus one can parse them here
            together.

    """

    def __init__(self, filename):
        self.filename = filename

        self.inp = None
        self.mol = None

        # these are attributes of Output
        self.is_ok = False
        self.energy_settings = {}

        self.parse_energy()

    def parse_energy(self):
        """Parse the whole output file

        Parse the whole output file including `inp` file and `mol` file and
        `orbit` information and calculation `results` in the `output` file.

        Note:
            Input:
                1) SO or SR,
                2) zff,
                3) dipole or quadrupole,
                4) CC or CI,
                5) X2C or DOSSSS
            Mol:
                1) molecule
                2) basis sets: default (e.g., dyall.acv4z) or explicit

        Returns:
            None
        """

        self.check_valid()

        if self.is_ok:
            self.parse_input()
            self.parse_mol()
            # self.parse_orbit()

            self.parse_results()

    def parse_input(self):
        """Parse input from output file

        Returns:
            input settings
        """

        start_line = re.compile(r"^Contents of the input file\s*")
        dash_line = re.compile(r"^-+")
        empty_line = re.compile(r"^$")
        comment_line = re.compile(r"^[#|!]+")
        inp_start_line = re.compile(r"^\*\*DIRAC\s*")
        endline = re.compile(r"^\s*\*END OF")

        with open(self.filename, "r") as f:
            context = f.readlines()

        for i, line in enumerate(context):
            if start_line.match(line):
                context = context[i + 1 :]

        inp_str = []
        is_start = False
        for i, line in enumerate(context):
            if is_start and not endline.match(line):
                inp_str.append(line)
            elif inp_start_line.match(line):
                is_start = True
                inp_str.append(line)
            elif dash_line.match(line) or empty_line.match(line) or comment_line.match(line):
                continue
            else:
                inp_str.append(line)
                break

        self.inp = Inp.from_string(inp_str)

    def parse_mol(self):
        """Parse mol from output file

        Returns:
            Mol settings
        """
        start_line = re.compile(r"^Contents of the molecule file\s*")
        dash_line = re.compile(r"^-+")
        empty_line = re.compile(r"^$")
        # inp_start_line = re.compile(r"^\*\*DIRAC\s*")
        nuclei_nb_line = re.compile(
            r"\s*(?P<nuclei>[-+]?(\d+(\.\d*)?|\d*\.\d+))\s+(?P<nb_atoms>\d+)"
        )
        basis_line = re.compile(r"^\b(?:LARGE|EXPLICIT)\b\s+BASIS\s+(\S+)\s*")
        endline = re.compile(r"^FINISH")
        atomic_calc_ptn = re.compile(r"This is an atomic calculation")

        basis_type = None

        with open(self.filename, "r") as f:
            context = f.readlines()

        if atomic_calc_ptn.search("".join(context), re.MULTILINE | re.DOTALL):
            is_atom = True
        else:
            is_atom = False

        mol_str = []
        for i, line in enumerate(context):
            if start_line.match(line):
                # context = context[i + 1:]
                context = context[i + 3 :]
                break
        for i, line in enumerate(context):
            if endline.match(line):
                break
            else:
                mol_str.append(line)

        self.mol = Mol.from_string(mol_str)
        self.mol.molecule.is_atom = is_atom

    def parse_results(self):
        """Parse calculation results from output file

        Returns:
            None
        """

        if not self.inp:
            self.parse_input()

        if not self.mol:
            self.parse_mol()

        # check if calculation is finished
        if self.inp.calc_method in ["CI", "CC"]:
            if self.inp.calc_method == "CC":
                self._parse_coupled_cluster()
            else:
                self._parse_configuration_interaction()
        else:
            self.energy_settings = {}

        calc_type = self.inp.calc_type or "null"
        hamiltonian = self.inp.hamiltonian or "null"
        calc_method = self.inp.calc_method or "null"
        basis_type = self.mol.basis_type or "null"

        # for a set of calculations, this task_type can be regarded as an id
        self._task_type = "-".join([calc_type, calc_method, hamiltonian]) + "@" + basis_type.strip()

    def parse_orbit(self):
        """

        Args:
            filename:

            For SO
            * Fermion symmetry E1
                * Closed shell, f = 1.0000
                  -2.47797334037  ( 2)
                * Open shell #1, f = 0.5000
                  -0.13783115994  ( 2)
                * Virtual eigenvalues, f = 0.0000
                   0.00374548932  ( 6)        0.00499273090  ( 2)        0.01402796501  ( 6)        0.02179466435  (10)        0.03009547324  ( 2)
                   0.03716487091  ( 6)        0.06363549163  (14)        0.07606125973  (10)        0.09203286658  ( 6)        0.12208961664  ( 2)
                   0.19147871966  (14)        0.21707019780  (10)        0.22905557723  ( 2)        0.22905883825  ( 4)        0.26586221717  (18)
                   0.39444173531  ( 2)        0.49420939940  (14)        0.55623366579  ( 2)        0.55624525507  ( 4)        0.60447384566  ( 4)
                   0.60447617309  ( 6)        0.77688901369  (18)        1.25239748512  ( 2)        1.25997676379  ( 6)        1.25998182039  ( 8)
                   1.36133214856  ( 2)        1.36137529052  ( 4)        1.75072060544  ( 4)        1.75074035505  ( 6)        3.49702377155  ( 2)
                   3.49722182837  ( 4)        3.76715748401  ( 2)        9.76930432334  ( 2)        9.77034247125  ( 4)       10.47049017426  ( 2)
                  28.40952020368  ( 2)       36.12853032246  ( 2)       36.14097956915  ( 4)       79.15864898987  ( 2)      236.13953875683  ( 2)
                 783.60786051767  ( 2)     3002.96623263569  ( 2)    13731.73008778683  ( 2)

            For SR
            * Boson symmetry A1
              * Closed shell, f = 1.0000
                -2.47787107822  ( 2)
              * Open shell #1, f = 0.5000
                -0.13782962396  ( 2)
              * Virtual eigenvalues, f = 0.0000
                 0.01000493192  ( 2)        0.01881578996  ( 2)        0.03354206140  ( 2)        0.06105795349  ( 4)        0.08825402911  ( 2)
                 0.11027718518  ( 2)        0.15768768851  ( 4)        0.20553440208  ( 4)        0.22554311803  ( 2)        0.38377353600  ( 2)
                 0.46966416528  ( 4)        0.55315302017  ( 2)        0.59582793164  ( 4)        0.67809753506  ( 6)        1.24252600686  ( 2)
                 1.24308361482  ( 4)        1.35864972273  ( 2)        1.74463236465  ( 4)        3.49480571226  ( 2)        3.75812114959  ( 2)
                 9.76801822932  ( 2)       10.46222624083  ( 2)       28.40185035569  ( 2)       36.13522290609  ( 2)       79.15135049148  ( 2)
               236.13236364858  ( 2)      783.60069042839  ( 2)     3002.95967127327  ( 2)    13731.72642240539  ( 2)

            * Boson symmetry B1
              * Virtual eigenvalues, f = 0.0000
                 0.01000493192  ( 2)        0.03354206140  ( 2)        0.06105795349  ( 2)        0.08825402911  ( 2)        0.15768768851  ( 4)
                 0.20553440208  ( 2)        0.22554311803  ( 2)        0.46966416528  ( 4)        0.55315302017  ( 2)        0.59582793164  ( 2)
                 0.67809753506  ( 4)        1.24308361482  ( 4)        1.35864972273  ( 2)        1.74463236465  ( 2)        3.49480571226  ( 2)
                 9.76801822932  ( 2)       36.13522290609  ( 2)

        Returns:

        """
        start_ptn = re.compile(r"\s+Eigenvalues\s+")
        sym_ptn = re.compile(r"^\* \b(?:Boson|Fermion)\b symmetry (.*)")
        open_ptn = re.compile(r"^\s+\*\s+Open shell #\d+, f = (\d+(\.\d*)?)")
        closed_ptn = re.compile(r"^\s+\*\s+Closed shell, f = (\d+(\.\d*)?)")
        virtual_ptn = re.compile(r"\s+\*\s+Virtual eigenvalues, f = (\d+(\.\d*)?)")
        number = re.compile(r"(?P<num>[-+]?(\d+(\.\d*)?|\d*\.\d+))\s+\(\s*(?P<degen>\d+)\)")
        number_line_ptn = re.compile(r"\s+[-+]?(\d+(\.\d+)?|\d*\.\d+) .*")
        endline_ptn = re.compile(r"^\* HOMO - LUMO")

        with open(self.filename, "r") as f:
            lines = f.readlines()

        for i, l in enumerate(lines):
            if start_ptn.match(l):
                lines = lines[i + 1 :]
                break
        else:
            warnings.warn(
                "No SCF calculation info founded, this is calculation"
                " without SCF, please check it carefully!"
            )

        cur_sym = None
        ao_type = None
        occ_frac = 0.00
        aos = []
        for l in lines:
            symmetry_match = sym_ptn.match(l)
            if symmetry_match:
                # add new branch
                cur_sym = symmetry_match.group(1)

            openshell_match = open_ptn.match(l)
            if openshell_match:
                occ_frac = float(openshell_match.group(1))
                ao_type = OrbitalType.OPEN_SHELL

            closeshell_match = closed_ptn.match(l)
            if closeshell_match:
                occ_frac = float(closeshell_match.group(1))
                # occ_frac = 1.0
                ao_type = OrbitalType.CLOSED_SHELL

            virtualshell_match = virtual_ptn.match(l)
            if virtualshell_match:
                occ_frac = float(virtualshell_match.group(1))
                # occ_frac = 0.0
                ao_type = OrbitalType.VIRTUAL_SHELL

            if number_line_ptn.match(l):
                ao_energy, ao_degen = [], []
                for m in number.finditer(l):
                    ao_energy.append(float(m.group("num")))
                    ao_degen.append(int(m.group("degen")))
                for e, d in zip(ao_energy, ao_degen):
                    ao = AtomicOrbital(cur_sym, ao_type, e, int(d), occ_frac)
                    aos.append(ao)

            if endline_ptn.match(l):
                break

        self.mos = MoleculeOrbitals(aos)
        return self.mos

    @property
    def task_type(self):
        """Get a task type from an output file

        Returns:
        """
        return self._task_type

    def check_valid(self):
        """
            ********** E N D   of   D I R A C  output  **********
        Returns:

        """
        end_pattern = re.compile(r"\s*\*+\s+E N D\s+of\s+D I R A C\s+output\s+\*+")
        with open(self.filename, "r") as f:
            context = f.read()
        if end_pattern.search(context):
            self.is_ok = True
        else:
            self.is_ok = False
        return self.is_ok

    def _parse_coupled_cluster(self, has_T=True):
        """Deal with RELCCSD calculations

        Notes:
            An example with RELCCSD calculations

             Overview of calculated energies
            @ SCF energy :                             -5880.437843209574567
            @ MP2 correlation energy :                    -0.486037237589513
            @ CCSD correlation energy :                   -0.451464950662248
            @ 4th order triples correction :              -0.016640269726526
            @ 5th order triples (T) correction :           0.000318657059527
            @ 5th order triples -T  correction :           0.000362317692980
            @ Total MP2 energy :                       -5880.923880447164265
            @ Total CCSD energy :                      -5880.889308160236396
            @ Total CCSD+T  energy :                   -5880.905948429963246
            @ Total CCSD(T) energy :                   -5880.905629772903922
            @ Total CCSD-T  energy :                   -5880.905586112270612

        Returns:
            None
        """
        start_pattern = re.compile(r"^\s*Overview of calculated energies")
        scf_pattern = re.compile(r"^@ SCF energy :\s+(?P<energy>[-+]?(\d+(\.\d*)?|\d*\.\d+))")
        mp2_pattern = re.compile(r"^@ Total MP2 energy :\s+(?P<energy>[-+]?(\d+(\.\d*)?|\d*\.\d+))")
        ccsd_pattern = re.compile(
            r"^@ Total CCSD energy :\s+(?P<energy>[-+]?(\d+(\.\d*)?|\d*\.\d+))"
        )
        ccsd_p_T_pattern = re.compile(
            r"^@ Total CCSD\(T\) energy :\s+(?P<energy>[-+]?(\d+(\.\d*)?|\d*\.\d+))"
        )
        end_line = re.compile(r"^$")

        with open(self.filename, "r") as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if start_pattern.match(line):
                    lines = lines[i:]
                    break

        energies = {}
        if has_T:
            patterns = [scf_pattern, mp2_pattern, ccsd_pattern, ccsd_p_T_pattern]
            keys = ["scf_e", "mp2_e", "ccsd_e", "ccsd_p_T_e"]
        else:
            patterns = [scf_pattern, mp2_pattern, ccsd_pattern]
            keys = ["scf_e", "mp2_e", "ccsd_e"]

        for i, line in enumerate(lines):
            if end_line.match(line):
                break

            for pattern, k in zip(patterns, keys):
                m = pattern.match(line)
                if m:
                    energies[k] = float(m.group("energy"))
                    break

        self.energy_settings = energies
        self._parse_cc_orbit()

    def _parse_cc_orbit(self):
        # search orbital info
        #
        #  Configuration in highest pointgroup
        #                                           E    E
        #  Spinor class : occupied                  13   11
        #  Spinor class : virtual                  220  222
        #
        #  Configuration in abelian subgroup
        orb_start_pattern = re.compile(r"^\s*Configuration in highest pointgroup")
        occ_pattern = re.compile(r"^\s*Spinor class : occupied\s+((?:\d+\s+)+)")
        vir_pattern = re.compile(r"^\s*Spinor class : virtual\s+((?:\d+\s+)+)")
        orb_end_pattern = re.compile(r"^\s*Configuration in abelian subgroup")

        with open(self.filename, "r") as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if orb_start_pattern.match(line):
                    lines = lines[i:]
                    break

        clc_orb = {}
        keys = ["occ", "vir"]
        patterns = [occ_pattern, vir_pattern]

        for i, line in enumerate(lines):
            if orb_end_pattern.match(line):
                break

            for pattern, k in zip(patterns, keys):
                m = pattern.match(line)
                if m:
                    clc_orb[k] = sum([int(nb) for nb in m.group(1).split()])
                    break
        self._cc_calc_orbit = clc_orb

    def _parse_configuration_interaction(self):
        """Deal with KRCI calculations

        Notes:
            CI example:
                &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                &&& KRCI calculation for symmetry no.    3
                &&& Number of CI roots for this symmetry       3
                &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

                Final CI energies  =  -5880.8481422675641  -5880.8380283826846  -5880.8380141759062

                root    1 ...... converged!

                root    2 ...... converged!

                root    3 ...... converged!

        Returns:
            None
        """
        # start_line = re.compile(r'Final CI energies\s+=\s+(?P<energy>([-+]?(\d+(\.\d*)?|\d*\.\d+)\s+)+)\s+')
        energy_line = re.compile(
            r"Final CI energies\s+=\s+(?P<energy>(?:[-+]?(?:\d+(?:\.\d*)?|\d*\.\d+)\s+)+)\s+"
        )
        conv_line = re.compile(r"(?:\s*root\s+\d+\s+\.+\s+(?:converged|unconverged)!\s*)+")
        sym_line = re.compile(
            r"\s*&+\s+KRCI calculation for symmetry no\.\s+(\d+)\s+&+\s+Number of CI roots for this symmetry\s+(\d+)\s+"
        )

        with open(self.filename, "r") as fin:
            context = fin.read()

        energy_l = energy_line.findall(context, re.MULTILINE | re.DOTALL)
        conv_l = conv_line.findall(context, re.MULTILINE | re.DOTALL)
        sym_l = sym_line.findall(context, re.MULTILINE | re.DOTALL)

        energies = {}

        all_e = []
        for e_str in energy_l:
            e_roots = [float(e) for e in e_str.split()]
            all_e.append(e_roots)

        if len(sym_l) != len(all_e):
            warnings.warn("The shape of energy and symmetry are not " "the same! Check results!")
            return {}

        for i, sym_rt in enumerate(sym_l):
            tag_sym = sym_rt[0]
            nb_root = int(sym_rt[1])
            if nb_root != len(all_e[i]):
                warnings.warn(
                    "The number of root energies and the number of "
                    "roots are not the same! Check results!"
                )
                return {}

            for r in range(1, nb_root + 1):
                key = "_".join(["sym", str(tag_sym), "root", str(r)])
                energies[key] = float(all_e[i][r - 1])
                # energies['_'.join(['sym',str(tag_sym),'ave.'])] = sum(all_e[i])/len(all_e[i])

        convergeds = []
        for line in conv_l:
            converged = []
            for wrd in line.split("\n"):
                if not len(wrd.strip()):
                    continue
                if "unconverged" in wrd:
                    converged.append(False)
                elif "converged" in wrd:
                    converged.append(True)
                else:
                    print(wrd)
                    raise RuntimeError("There is no converged info")
            convergeds.append(converged)

        self.energy_settings = {"ci_e": energies, "ci_converged": convergeds}

    def _parse_hartree_fock(self):
        """Deal with HF calculations

        Returns:
            None
        """

        pass

    def _parse_visual_module(self):
        """Deal with visual module calculations

        Returns:
            None
        """

        pass

    @property
    def calc_orbit(self):
        if self.inp.calc_method == "CC":
            return self._cc_calc_orbit
        elif self.inp.calc_method == "CI":
            return self.inp.ci_calc_orbit
        else:
            return {"occ": 0, "vir": 0}

    def as_dict(self) -> dict:
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "energy_settings": self.energy_settings,
            "filename": self.filename,
            "task_type": self.task_type,
            "mol": self.mol.as_dict(),
            "inp": self.inp,
            "is_ok": self.is_ok,
        }
        return jsanitize(d, strict=True)
