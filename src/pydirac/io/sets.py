import abc
import glob
import os
import shutil
import warnings
from copy import deepcopy
from pathlib import Path

import importlib_resources
from monty.json import MSONable
from monty.serialization import loadfn

from pydirac.core.molecule import Molecule
from pydirac.core.settings import Settings
from pydirac.io.inputs import DiracInput, Inp, Mol
from pydirac.io.outputs import Output

__all__ = ["DiracInput", "AtomicDHFSet", "AtomicCCSet", "AtomicCISet"]

dir_path = os.path.dirname(__file__)

MODULE_DIR = Path(__file__).resolve().parent

"""
This module defines the DiracInputSet abstract base class and a concrete
implementation for the parameters developed and tested by the author.
The basic concept behind an input set is to specify a scheme to generate
a consistent set of DIRAC inputs from a structure without further user
intervention. This ensures comparability across runs.
"""


class DiracInputSet(MSONable, metaclass=abc.ABCMeta):
    """
    Abstract base class representing a set of Dirac input parameters with a
    structure supplied as initialization parameters.

    Attributes:
        inp (property): Inp object.
        mol (property): Mol object.
    """

    @property
    @abc.abstractmethod
    def inp(self):
        """Inp object"""
        pass

    @property
    @abc.abstractmethod
    def mol(self):
        """Mol object
        Returns:
        """
        pass

    def get_dirac_input(self):
        """Return a :class:`~pydirac.inputs.DiracInput` instance from the input and
        molecule objects.

        Returns:
            DiracInput: DiracInput instance.
        """
        return DiracInput(inp=self.inp, mol=self.mol)

    def write_input(self, output_dir, make_dir_if_not_present=True):
        """
        Write a set of DIRAC input to a directory.

        Args:
            output_dir (str): Path to the output directory.
            make_dir_if_not_present (bool): If True, create the directory if not present.

        Returns:
            None
        """
        dinput = self.get_dirac_input()
        dinput.write_input(output_dir, make_dir_if_not_present=make_dir_if_not_present)


def _load_yaml_config(fname):
    """
    Load a YAML configuration file and merge it with its parent configurations.

    Parameters
    ----------
    fname : str
        The name of the YAML configuration file to load.

    Returns
    -------
    dict
        The merged YAML configuration dictionary.

    """
    config = loadfn(str(importlib_resources.files("pydirac.data.template") / f"{fname}.yaml"))
    if "PARENT" in config:
        parent_config = _load_yaml_config(config["PARENT"])
        for k, v in parent_config.items():
            if k not in config:
                config[k] = v
            elif isinstance(v, dict):
                v_new = config.get(k, {})
                v_new.update(v)
                config[k] = v_new
    return config


class DictSet(DiracInputSet):
    """
    Concrete implementation of :class:`~pydirac.input_sets.DiracInputSet` that generates
    DIRAC input sets from a dictionary of settings and a molecule.

    Parameters
    ----------
    molecule : pymatgen.Molecule
        The molecule object.
    config_dict : dict
        A dictionary of input and molecule settings.
    files_to_transfer : dict, optional
        A dictionary of files to transfer.
    user_inp_settings : dict, optional
        A dictionary of user-defined input settings.
    user_mol_settings : dict, optional
        A dictionary of user-defined molecule settings.
    use_structure_charge : bool, optional
        Whether to use the structure charge.

    Attributes
    ----------
    inp : Inp
        Inp object.
    molecule : pymatgen.Molecule
        The molecule object.
    mol : Mol
        Mol object.

    Methods
    -------
    write_input(output_dir, make_dir_if_not_present=True)
        Write a set of DIRAC input to a directory.

    """

    def __init__(
        self,
        molecule,
        config_dict,
        files_to_transfer=None,
        user_inp_settings=None,
        user_mol_settings=None,
        use_structure_charge=False,
    ):
        """
        Constructor for :class:`~pydirac.input_sets.DictSet`.

        Parameters
        ----------
        molecule : pymatgen.Molecule
            The molecule object.
        config_dict : dict
            A dictionary of input and molecule settings.
        files_to_transfer : dict, optional
            A dictionary of files to transfer.
        user_inp_settings : dict, optional
            A dictionary of user-defined input settings.
        user_mol_settings : dict, optional
            A dictionary of user-defined molecule settings.
        use_structure_charge : bool, optional
            Whether to use the structure charge.

        """
        self._molecule = molecule
        self._config_dict = deepcopy(config_dict)
        self.files_to_transfer = files_to_transfer or {}
        self.user_inp_settings = user_inp_settings or {}
        self.user_mol_settings = user_mol_settings or {}
        self.use_molecule_charge = use_structure_charge

    @property
    def inp(self) -> Inp:
        """
        Inp object.

        Returns
        -------
        Inp
            Inp instance.

        """
        settings = dict(self._config_dict["inp"])
        inp = Inp(settings)
        return inp

    @property
    def molecule(self) -> Molecule:
        return self._molecule

    @property
    def mol(self) -> Mol:
        """
        Mol object.

        Returns
        -------
        Mol
            Mol instance.

        """
        mol = Mol(self.molecule)
        mol_settings = dict(self._config_dict["mol"])
        mol.basis_type = mol_settings.get("basis_type", "dyall.v2z")

        if "basis_type" in self.user_mol_settings:
            mol.basis_type = self.user_mol_settings["basis_type"]
        if "molecule" in self.user_mol_settings:
            mol.basis_type = self.user_mol_settings["molecule"]
        if "basis_lib" in self.user_mol_settings:
            mol.basis_lib = self.user_mol_settings["basis_lib"]
        if "comment" in self.user_mol_settings:
            mol.comment = self.user_mol_settings["comment"]
        return mol

    def __str__(self):
        """
        Return the class name.

        Returns
        -------
        str
            The class name.

        """
        return self.__class__.__name__

    def write_input(self, output_dir, make_dir_if_not_present=True):
        """
        Write a set of DIRAC input to a directory.

        Parameters
        ----------
        output_dir : str
            The path of the directory to write the input files to.
        make_dir_if_not_present : bool, optional
            Whether to create the directory if it does not exist.

        Returns
        -------
        None

        """
        super().write_input(output_dir=output_dir, make_dir_if_not_present=make_dir_if_not_present)
        for k, v in self.files_to_transfer.items():
            shutil.copy2(v, str(Path(output_dir) / k))


class AtomicDHFSet(DictSet):
    """
    A class representing a set of input parameters for an atomic Dirac-Hartree-Fock (DHF) calculation.

    Parameters
    ----------
    molecule : Molecule
        The molecule to compute the DHF calculation for.
    hamiltonian_mode : str, optional
        The type of Hamiltonian to use. Must be one of "2C" (two-component) or "4C" (four-component).
        Defaults to "4C".
    is_spinfree : bool, optional
        Whether to use a spin-free Hamiltonian. Defaults to False.
    is_ff : bool, optional
        Whether to include a frozen-core approximation. Defaults to False.
    ff_mode : str, optional
        The type of frozen-core approximation to use. Must be one of "D" (dipole approximation) or "Q"
        (quadrupole approximation). Defaults to "D".
    zff : float, optional
        The factor by which to scale the frozen-core potential. Only used if is_ff is True. Defaults to None.
    **kwargs
        Additional keyword arguments to be passed to the superclass constructor.

    Attributes
    ----------
    SUPPORTED_FF_MODES : tuple of str
        The supported frozen-core approximation modes.
    SUPPORTED_HAMILTONIAN_MODES : tuple of str
        The supported Hamiltonian modes.

    Raises
    ------
    ValueError
        If the ff_mode parameter is not one of the supported modes.

    Notes
    -----
    This class inherits from DictSet and uses the configuration file "AtomicDHFSet" to populate its input dictionary.

    """

    CONFIG = _load_yaml_config("AtomicDHFSet")

    SUPPORTED_FF_MODES = ("D", "Q")
    SUPPORTED_HAMILTONIAN_MODES = ("2C", "4C")

    def __init__(
        self,
        molecule,
        hamiltonian_mode="4C",
        is_spinfree=False,
        is_ff=False,
        ff_mode="D",
        zff=None,
        **kwargs,
    ):
        super().__init__(molecule, AtomicDHFSet.CONFIG, **kwargs)
        self.hamiltonian_mode = hamiltonian_mode
        self.is_spinfree = is_spinfree
        self.is_ff = is_ff
        self.ff_mode = ff_mode
        self.zff = zff
        self.kwargs = kwargs

        if self.ff_mode not in AtomicDHFSet.SUPPORTED_FF_MODES:
            raise ValueError(
                "{} not one of the support modes : {}".format(
                    self.ff_mode, AtomicDHFSet.SUPPORTED_FF_MODES
                )
            )

    @property
    def inp(self):
        """
        The input dictionary for the DHF calculation.

        Returns
        -------
        dict
            The input dictionary.
        """
        parent_inp = super().inp

        # update open shell setup
        if self.molecule.is_atom:
            elec_dict = self.molecule.atomic_info.get_dirac_shell_str(return_dict=True)
            for k in ["CLOSED SHELL", "OPEN SHELL"]:
                if k in elec_dict:
                    parent_inp[parent_inp.wf_tag]["SCF"].update({k: elec_dict[k]})
        else:
            # default is closed shell
            pass

        # update hamiltonian, 4c and spinfree
        parent_inp.update(
            {
                "HAMILTONIAN": {
                    "DOSSSS": self.hamiltonian_mode == "4C",
                    "X2C": self.hamiltonian_mode == "2C",
                    "NOSPIN": self.is_spinfree,
                }
            }
        )

        # update operator
        if self.is_ff:
            if self.zff is None:
                zff = "zff"
            else:
                zff = float(self.zff)
            if "HAMILTONIAN" in parent_inp:
                if self.ff_mode.upper() == "Q":
                    operator_dict = {
                        "OPERATOR": [
                            " Theta quadru-field",
                            " DIAGONAL",
                            " ZZTHETA",
                            " COMFACTOR",
                            f" {zff}",
                        ]
                    }
                else:
                    operator_dict = {"OPERATOR": [" ZDIPLEN", " COMFACTOR", f" {zff}"]}
                parent_inp["HAMILTONIAN"].update(operator_dict)
        return parent_inp


class AtomicCCSet(AtomicDHFSet):
    """
    Specification of reference determinant, type of calculation, and general settings.

    Default is energy calculation.

    Parameters
    ----------
    molecule : AtomicStruct or Molecule
        The atomic or molecular structure.
    e_min : float, optional
        Minimum energy of the active space. Default is -20.0.
    e_max : float, optional
        Maximum energy of the active space. Default is 25.0.
    e_error : float, optional
        Convergence threshold for energy. Default is 0.01.
    maxit : int, optional
        Maximum number of iterations. Default is 60.
    ntol : int, optional
        Convergence threshold for the norm of the residual. Default is 10.
    nelec : tuple of int, optional
        Number of electrons for closed-shell systems. Default is None.
    nelec_open : tuple of int, optional
        Number of electrons for open-shell systems. Default is None.
    nel_f1 : tuple of int, optional
        Number of frozen core electrons in the first fragment. Default is None.
    nel_f2 : tuple of int, optional
        Number of frozen core electrons in the second fragment. Default is None.
    no_T : bool, optional
        If True, does not include the T amplitudes. Default is None.
    print_level : int, optional
        Level of detail printed. Default is 1.
    prev_inp : str or Inp, optional
        Path to previous input file or the previous input file object. Default is None.
    prev_mol : str or Mol, optional
        Path to previous molecular file or the previous molecular file object. Default is None.
    **kwargs
        Additional arguments to pass to AtomicDHFSet.

    Attributes
    ----------
    e_min : float
        Minimum energy of the active space.
    e_max : float
        Maximum energy of the active space.
    e_error : float
        Convergence threshold for energy.
    maxit : int
        Maximum number of iterations.
    ntol : int
        Convergence threshold for the norm of the residual.
    nelec : tuple of int or None
        Number of electrons for closed-shell systems.
    nelec_open : tuple of int or None
        Number of electrons for open-shell systems.
    nel_f1 : tuple of int or None
        Number of frozen core electrons in the first fragment.
    nel_f2 : tuple of int or None
        Number of frozen core electrons in the second fragment.
    no_T : bool or None
        If True, does not include the T amplitudes.
    print_level : int
        Level of detail printed.
    prev_inp : Inp or None
        Previous input file object or None.
    prev_mol : Mol or None
        Previous molecular file object or None.
    kwargs : dict
        Additional arguments passed to AtomicDHFSet.

    """

    def __init__(
        self,
        molecule,
        e_min=-20.0,
        e_max=25.0,
        e_error=0.01,
        maxit=60,
        ntol=10,
        nelec=None,
        nelec_open=None,
        nel_f1=None,
        nel_f2=None,
        no_T=None,
        print_level=1,
        prev_inp=None,
        prev_mol=None,
        **kwargs,
    ):
        super().__init__(molecule, **kwargs)
        # self._config_dict['inp'].update()
        if isinstance(prev_inp, str):
            prev_inp = Inp.from_file(prev_inp)
        if isinstance(prev_mol, str):
            prev_mol = Mol.from_file(prev_mol)

        self.e_min = e_min
        self.e_max = e_max
        self.e_error = e_error
        self.maxit = maxit
        self.ntol = ntol
        self.nelec = nelec
        self.nelec_open = nelec_open
        self.nel_f1 = nel_f1
        self.nel_f2 = nel_f2
        self.no_T = no_T
        self.print_level = print_level
        self.prev_inp = prev_inp
        self.prev_mol = prev_mol
        self.kwargs = kwargs

    @property
    def inp(self):
        """
        Generates an input file object with the settings required for a CC calculation.

        Returns
        -------
        Inp
            Input file object.
        """
        parent_inp = super().inp
        inp = self.prev_inp if self.prev_inp is not None else parent_inp

        s = Settings(inp)
        # no matter whether previous inp has 4index
        s.DIRAC["4INDEX"] = True
        s[inp.wf_tag].RELCCSD = True
        # s.integrals.readinp.uncontract = True
        s.INTEGRALS.READINP = "UNCONTRACT"
        if "active" not in s.MOLTRA:
            s.MOLTRA.ACTIVE = f"{self.e_min} {self.e_max} {self.e_error}"

        # RELCC setup
        RELCC = Settings()
        RELCC.ENERGY = True
        RELCC.PRINT = self.print_level
        RELCC.CCENER = True
        RELCC.MAXIT = self.maxit
        RELCC.NTOL = self.ntol
        if "NOSDT" not in s.RELCC and self.no_T is True:
            RELCC.NOSDT = self.no_T

        nelec = self.nelec_open or self.nelec
        if nelec and len(nelec):
            RELCC.NELEC = " ".join([str(i) for i in nelec])
        if self.nel_f1 and len(self.nel_f1):
            RELCC.NEL_F1 = " ".join([str(i) for i in self.nel_f1])
        if self.nel_f2 and len(self.nel_f2):
            RELCC.NEL_F2 = " ".join([str(i) for i in self.nel_f1])
        s.RELCC = RELCC
        # update RELCC setup
        inp.update(s.as_dict())
        inp.update(self.user_inp_settings)
        return inp

    @property
    def mol(self):
        """
        Returns the molecular file object.

        Returns
        -------
        Mol
            Molecular file object.
        """
        parent_mol = super().mol
        this_mol = self.prev_mol if self.prev_mol is not None else parent_mol
        this_mol.update(self.user_mol_settings)
        return this_mol

    def override_from_prev_calc(self, prev_calc_dir=".", no_scf=False):
        """
        Update the input set to include settings from a previous calculation.

        Parameters
        ----------
        prev_calc_dir : str, optional
            Path to the previous calculation directory. Default is ".".
        no_scf : bool, optional
            If True, does not perform SCF calculations. Default is False.

        Returns
        -------
        AtomicCCSet
            Input set with updated settings.
        """
        output = get_diracrun_output(prev_calc_dir)

        self.prev_inp = output.inp
        self.prev_mol = output.mol
        # TODO: we do not need to override a DHF calculation, only one maybe .SCF
        if no_scf:
            self.prev_inp[self.prev_inp.wf_tag]["SCF"].update({"_en": False})

        # TODO, the only need to change is the number of active electrons
        # specifiec by .NELEC keyword or NEL_F1
        return self

    @classmethod
    def from_prev_calc(cls, prev_calc_dir, ff_mode="D", no_scf=False, **kwargs):
        """
        Create an AtomicCCSet object using settings from a previous calculation.

        Parameters
        ----------
        prev_calc_dir : str
            Path to the previous calculation directory.
        ff_mode : str, optional
            Fragmentation mode. Default is "D".
        no_scf : bool, optional
            If True, does not perform SCF calculations. Default is False.
        **kwargs
            Additional arguments to pass to AtomicDHFSet.

        Returns
        -------
        AtomicCCSet
            Input set with updated settings from previous calculation.
        """
        input_set = cls(_dummy_structure, ff_mode=ff_mode, **kwargs)
        return input_set.override_from_prev_calc(prev_calc_dir, no_scf)


class AtomicCISet(AtomicDHFSet):
    """
    Specifies the reference determinant and type of calculation for CI calculations.

    Inherits from AtomicDHFSet.

    Parameters
    ----------
    molecule : str or pymatgen.Molecule
        A string or pymatgen.Molecule object representing the molecule to calculate.
    prev_inp : str or Inp, optional
        A string or Inp object representing the previous input file to be used as a template.
        Default is None.
    prev_mol : str or Mol, optional
        A string or Mol object representing the previous molecular file to be used as a template.
        Default is None.
    **kwargs
        Additional keyword arguments to pass to AtomicDHFSet.

    Attributes
    ----------
    prev_inp : Inp
        The previous input file to be used as a template.
    prev_mol : Mol
        The previous molecular file to be used as a template.
    kwargs : dict
        Additional keyword arguments passed to AtomicDHFSet.

    Notes
    -----
    The default calculation is energy calculation.
    """

    def __init__(self, molecule, prev_inp=None, prev_mol=None, **kwargs):
        super().__init__(molecule, **kwargs)
        if isinstance(prev_inp, str):
            prev_inp = Inp.from_file(prev_inp)
        if isinstance(prev_mol, str):
            prev_mol = Mol.from_string(prev_mol)

        self.prev_inp = prev_inp
        self.prev_mol = prev_mol
        self.kwargs = kwargs

    @property
    def inp(self):
        """
        Generates an input file object with the settings required for a CI calculation.

        Returns
        -------
        Inp
            Input file object.
        """
        inp = super().inp

        # Step 1: the default KRCI setup
        s = Settings(inp)
        # s.dirac.analyze = True
        # wf_tag = inp.wf_tag or 'WAVE FUNCTIONS'
        # s.dirac[wf_tag] = True
        # s.analyze.mulpop._en = True
        # s.analyze.mulpop.vecpop = '1..oo'
        # s.integrals.readinp = 'uncontract'

        # wave_func = Settings()
        # wave_func['KR CI'] = True
        # wave_func.resolve = True
        # s[wf_tag] = wave_func

        # krci = Settings()
        # krci['CI PROGRAM'] = 'LUCIAREL'
        # # TODO: for different elements, the krci setup is different
        # # krci.inactive = 1
        # # krci['GAS SHELLS'] = [3, '0 2 / 1', '1 3 / 3', '3 3 / 10']
        # krci['MAX CI'] = 60
        # krci['NOOCCN'] = True
        # krci.nooccn = True
        # krci.rstrci = 0
        # # TODO: for different element, the number of roots is also different
        # krci['CIROOTS_id_0'] = '3  3'
        # krci['CIROOTS_id_1'] = '4  3'
        # s[wf_tag].krcicalc = krci
        inp.update(s.as_dict())

        # Step 2: if the previous inp exist, we use the previous one
        if self.prev_inp:
            # pre_s = Settings(self.prev_inp)
            # print(pre_s)
            # inp.update(pre_s.as_dict())
            inp = self.prev_inp
        else:
            warnings.warn("For KRCI module, more detail need to be considered")
            return None

        # Step 3: user setting
        inp.update(self.user_inp_settings)
        return inp

    @property
    def mol(self):
        """
        Returns the molecular file object.

        Returns
        -------
        Mol
            Molecular file object.
        """
        parent_mol = super().mol
        this_mol = self.prev_mol if self.prev_mol is not None else parent_mol
        this_mol.update(self.user_mol_settings)
        return this_mol

    def override_from_prev_calc(self, prev_calc_dir=".", no_scf=False):
        """
        Update the input set to include settings from a previous calculation.

        Parameters
        ----------
        prev_calc_dir : str, optional
            Path to the previous calculation directory. Default is ".".
        no_scf : bool, optional
            If True, does not perform SCF calculations. Default is False.

        Returns
        -------
        AtomicCISet
            Input set with updated settings.
        """
        output = get_diracrun_output(prev_calc_dir)

        self.prev_inp = output.inp
        self.prev_mol = output.mol
        if no_scf:
            self.prev_inp[self.prev_inp.wf_tag]["SCF"].update({"_en": False})
        return self

    @classmethod
    def from_prev_calc(cls, prev_calc_dir, ff_mode="D", no_scf=False, **kwargs):
        """
        Create an AtomicCISet object using settings from a previous calculation.

        Parameters
        ----------
        prev_calc_dir : str
            Path to the previous calculation directory.
        ff_mode : str, optional
            Fragmentation mode. Default is "D".
        no_scf : bool, optional
            If True, does not perform SCF calculations. Default is False.
        **kwargs
            Additional arguments to pass to AtomicDHFSet.

        Returns
        -------
        AtomicCISet
            Input set with updated settings from previous calculation.
        """
        input_set = cls(_dummy_structure, ff_mode=ff_mode, **kwargs)
        return input_set.override_from_prev_calc(prev_calc_dir, no_scf)

    @staticmethod
    def get_info_for_s_block(out: Output, e_min, e_max, calc_type="quadrupole"):
        """
        Get information for the s block.

        Parameters
        ----------
        out : Output
            The output file to parse.
        e_min : float
            The minimum energy.
        e_max : float
            The maximum energy.
        calc_type : str, optional
            The type of calculation. Default is "quadrupole".

        Returns
        -------
        nb_active_elec: int
            The number of electrons in the GAS.
        gas_list: list
            A list for GAS setup.
        root_list: list
            Root for different symmetries.
        """
        atomic_info = out.mol.molecule.atomic_info
        mos = out.mos
        assert atomic_info.block == "s"

        (
            nb_closed_elec,
            nb_open_elec,
            nb_total_elec,
            nb_closed_shell,
            nb_open_shell,
            nb_vir_shell,
        ) = mos.get_ao_and_elec(e_min, e_max)

        assert nb_closed_elec == nb_closed_shell
        assert nb_open_elec <= nb_open_shell
        period = atomic_info.period

        gas_list = []
        root_list = []

        if nb_open_shell == 0:
            # for nobel gas elements
            if period == 1:
                # He
                nb_active_elec = 2
                gas1 = "0 2 / 1"
                gas2 = f"2 2 / {nb_vir_shell // 2}"
                gas_list.extend([gas1, gas2])
            elif period == 2:
                # Be
                nb_active_elec = 4
                gas1 = "2 4 / 2"
                gas2 = f"4 4 / {nb_vir_shell // 2}"
                gas_list.extend([gas1, gas2])
            elif period == 3:
                # Mg
                nb_active_elec = 10
                gas1 = "0 2 / 1"
                gas2 = "8 10 / 4"
                gas3 = f"10 10 / {nb_vir_shell // 2}"
                gas_list.extend([gas1, gas2, gas3])
                pass
            else:
                # Ca Sr Ba Ra
                # including d electrons
                nb_active_elec = 20
                gas1 = "10 12 / 6 ! (n-1)s2 (n-2)d10"
                gas2 = "18 20 / 4 ! (n-1)p6 (n)s2"
                gas3 = f"20 20 / {nb_vir_shell // 2} ! virtual orbitals"
                gas_list.extend([gas1, gas2, gas3])
        else:
            # for s-block open-shell elements
            if period > 3:
                # K Rb Cs Fr
                # including d electrons
                nb_active_elec = 19
                gas1 = "10 12 / 6 ! (n-1)s2 (n-2)d10"
                gas2 = "17 19 / 4 ! (n-1)p6 (n)s1"
                gas3 = f"19 19 / {nb_vir_shell // 2} ! virtual orbitals"
                gas_list.extend([gas1, gas2, gas3])
            elif period == 3:
                # for Na
                nb_active_elec = 9
                gas1 = "0 2 / 1"
                gas2 = "7 9 / 4"
                gas3 = f"9 9 / {nb_vir_shell // 2}"
                gas_list.extend([gas1, gas2, gas3])
            elif period == 2:
                # for Li
                nb_active_elec = 3
                gas1 = "1 3 / 2 ! 1s"
                gas2 = f"3 3 / {nb_vir_shell // 2}"
                gas_list.extend([gas1, gas2])
            else:
                # H
                assert period == 1
                nb_active_elec = 1
                gas1 = "0 1 / 1 ! 1s"
                gas2 = f"1 1 / {nb_vir_shell // 2}"
                gas_list.extend([gas1, gas2])

        # specify many roots should we compute
        if nb_open_elec == 1:
            # IA: H Li Na K Rb Cs Fr
            ciroot1 = "3 1"
            root_list.extend([ciroot1])
        elif nb_open_elec == 0:
            # IIA: He Be Mg Ca Sr Ba Ra
            ciroot1 = "1 1"
            # ciroot2 = '2 1'
            root_list.extend([ciroot1])
        else:
            raise RuntimeError(
                "The number of electrons in open shell " "is wrong: {}".format(nb_open_elec)
            )

        return nb_active_elec, gas_list, root_list

    @staticmethod
    def get_info_for_p_block(out: Output, e_min, e_max, calc_type="quadrupole"):
        """
        Get information for the p block.

        Parameters
        ----------
        out : Output
            The output file to parse.
        e_min : float
            The minimum energy.
        e_max : float
            The maximum energy.
        calc_type : str, optional
            The type of calculation. Default is "quadrupole".

        Returns
        -------
        nb_active_elec: int
            The number of electrons in the GAS.
        gas_list: list
            A list for GAS setup.
        root_list: list
            Root for different symmetries.
        """
        atomic_info = out.mol.molecule.atomic_info
        mos = out.mos
        assert atomic_info.block == "p"

        (
            nb_closed_elec,
            nb_open_elec,
            nb_total_elec,
            nb_closed_shell,
            nb_open_shell,
            nb_vir_shell,
        ) = mos.get_ao_and_elec(e_min, e_max)

        assert nb_closed_elec == nb_closed_shell
        assert nb_open_elec <= nb_open_shell
        period = atomic_info.period

        gas_list = []
        root_list = []

        if nb_open_shell == 0:
            # for nobel gas elements
            # if period == 1:
            #     # He
            #     nb_active_elec = 2
            #     gas1 = '0 2 / 1'
            #     gas2 = '2 2 / {0}'.format(nb_vir_shell // 2)
            #     gas_list.extend([gas1, gas2])
            if period == 2:
                # Ne
                nb_active_elec = 10
                gas1 = "2 4 / 2"
                gas2 = "8 10 / 3"
                gas3 = f"10 10 / {nb_vir_shell // 2}"
                gas_list.extend([gas1, gas2, gas3])
            else:
                # Ar, Kr, Xe, Rn, Og
                # including d electrons
                nb_active_elec = 18
                gas1 = "10 12 / 6 ! (n-1)s2 (n-2)d10"
                gas2 = "16 18 / 3 ! (n-1)p6"
                gas3 = f"18 18 / {nb_vir_shell // 2} ! virtual orbitals"
                gas_list.extend([gas1, gas2, gas3])
        else:
            # for p-block open-shell elements
            if period > 3:
                # including d electrons
                nb_active_elec = 10 + 2 + nb_open_elec
                gas1 = "10 12 / 6 ! (n-1)s2 (n-2)d10"
                gas2 = f"{nb_open_elec + 10} {nb_open_elec + 12} / 3 ! (n-1)p6"
                gas3 = "{} {} / {} ! virtual orbitals".format(
                    nb_open_elec + 12, nb_open_elec + 12, nb_vir_shell // 2
                )
                gas_list.extend([gas1, gas2, gas3])
            elif period == 3:
                # for Al Si P S Cl
                nb_active_elec = 2 + nb_open_elec
                gas1 = "0 2 / 1"
                gas2 = f"{nb_open_elec} {nb_open_elec + 2} / 3"
                gas3 = f"{nb_open_elec + 2} {nb_open_elec + 2} / {nb_vir_shell // 2}"
                gas_list.extend([gas1, gas2, gas3])
            else:
                # without d electrons
                assert period == 2
                # for B C N O F
                nb_active_elec = nb_total_elec
                gas1 = "2 4 / 2 ! 1s2s"
                gas2 = f"{nb_total_elec - 2} {nb_total_elec} / 3 ! 2p"
                gas3 = f"{nb_total_elec} {nb_total_elec} / {nb_vir_shell // 2}"
                gas_list.extend([gas1, gas2, gas3])

        # specify many roots should we compute
        if nb_open_elec == 1:
            # IIIA: B Al Ga In Tl Nh
            ciroot1 = "3 3"
            ciroot2 = "4 3"
            root_list.extend([ciroot1, ciroot2])
        elif nb_open_elec == 2:
            # IVA: C Si Ge Sn Pb Fl
            if calc_type == "dipole":
                # DIRAC bugs, only one root does not get exact gs state
                ciroot1 = "1 2"
                # ciroot2 = '2 1'
                root_list.extend([ciroot1])
            elif calc_type == "quadrupole":
                # ciroot1 = '1 9'
                # ciroot2 = '2 6'
                ciroot1 = "1 2"
                ciroot2 = "2 2"
                root_list.extend([ciroot1, ciroot2])
            else:
                raise TypeError('Calc_type wrong! "dipole" or "quadrupole"')
        elif nb_open_elec == 3:
            # VA: N P As Sb Bi Mc
            if calc_type == "dipole":
                ciroot1 = "3 2"
                ciroot2 = "4 2"
                root_list.extend([ciroot1, ciroot2])
            elif calc_type == "quadrupole":
                # ciroot1 = '3 10'
                # # ciroot2 = '4 10'
                ciroot1 = "3 2"
                # ciroot2 = '4 2'
                root_list.extend([ciroot1])
            else:
                raise TypeError('Calc_type wrong! "dipole" or "quadrupole"')
        elif nb_open_elec == 4:
            # # VIA: O S Se Te Po Lv
            if calc_type == "dipole":
                ciroot1 = "1 3"
                ciroot2 = "2 2"
                root_list.extend([ciroot1, ciroot2])
            elif calc_type == "quadrupole":
                # ciroot1 = '1 9'
                # ciroot2 = '2 6'
                ciroot1 = "1 3"
                ciroot2 = "2 2"
                root_list.extend([ciroot1, ciroot2])
            else:
                raise TypeError('Calc_type wrong! "dipole" or "quadrupole"')
        elif nb_open_elec == 5:
            # VIIA: F Cl Br I At Ts
            ciroot1 = "3 3"
            ciroot2 = "4 3"
            root_list.extend([ciroot1, ciroot2])
        elif nb_open_elec == 0:
            # closed-shell elementss
            ciroot1 = "1 1"
            # ciroot2 = '2 1'
            root_list.extend([ciroot1])
        else:
            raise RuntimeError(
                "The number of electrons in open shell " "is wrong: {}".format(nb_open_elec)
            )

        return nb_active_elec, gas_list, root_list

    @staticmethod
    def get_info_for_d_block(out: Output, e_min, e_max, calc_type="quadrupole"):
        """
        Get information for the d block.

        Parameters
        ----------
        out : Output
            The output file to parse.
        e_min : float
            The minimum energy.
        e_max : float
            The maximum energy.
        calc_type : str, optional
            The type of calculation. Default is "quadrupole".

        Returns
        -------
        nb_active_elec: int
            The number of electrons in the GAS.
        gas_list: list
            A list for GAS setup.
        root_list: list
            Root for different symmetries.
        """
        atomic_info = out.mol.molecule.atomic_info
        mos = out.mos
        assert atomic_info.block == "d"

        (
            nb_closed_elec,
            nb_open_elec,
            nb_total_elec,
            nb_closed_shell,
            nb_open_shell,
            nb_vir_shell,
        ) = mos.get_ao_and_elec(e_min, e_max)

        assert nb_closed_elec == nb_closed_shell
        assert nb_open_elec <= nb_open_shell
        period = atomic_info.period

        gas_list = []
        root_list = []

        if nb_open_shell == 0:
            # IIB
            # Zn Cd Hg Cn
            # including d electrons
            nb_active_elec = 12
            gas1 = "8 10 / 5 ! (n-1)d10"
            gas2 = "10 12 / 1 ! (n)s2"
            gas3 = f"12 12 / {nb_vir_shell // 2} ! virtual orbitals"
            gas_list.extend([gas1, gas2, gas3])
        else:
            # for d-block open-shell elements
            if nb_open_shell == 2:
                # s is open
                # Cu Ag Au Rg
                assert nb_open_elec == 1
                nb_active_elec = 11
                gas1 = "8 10 / 5 ! (n-1)d10"
                gas2 = "9 11 / 1 ! (n)s1"
                gas3 = f"11 11 / {nb_vir_shell // 2} ! virtual orbitals"
                gas_list.extend([gas1, gas2, gas3])
            elif nb_open_shell == 10:
                # d is open
                nb_active_elec = 2 + nb_open_elec
                # normally there is no this case, since in DIRAC (n-1)d locates
                # before (n)s shells
                gas1 = "0 2 / 1 ! (n)s2"
                gas2 = f"{nb_open_elec} {nb_open_elec + 2} / 5 ! (n-1)d10"
                gas3 = "{0} {0} / {1} ! virtual orbitals".format(
                    nb_open_elec + 2, nb_vir_shell // 2
                )
                gas_list.extend([gas1, gas2, gas3])
            elif nb_open_shell == 12:
                # s + d are open
                nb_active_elec = nb_open_elec
                gas1 = f"{nb_open_elec - 2} {nb_open_elec} / 6 ! (n-1)d10 (n)s2"
                gas2 = "{0} {0} / {1} ! virtual orbitals".format(nb_open_elec, nb_vir_shell // 2)
                gas_list.extend([gas1, gas2])
            else:
                raise RuntimeError(f"Unknown number of open shell: {nb_open_shell}")

        # specify many roots should we compute
        if nb_open_elec == 1:
            if nb_open_shell == 2:
                # only s is open
                ciroot1 = "3 1"
                # ciroot2 = '2 1'
                root_list.extend([ciroot1])
            elif nb_open_shell == 12:
                raise NotImplementedError(
                    "For two open-shell d elements, " "there is no implement!"
                )
            else:
                raise RuntimeError("the number of open shell {} " "is error!".format(nb_open_shell))
        elif nb_open_elec == 0:
            # closed-shell elementss
            ciroot1 = "1 1"
            # ciroot2 = '2 1'
            root_list.extend([ciroot1])
        else:
            raise NotImplementedError(
                "For open-shell d elements and more than "
                "one open-shell element, it has not "
                "been implmented yet!"
            )

        return nb_active_elec, gas_list, root_list

    @staticmethod
    def get_info_for_f_block(out: Output, e_min, e_max, calc_type="quadrupole"):
        """
        Get information for the f block.

        Parameters
        ----------
        out : Output
            The output file to parse.
        e_min : float
            The minimum energy.
        e_max : float
            The maximum energy.
        calc_type : str, optional
            The type of calculation. Default is "quadrupole".

        Returns
        -------
        nb_active_elec: int
            The number of electrons in the GAS.
        gas_list: list
            A list for GAS setup.
        root_list: list
            Root for different symmetries.
        """
        atomic_info = out.mol.molecule.atomic_info
        mos = out.mos
        assert atomic_info.block == "f"

        (
            nb_closed_elec,
            nb_open_elec,
            nb_total_elec,
            nb_closed_shell,
            nb_open_shell,
            nb_vir_shell,
        ) = mos.get_ao_and_elec(e_min, e_max)

        assert nb_closed_elec == nb_closed_shell
        assert nb_open_elec <= nb_open_shell
        period = atomic_info.period

        gas_list = []
        root_list = []

        raise NotImplementedError("For open-shell f elements, it has not " "been implmented yet!")

    @classmethod
    def from_prev_dhf_calc(cls, filename_input, e_min=-10, e_max=10.0):
        """
        Create MRCI input based on the previous SCF calculation.

        Parameters
        ----------
        filename_input : str
            The name of the input file.
        e_min : float, optional
            The minimum energy. Default is -10.
        e_max : float, optional
            The maximum energy. Default is 10.0.

        Returns
        -------
        cls
            The MRCI input object.
        """

        out = Output(filename_input)
        out.parse_orbit()
        atom = out.mol.molecule
        # print(atom.is_atom)

        if not atom.is_atom:
            warnings.warn("This is not an atomic calculation!")
            return

        # this is for p block elements
        if atom.atomic_info.block == "p":
            nb_active_elec, gas_list, root_list = AtomicCISet.get_info_for_p_block(
                out, e_min=e_min, e_max=e_max, calc_type="quadrupole"
            )
        elif atom.atomic_info.block == "s":
            nb_active_elec, gas_list, root_list = AtomicCISet.get_info_for_s_block(
                out, e_min=e_min, e_max=e_max, calc_type="quadrupole"
            )
        elif atom.atomic_info.block == "d":
            nb_active_elec, gas_list, root_list = AtomicCISet.get_info_for_d_block(
                out, e_min=e_min, e_max=e_max, calc_type="quadrupole"
            )
        elif atom.atomic_info.block == "f":
            nb_active_elec, gas_list, root_list = AtomicCISet.get_info_for_f_block(
                out, e_min=e_min, e_max=e_max, calc_type="quadrupole"
            )
        else:
            raise RuntimeError("Unknown block element!")

        # print('total number electrons is {0}'.format(atom.atomic_info.Z))
        inactivate = (atom.atomic_info.Z - nb_active_elec) // 2
        nb_gas_shell = len(gas_list)
        nb_root_sym = len(root_list)

        krci_setup = {
            "CI PROGRAM": "LUCIAREL",
            "INACTIVE": str(inactivate),
            "GAS SHELLS": str(nb_gas_shell) + "\n" + "\n".join(gas_list),
            "MAX CI": 120,
            "MXCIVE": 60,
            "ANALYZ": True,
            "RSTRCI": "rstr",
            "CHECKP": True,
            "NOOCCN": True,
        }
        for i in range(nb_root_sym):
            krci_setup[f"CIROOTS_id_{i}"] = root_list[i]

        inp = out.inp
        if "HAMILTONIAN" in inp and "OPERATOR" in inp["HAMILTONIAN"]:
            inp["HAMILTONIAN"]["OPERATOR"][-1] = " zff"
        inp[inp.wf_tag]["KR CI"] = True
        inp[inp.wf_tag]["KRCICALC"] = krci_setup
        ci_obj = cls(atom)
        ci_obj.prev_inp = inp
        ci_obj.prev_mol = out.mol
        return ci_obj


_dummy_structure = Molecule([1], [[0.0, 0.0, 0.0]])


def get_diracrun_output(path) -> Output:
    """
    Get the output from a previous DIRAC calculation.

    Parameters
    ----------
    path : str
        The path to the previous calculation directory.

    Returns
    -------
    Output
        The output object from the previous calculation.
    """
    path = Path(path)
    outputs = list(glob.glob(str(path / "*.out")))
    if len(outputs) == 0:
        raise ValueError(f"Unable to get *.out from prev calculation in ${path}")
    outfile = sorted(outputs)[-1]
    return Output(outfile)


if __name__ == "__main__":
    import tempfile

    def main2():
        module_dir = os.path.dirname(__file__)
        output = str(Path(module_dir) / "../tests/data/Rn_q_so/d-aug-dyall.acv3z_+0.00001")

        # Use a temporary directory to store the output files
        with tempfile.TemporaryDirectory() as tmp_dir:
            old_cc = AtomicCCSet.from_prev_calc(
                output, no_scf=True, user_mol_settings={"basis_type": "s-aug-ANO-RCC"}
            )
            old_cc.write_input(output_dir=tmp_dir)
            new_cc = AtomicCCSet(
                molecule=Molecule([7], [[0.0, 0.0, 0.0]]),
                hamiltonian_mode="2C",
                is_spinfree=True,
                is_ff=True,
                ff_mode="D",
                user_mol_settings={"basis_type": "s-aug-ANO-RCC"},
            )
            new_cc.write_input(output_dir=tmp_dir)

    main()
    main2()
