import os
import re
import warnings
import subprocess
from monty.json import MSONable, jsanitize, MontyDecoder
from monty.os import cd

from pydirac import SETTINGS
from pydirac.core.settings import Settings
from pydirac.core.molecule import Molecule
from pydirac.core.periodic import Element

__all__ = ["DiracInput", "Inp", "Mol"]


class DiracInput(dict, MSONable):
    """
    Class to contain a set of DIRAC input objects corresponding to a run.

    Parameters
    ----------
    inp : Inp
        Inp object.
    mol : Mol
        Mol object.
    optional_files : dict, optional
        Other input files supplied as a dictionary of {filename: object}. The object
        should follow standard conventions in implementing an `as_dict()` and
        `from_dict()` method. Defaults to None.

    Attributes
    ----------
    keys : List[str]
        List of keys in the DiracInput dictionary.

    """

    def __init__(self, inp, mol, optional_files=None, **kwargs):
        """
        Parameters
        ----------
        inp : Inp
            Inp object.
        mol : Mol
            Mol object.
        optional_files : dict, optional
            Other input files supplied as a dictionary of {filename: object}. The object
            should follow standard conventions in implementing an `as_dict()` and
            `from_dict()` method. Defaults to None.
        **kwargs
            Keyword arguments to be passed to the base class.

        """
        super().__init__(**kwargs)
        self.update({"inp": inp, "mol": mol})
        if optional_files is not None:
            self.update(optional_files)

    def __str__(self):
        """
        Returns a string representation of the DiracInput object.

        Returns
        -------
        str
            String representation of the DiracInput object.

        """
        output = []
        for k, v in self.items():
            output.append(k)
            output.append(str(v))
            output.append("")
        return "\n".join(output)

    def as_dict(self) -> dict:
        """
        Returns a MSONable dict representation of the DiracInput object.

        Returns
        -------
        dict
            MSONable dict representation of the DiracInput object.

        """
        d = {k: v.as_dict() for k, v in self.items()}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Returns a DiracInput object from a dictionary.

        Parameters
        ----------
        d : dict
            Dict representation.

        Returns
        -------
        DiracInput
            DiracInput object.

        """
        dec = MontyDecoder()
        sub_d = {"optional_files": {}}
        for k, v in d.items():
            # TODO: here, inp and mol is not filename but suffix
            if k in ["inp", "mol"]:
                sub_d[k.lower()] = dec.process_decoded(v)
            elif k not in ["@module", "@class"]:
                sub_d["optional_files"][k] = dec.process_decoded(v)
        return cls(**sub_d)

    def write_input(self, output_dir=".", make_dir_if_not_present=True):
        """
        Write DIRAC input to a directory.

        Parameters
        ----------
        output_dir : str, optional
            Directory to write to. Defaults to current directory (".").
        make_dir_if_not_present : bool, optional
            Create the directory if not present. Defaults to True.

        Returns
        -------
        None

        """
        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for k, v in self.items():
            if v is not None:
                with open(os.path.join(output_dir, k), "wt") as f:
                    f.write(v.__str__())

    @staticmethod
    def from_directory(input_dir, optional_files=None):
        """
        Read in a set of DIRAC input from a directory. Note that only the standard
        inp, mol files are read unless optional_files is specified.

        Parameters
        ----------
        input_dir : str
            Directory to read DIRAC input from.
        optional_files : dict, optional
            Optional files to read in as well as a dictionary of {filename: Object type}.
            Object type must have a static method from_file. Defaults to None.

        Returns
        -------
        DiracInput
            DiracInput object.

        """
        sub_d = {}
        for fname, ftype in [("PYDIRAC.inp", Inp), ("PYDIRAC.mol", Mol)]:
            try:
                fullpath = os.path.join(input_dir, fname)
                sub_d[fname.lower()] = ftype.from_file(fullpath)
            except FileNotFoundError:
                sub_d[fname.lower()] = None

        sub_d["optional_files"] = {}
        if optional_files is not None:
            for fname, ftype in optional_files.items():
                sub_d["optional_files"][fname] = ftype.from_file(
                    os.path.join(input_dir, fname)
                )
        return DiracInput(**sub_d)

    def run_dirac(
        self,
        run_dir=".",
        dirac_cmd: list = None,
        output_file="dirac.out",
        err_file="dirac.err",
    ):
        """
        Write input files and run DIRAC.

        Parameters
        ----------
        run_dir : str, optional
            Where to write input files and do the run.
        dirac_cmd : list, optional
            Args to be supplied to run DIRAC.
        output_file : str, optional
            File to write output. Defaults to "dirac.out".
        err_file : str, optional
            File to write error. Defaults to "dirac.err".

        Returns
        -------
        None

        Raises
        ------
        RuntimeError
            If `dirac_cmd` is not supplied or `DIRAC_EXE` is not set in `.pydirac.yaml`.

        """
        self.write_input(output_dir=run_dir)
        dirac_cmd = dirac_cmd or SETTINGS.get("DIRAC_EXE")
        dirac_cmd = [os.path.expanduser(os.path.expandvars(t)) for t in dirac_cmd]
        if not dirac_cmd:
            raise RuntimeError(
                "You need to supply dirac_cmd or set the DIRAC_EXE in .pydirac.yaml to run DIRAC."
            )
        with cd(run_dir):
            with open(output_file, "w") as f_std, open(
                err_file, "w", buffering=1
            ) as f_err:
                subprocess.check_call(dirac_cmd, stdout=f_std, stderr=f_err)


class Inp(dict, MSONable):
    """Class to represent DIRAC input file.

    Parameters
    ----------
    params : dict, optional
        Dictionary containing the input file data.

    Attributes
    ----------
    _allowed_duplicated_list : list
        A list of keywords that are allowed to be duplicated.
    _top : list
        A list of top-level keywords in the input file.

    Methods
    -------
    as_dict() -> dict
        Convert the instance to a dictionary.
    from_dict(d: dict) -> Inp
        Initialize the instance from a dictionary.
    add_sub_node(settinging_obj, dir_node, subdir_node, keyword_node, value_list) -> bool
        Add sub node to parent node.
    from_string(lines: list) -> Inp
        Initialize the instance from a list of strings.
    get_string() -> str
        Transform all contents of ``input`` branch of ``settings`` into string with blocks,
        subblocks, keys and values.
    __str__() -> str
        Return the string representation of the instance.

    Properties
    ----------
    use_wavefunc : bool
        True if wavefunction is to be used in the calculation, False otherwise.
    has_hamiltonian : bool
        True if Hamiltonian is present in the input file, False otherwise.
    hamiltonian : dict
        Dictionary representing the Hamiltonian data in the input file.
    calc_type : str
        The calculation type (e.g., 'energy', 'gradient', 'hessian', etc.) in the input file.
    calc_method : str
        The calculation method (e.g., 'HF', 'DFT', 'MP2', etc.) in the input file.
    electric_field : list
        List of electric field strengths and orientations in the input file.
    ci_calc_orbit : dict
        Dictionary representing the CI calculation orbit in the input file.
    wf_tag : str
        The wavefunction tag in the input file.

    Raises
    ------
    AssertionError
        If the directory node is None or if the value list is not a list.
        RuntimeError: If a duplicated key is encountered.

    Examples
    --------
    >>> inp = Inp({'dirac': {}})
    >>> inp
    {'dirac': {}}

    Notes
    -----
    The Inp class is used to represent DIRAC input files.

    """

    _allowed_duplicated_list = ["CIROOTS"]
    _top = ["dirac", "analyze", "hamiltonian", "integrals", "general"]

    def __init__(self, params=None):
        """Initialize the Inp class instance."""
        super(Inp, self).__init__()
        if params:
            pass
        self.update(params)
        self._parse_input()

    def as_dict(self) -> dict:
        """Convert the instance to a dictionary."""
        d = dict(self)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        """Initialize the instance from a dictionary."""
        return Inp({k: v for k, v in d.items() if k not in ("@module", "@class")})

    @staticmethod
    def add_sub_node(settinging_obj, dir_node, subdir_node, keyword_node, value_list):
        """Add sub node to parent node."""
        assert type(dir_node) is not None
        assert type(value_list) == list

        # whether add value to parent node
        # if yes: value_list has been added to parent node
        # if no: nothing has been done
        is_set = False

        if subdir_node is not None and keyword_node is not None:
            # -------------------------------------------------------
            # Case 1, for example:
            # -------------------------------------------------------
            # **ANALYZE
            # .MULPOP
            # *MULPOP
            # .VECPOP
            # 1..oo
            if len(value_list) > 0:
                settinging_obj[dir_node][subdir_node][keyword_node] = value_list
            else:
                settinging_obj[dir_node][subdir_node][keyword_node] = True
            is_set = True
        elif subdir_node is None and keyword_node is not None:
            # -------------------------------------------------------
            # Case 2, for example:
            # -------------------------------------------------------
            # **DIRAC
            # .TITLE
            # B, DOSSSS, KRCI
            if len(value_list) > 0:
                settinging_obj[dir_node][keyword_node] = value_list
            else:
                settinging_obj[dir_node][keyword_node] = True
            is_set = True
        elif subdir_node is not None and keyword_node is None:
            # if subdir_node is not None and keyword_node is None
            # -------------------------------------------------------
            # Case 3, for example:
            # -------------------------------------------------------
            # **INTEGRALS
            # *READINP
            # .UNCONTRACT
            #
            # when it's processing .UNCONTRACT, it's the case 3
            pass
        else:
            # if subdir_node is None and keyword_node is None
            # -------------------------------------------------------
            # Case 4, for example:
            # -------------------------------------------------------
            # **DIRAC
            #
            # when processing the part in which there is only directory
            pass
        return is_set

    @classmethod
    def from_string(cls, lines):
        """
        Initialize the instance from a list of strings.

        Parameters
        ----------
        lines : list
            List of strings representing the input file.

        Returns
        -------
        Inp
            Inp instance initialized from the input file.
        """
        setting = {}

        old_dir = None
        old_subdir = None
        old_dotkey = None
        curr_value_list = []

        for line in lines:
            if line.startswith("#"):
                # this is a comment line
                continue

            if line.startswith("**"):
                if old_dir is not None:
                    is_set = Inp.add_sub_node(
                        setting, old_dir, old_subdir, old_dotkey, curr_value_list
                    )
                    if is_set:
                        curr_value_list = []
                if len(curr_value_list) > 0:
                    Inp.add_sub_node(
                        setting, old_dir, old_subdir, old_dotkey, curr_value_list
                    )

                # this is a directory
                dir_name = line.lstrip("**").rstrip()
                old_dir = dir_name
                # setting[old_dir] = Settings()
                setting[old_dir] = {}

                ## clear current value for subdir and keyword
                old_subdir = None
                old_dotkey = None
                curr_value_list.clear()

            elif line.startswith("*"):
                # this is a sub-directory
                if old_dir is not None:
                    is_set = Inp.add_sub_node(
                        setting, old_dir, old_subdir, old_dotkey, curr_value_list
                    )
                    if is_set:
                        curr_value_list = []

                if line.startswith("*END"):
                    continue

                sub_dir_name = line.lstrip("*").rstrip()
                old_subdir = sub_dir_name

                assert old_dir is not None
                if sub_dir_name in setting[old_dir]:
                    ## if sub_dir_name in old_dir word list
                    setting[old_dir][sub_dir_name] = {}
                    setting[old_dir][sub_dir_name]["_en"] = True
                else:
                    ## this is a real second level directory
                    setting[old_dir][old_subdir] = {}

                ## clear current value for keyword
                old_dotkey = None
                curr_value_list = []

            elif line.startswith("."):
                # this is dot keyword
                curr_dotkey = line.lstrip(".").rstrip()
                if old_dir is not None:
                    ## first we probabily have two same dotkey, e.g., CIROOTS
                    is_set = Inp.add_sub_node(
                        setting, old_dir, old_subdir, old_dotkey, curr_value_list
                    )
                    if is_set:
                        curr_value_list = []

                    if (
                        old_subdir is not None
                        and curr_dotkey in setting[old_dir][old_subdir].keys()
                        and curr_dotkey.upper() in Inp._allowed_duplicated_list
                    ):
                        # first we need to clear curr_value_list and add this to the
                        # previous one
                        if is_set:
                            # curr_value_list = setting[old_dir][old_subdir][
                            #     curr_dotkey]
                            # curr_value_list.append('.' + curr_dotkey)

                            ## get id:
                            _tmp_k = []
                            for _k in setting[old_dir][old_subdir].keys():
                                if curr_dotkey == _k:
                                    continue
                                if curr_dotkey in _k:
                                    _tmp_k.append(int(_k.split("_id_")[-1]))
                            if len(_tmp_k):
                                new_id = max(_tmp_k) + 1
                            else:
                                new_id = 1
                            curr_dotkey = curr_dotkey + "_id_" + str(new_id)

                        old_dotkey = curr_dotkey
                        continue

                    # if old_subdir is None and curr_dotkey in setting[old_dir]:
                    #     print('{}, {}, {}'.format(old_dir, old_subdir,
                    #                               curr_dotkey))
                    #     curr_value_list.append('.' + curr_dotkey)
                    #     continue

                old_dotkey = curr_dotkey
                assert not (old_dir is None and old_subdir is None)

                ## clear value for value list
                curr_value_list = []
            else:
                # this is value of keyword
                curr_value = line.rstrip()
                if len(curr_value):
                    curr_value_list.append(curr_value)
        return Inp.from_dict(setting)

    @classmethod
    def from_file(cls, filename):
        """
        Parse DIRAC input file to restore python Settings object
        """
        with open(filename, "r") as f:
            lines = f.readlines()
        return Inp.from_string(lines)

    def get_string(self):
        """Transform all contents of ``input`` branch of ``settings``
        into string with blocks, subblocks, keys and values.

        On the highest level alphabetic order of iteration is modified:
        keys occuring in class attribute ``_top`` are printed first.
        See :ref:`dirac-input` for details.
        """
        is_empty = lambda x: isinstance(x, dict) and len(x) == 0

        def parse_key(key, value):
            id_ptn = re.compile("^(.*)_id_\d+$")
            m = id_ptn.match(key)
            if m:
                key = m.group(1)
                if key.upper() not in self._allowed_duplicated_list:
                    raise RuntimeError("Duplicated key does not support!")
            if value is False:
                return ""
            ret = "." + key.upper() + "\n"
            if not (value is True or is_empty(value)):
                if isinstance(value, list):
                    for i in value:
                        ret += str(i) + "\n"
                else:
                    ret += str(value) + "\n"
            return ret

        def parse_block(block):
            enabler = "_en"
            ret = "**" + block.upper() + "\n"
            s = self[block]
            for k, v in s.items():
                if not isinstance(v, dict) or is_empty(v):
                    ret += parse_key(k, v)
            for k, v in s.items():
                if isinstance(v, dict) and enabler in v:
                    ret += parse_key(k, v[enabler])
            for k, v in s.items():
                if isinstance(v, dict) and len(v) > 0:
                    ret += "*" + k.upper() + "\n"
                    for kk, vv in v.items():
                        if kk != enabler:
                            ret += parse_key(kk, vv)
            return ret

        inp = ""
        for block in self._top:
            if block in self:
                inp += parse_block(block)
        for block in self:
            if block.lower() not in self._top:
                inp += parse_block(block)
        inp += "*END OF INPUT\n"
        return inp

    def __str__(self):
        """
        Return a string representation of the object.

        Returns
        -------
        str
            The string representation of the object.
        """
        return self.get_string()

    @property
    def use_wavefunc(self):
        """
        Return whether or not to use a wave function.

        Returns
        -------
        bool
            Whether or not to use a wave function.
        """
        return self._use_wavefunc

    @property
    def has_hamiltonian(self):
        """
        Return whether or not a Hamiltonian is present.

        Returns
        -------
        bool
            Whether or not a Hamiltonian is present.
        """
        return self._has_hamiltonian

    @property
    def hamiltonian(self):
        """
        Return the Hamiltonian.

        Returns
        -------
        str
            The Hamiltonian.
        """
        return self._hamiltonian

    @property
    def calc_type(self):
        """
        Return the calculation type.

        Returns
        -------
        str
            The calculation type.
        """
        return self._calc_type

    @property
    def calc_method(self):
        """
        Return the calculation method.

        Returns
        -------
        str
            The calculation method.
        """
        return self._calc_method

    @property
    def electric_field(self):
        """
        Return the strength of the electric field.

        Returns
        -------
        float
            The strength of the electric field.
        """
        return self._electric_field

    @property
    def ci_calc_orbit(self):
        """
        Return the CI calculation orbit.

        Returns
        -------
        dict
            The CI calculation orbit.
        """
        if hasattr(self, "_ci_calc_orbit"):
            return self._ci_calc_orbit
        return {"occ": 0, "vir": 0}

    def _find_real_key(self, k, keys):
        """
        Find the real key from a list of keys.

        Parameters
        ----------
        k : str
            The key to find.
        keys : list
            The list of keys to search.

        Returns
        -------
        str or None
            The real key if found, otherwise None.
        """
        for _k in keys:
            if k in _k:
                return _k
        else:
            return None

    @property
    def wf_tag(self):
        """
        Return the wave function tag.

        Returns:
        str
             the wave function tag.
        """
        return self._wf_tag

    def _parse_input(self):
        """Parse input from output file

        Returns
        -------
            input settings
        """

        inp_settings = Settings(self)

        self._wf_tag = None
        for dk in inp_settings.dirac:
            if "WAVE F" in dk:
                self._use_wavefunc = True
                self._wf_tag = self._find_real_key("WAVE F", self)
                break
        else:
            self._use_wavefunc = False
            self._has_hamiltonian = False
            self._hamiltonian = None
            self._calc_type = None
            self._calc_method = None
            self._electric_field = None
            return

        # check Hamiltonian
        if "HAMILTONIAN" in inp_settings:
            self._has_hamiltonian = True
            if "DOSSSS" in inp_settings.hamiltonian:
                # this is 4-component Hamiltonian
                self._hamiltonian = "4C"
            elif "X2C" in inp_settings.hamiltonian:
                # this is 2-component Hamiltonian
                self._hamiltonian = "2C"
            elif "X2Cmmf" in inp_settings.hamiltonian:
                # need to do
                self._hamiltonian = "2C"
            elif "NONREL" in inp_settings.hamiltonian:
                # self._hamiltonian = "NR-2C"
                self._hamiltonian = "2C-NR"
            else:
                warnings.warn(
                    "We do not know how many components, " "please check inp file"
                )
                self._hamiltonian = "?C"

            if not "NONREL" in inp_settings.hamiltonian:
                if "NOSPIN" in inp_settings.hamiltonian:
                    self._hamiltonian += "-SR"
                else:
                    self._hamiltonian += "-DC"

            if "Gaunt" in inp_settings.hamiltonian:
                self._hamiltonian = self._hamiltonian + "-Gaunt"
        else:
            self._has_hamiltonian = False

        # check if CC or CI calculation
        if (
            self._has_hamiltonian
            and "RELCC" in inp_settings
            and "RELCCSD" in inp_settings[self._wf_tag]
        ):
            self._calc_method = "CC"
        elif (
            self._has_hamiltonian
            and "KRCICALC" in inp_settings[self._wf_tag]
            and "KR CI" in inp_settings[self._wf_tag]
        ):
            self._calc_method = "CI"
            # here, we can obtain calc_orbit info from input directly
            gas_tag = None
            for k in inp_settings[self._wf_tag]["KRCICALC"].keys():
                if "GAS" in k:
                    gas_tag = k
                    break
            gas_res = inp_settings[self._wf_tag]["KRCICALC"][gas_tag]
            # the first line is the number of GAS shell
            gas_nb = int(gas_res[0])
            assert gas_nb == len(gas_res[1:])

            nb_e = 0
            t_orb = []
            gas_ptn = re.compile(r"\s*\d+\s+(?P<nb_e>\d+)\s*/\s*(?P<nb_orb>\d+)")
            for l in gas_res[1:]:
                m = gas_ptn.match(l)
                if m:
                    t_orb.append(int(m.group("nb_orb")))
            m = gas_ptn.match(gas_res[-1])
            if m:
                nb_e = int(m.group("nb_e"))
                tot_orb = sum(t_orb)
                occ_orb = nb_e
                vir_orb = tot_orb * 2 - occ_orb
                self._ci_calc_orbit = {"occ": occ_orb, "vir": vir_orb}
            else:
                self._ci_calc_orbit = {"occ": 0, "vir": 0}

        elif self._has_hamiltonian and "SCF" in inp_settings[self._wf_tag]:
            if inp_settings[self._wf_tag]["SCF"]["_en"]:
                self._calc_method = "SCF"
            else:
                self._calc_method = None
        # TODO: visual module calculation
        else:
            self._calc_method = None

        # check dipole or quadrupole
        if self._has_hamiltonian and "OPERATOR" in inp_settings.hamiltonian:
            if " ZZTHETA" in inp_settings.hamiltonian.operator:
                self._calc_type = "Q"
            elif " ZDIPLEN" in inp_settings.hamiltonian.operator:
                self._calc_type = "D"
            else:
                warnings.warn("this is not finite-field calculations")
                self._calc_type = "non_finite_field"
        else:
            if not self._has_hamiltonian:
                warnings.warn("Hamiltonian is invalid")
            elif not "OPERATOR" in inp_settings.hamiltonian:
                # warnings.warn('".OPERATOR" is not in Hamiltonian, this is not finit-field calculation!')
                pass
            self._calc_type = "non_finite_field"

        # check electric field
        if self._calc_type in "QD":
            # this is finit-field
            pos = inp_settings.hamiltonian.operator.index(" COMFACTOR")
            self._electric_field = inp_settings.hamiltonian.operator[pos + 1]
        else:
            # warnings.warn('Did not find electric field')
            self._electric_field = "null"

    def find_case(self, key):
        """Check if this instance contains a key consisting of the same
        letters as *key*, but possibly with different case. If found, return
        such a key. If not, return *key*.

        Args:
            key:

        Returns:

        """
        if not isinstance(key, str):
            return key
        lowkey = key.lower()
        for k in self:
            if k.lower() == lowkey:
                return k
        return key

    def __contains__(self, name):
        """Like regular ``__contains`__``, but ignore the case."""
        return dict.__contains__(self, self.find_case(name))

    def __getitem__(self, name):
        """Like regular ``__getitem__``, but ignore the case."""
        return dict.__getitem__(self, self.find_case(name))

    def __setitem__(self, name, value):
        """Like regular ``__setitem__``, but ignore the case and if the
        value is a dict, convert it to |Settings|.
        """
        dict.__setitem__(self, self.find_case(name), value)

    def __delitem__(self, name):
        """Like regular ``__detitem__``, but ignore the case."""
        return dict.__delitem__(self, self.find_case(name))

    def write_file(self, filename):
        with open(filename, "wt") as f:
            f.write(self.__str__())


class Mol(MSONable):
    """Class to represent DIRAC Mol file"""

    _basis_info_line = 2

    def __init__(
        self,
        molecule: Molecule,
        comment: str = None,
        basis_type: str = None,
        basis_lib: str = "BASIS",
    ):
        self.molecule = molecule
        self.comment = comment
        self.basis_type = basis_type
        self.basis_lib = basis_lib

    def as_dict(self) -> dict:
        """Return object as dict."""
        d = {
            "molecule": self.molecule.as_dict(),
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
        }

        return d

    @classmethod
    def from_file(cls, filename):
        """Initialize from file

        Parameters
        ----------
        filename: str
            Filename
        """
        # TODO, more info from mol file including basis sets.
        with open(filename, "r") as f:
            lines = f.readlines()
        return Mol.from_string(lines)

    @classmethod
    def from_string(cls, context):
        """
        Construct a `Mol` object by parsing the contents of a MOL file.

        Parameters
        ----------
        context : list of str
            The contents of the MOL file as a list of strings.

        Returns
        -------
        Mol
            A `Mol` object representing the molecule.

        """

        start_line = re.compile(r"^Contents of the molecule file\s*")
        dash_line = re.compile(r"^-+")
        empty_line = re.compile(r"^$")
        # inp_start_line = re.compile(r"^\*\*DIRAC\s*")
        nuclei_nb_line = re.compile(
            r"^\s*(?P<nuclei>[-+]?(\d+(\.\d*)?|\d*\.\d+))\s+(?P<nb_atoms>\d+)\s*$"
        )
        basis_line = re.compile(r"^\b(?:LARGE|EXPLICIT)\b\s+BASIS\s+(\S+)\s*")
        endline = re.compile(r"^FINISH")

        nuclei_id = 0
        nb_atoms = 0
        basis_type = None

        for i, line in enumerate(context):
            if endline.match(line):
                break
            elif dash_line.match(line) or empty_line.match(line):
                continue
            elif nuclei_nb_line.match(line):
                m = nuclei_nb_line.match(line)
                nuclei_id = int(round(float(m.group("nuclei"))))
                nb_atoms = int(m.group("nb_atoms"))
                # print(self.nuclei_id, self.nb_atoms)
            elif basis_line.match(line):
                m = basis_line.match(line)
                if m:
                    basis_type = m.group(1)
                    break
            else:
                # we do not interest with these lines
                continue

        if basis_type is None:
            basis_type = context[Mol._basis_info_line].strip()
        basis_type = basis_type
        # TODO: we are doing atomic calculation
        assert nb_atoms == 1
        molecule = Molecule(atoms=[nuclei_id], coordinates=[[0.0, 0.0, 0.0]])
        return cls(basis_type=basis_type, molecule=molecule)

    @staticmethod
    def get_string(
        atom_info: Element,
        basis_type: str,
        basis_lib: str = "BASIS",
    ) -> str:
        """
        Return a string representation of a molecule in MOL format with the specified
        atom information and basis.

        Parameters
        ----------
        atom_info : Element
            The atomic information for the molecule.
        basis_type : str
            The type of basis.
        basis_lib : str, optional
            The type of basis library. Should be either "BASIS" or "EXPLICIT".

        Returns
        -------
        str
            A string representation of the molecule in MOL format.

        Raises
        ------
        TypeError
            If the `basis_lib` parameter is not "BASIS" or "EXPLICIT".

        """
        atom_index = atom_info.Z
        atom_type = atom_info.symbol

        if basis_lib not in ["EXPLICIT", "BASIS"]:
            raise TypeError(
                'Basis type should be "BASIS" or "EXPLICIT" '
                "for builtin basis or custom basis."
            )

        if basis_lib == "EXPLICIT":
            with open("basis/{0}.dat".format(atom_type), "r") as f:
                basis_info = f.read()
            template = Mol.get_mol_by_custom_basis(
                atom_type, atom_index, basis_lib, basis_info
            )
        elif basis_lib == "BASIS":
            template = Mol.get_mol_by_default_basis(atom_type, atom_index, basis_type)
        else:
            raise TypeError(
                'Basis type should be "BASIS" or "EXPLICIT" '
                "for builtin basis or custom basis."
            )
        return template

    @staticmethod
    def get_mol_by_custom_basis(atom_type, atom_index, basis_choice, basis_info):
        """
        Return a string representation of a molecule in MOL format with a custom basis.

        Parameters
        ----------
        atom_type : str
            The type of atom.
        atom_index : int
            The index of the atom.
        basis_choice : str
            The type of basis.
        basis_info : str
            The custom basis information.

        Returns
        -------
        str
            A string representation of the molecule in MOL format with the custom basis.

        """
        template = """INTGRL
{atom} atom
{basis}
C   1    2  X  Y
    {elec}.    1
{atom}    0.0000000000        0.0000000000        0.0000000000
{basis_info}
""".format(
            **{
                "atom": atom_type,
                "elec": atom_index,
                "basis": basis_choice,
                "basis_info": basis_info,
            }
        )
        return template

    @staticmethod
    def get_mol_by_default_basis(atom_type, atom_index, basis_type):
        """
        Return a string representation of a molecule in MOL format, with the
        specified atom type, atom index, and basis type.

        Parameters
        ----------
        atom_type : str
            The type of atom.
        atom_index : int
            The index of the atom.
        basis_type : str
            The type of basis.

        Returns
        -------
        str
            A string representation of the molecule in MOL format.

        """
        template = """INTGRL
{atom} atom
{basis}
C   1    2  X  Y
    {elec}.    1
{atom}    0.0000000000        0.0000000000        0.0000000000
LARGE BASIS {basis}
FINISH

""".format(
            **{"atom": atom_type, "elec": atom_index, "basis": basis_type}
        )
        return template

    def __str__(self):
        """
        Return a string representation of the molecule.

        Returns
        -------
        str
            A string representation of the molecule.

        """

        basis_type = self.basis_type or "null"
        if self.basis_lib == "EXPLICIT":
            warnings.warn("No implement {0} yet".format(self.basis_lib))
            pass
        return Mol.get_string(self.molecule.atomic_info, basis_type)

    def write_file(self, filename=None):
        """
        Write the molecule to a file in MOL format.

        Parameters
        ----------
        filename : str, optional
            The name of the file to write to. If not provided, the file will be
            named after the atom symbol and basis type.

        Returns
        -------
        None

        Raises
        ------
        NotImplementedError
            If the molecule contains more than one atom.

        """
        if self.molecule.is_atom:
            fname = (
                filename
                or self.molecule.atomic_info.symbol + "_" + self.basis_type + ".mol"
            )
            with open(fname, "w") as f:
                basis_type = self.basis_type or "null"
                f.write(Mol.get_string(self.molecule.atomic_info, basis_type))
        else:
            raise NotImplementedError(
                "For many-atoms molecule, " "it has not been implemented yet! "
            )

    def update(self, new_dict):
        """
        Update the attributes of the molecule with the values in `new_dict`.

        Parameters
        ----------
        new_dict : dict
            A dictionary of attribute-value pairs to update.

        Returns
        -------
        None

        """
        for k, v in new_dict.items():
            setattr(self, k, v)
