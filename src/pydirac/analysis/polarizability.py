import glob
import os
import warnings

import numpy as np
from scipy.linalg import lstsq
from scipy.special import factorial
from monty.os import cd

from pydirac.analysis.utility import get_keyword, get_orbital_info
from pydirac.io.outputs import Output

__all__ = [
    "FiniteFieldLstsqSolver",
    "PolarizabilityCalculator",
    "check_results",
    "get_polarizability",
    "get_polarizability_from_output_list",
    "do_one_basis",
    "is_valid",
]


def check_results(x_svd, errors, nb_coeffs):
    """
    Check the results of the least squares fitting for validity.

    Parameters
    ----------
    x_svd : ndarray
        Solution to the least squares problem.
    errors : ndarray
        Standard errors of the solution.
    nb_coeffs : int
        The number of fitted coefficients required.

    Returns
    -------
    bool
        True if the results are valid, False otherwise.
    """
    if nb_coeffs > 1:
        if np.isclose(x_svd[1], 0.0):
            return np.allclose(x_svd[1:], 0.0) and np.allclose(errors[1:], 0.0)
        if x_svd[1] < 0 or errors[1] / x_svd[1] < 1:
            return False
    return True


class FiniteFieldLstsqSolver:
    """
    A class to solve the least squares problem for polarizability calculations using finite field method.
    """

    def __init__(
        self, energy, field, nb_coeffs, rcond, threshold, calc_type, check_hyper_polar
    ):
        """
        Initialize the FiniteFieldLstsqSolver instance.

        Parameters
        ----------
        energy : array_like
            Array of energy values corresponding to different field strengths.
        field : array_like
            Array of field strengths.
        nb_coeffs : int
            Number of coefficients to use for the least squares fitting. Must be between 1 and 3.
        rcond : float
            Cut-off ratio for small singular values in the least squares solution.
        threshold : float
            Threshold value for field strength to include in the fitting.
        calc_type : {'dipole', 'quadrupole'}
            Type of polarizability calculation to perform.
        check_hyper_polar : bool
            Whether to check for hyperpolarizability during fitting.
        """
        assert calc_type in ["dipole", "quadrupole"]
        self.energy = np.asarray(energy)
        self.field = np.asarray(field)
        self.nb_coeffs = nb_coeffs
        self.rcond = rcond
        self.threshold = threshold
        self.calc_type = calc_type
        self.check_hyper_polar = check_hyper_polar

    def calc_dipole_mat(self, fields):
        """
        Calculate the matrix for dipole polarizability fitting.

        Parameters
        ----------
        fields : array_like
            Array of field strengths.

        Returns
        -------
        ndarray
            The matrix for dipole polarizability fitting.
        """
        f_list = [
            -1 / factorial(2 * i + 2) * fields ** (2 * i + 2)
            for i in range(self.nb_coeffs)
        ]
        return np.column_stack(f_list)

    def calc_quadrupole_mat(self, fields):
        """
        Calculate the matrix for quadrupole polarizability fitting.

        Parameters
        ----------
        fields : array_like
            Array of field strengths.

        Returns
        -------
        ndarray
            The matrix for quadrupole polarizability fitting.
        """
        if self.nb_coeffs == 1:
            f_list = [1 / 2 * fields]
        elif self.nb_coeffs == 2:
            f_list = [1 / 2 * fields, -1 / 8 * fields**2]
        elif self.nb_coeffs == 3:
            f_list = [1 / 2 * fields, -1 / 8 * fields**2, 1 / 24 * fields**3]
        else:
            raise RuntimeError(
                f"nb_coeffs={self.nb_coeffs} should be less than 4 and larger than 0"
            )
        return np.column_stack(f_list)

    def calc_mat(self, fields):
        """
        Calculate the fitting matrix based on the calculation type.

        Parameters
        ----------
        fields : array_like
            Array of field strengths.

        Returns
        -------
        ndarray
            The fitting matrix.
        """
        if self.calc_type == "dipole":
            return self.calc_dipole_mat(fields)
        elif self.calc_type == "quadrupole":
            return self.calc_quadrupole_mat(fields)
        else:
            raise NotImplementedError()

    def preprocess_data(self):
        """
        Preprocess the data by filtering and sorting the field and energy arrays.

        Returns
        -------
        sorted_field : ndarray
            Sorted array of filtered field strengths.
        sorted_energy : ndarray
            Sorted array of filtered energy values.
        e_ref : float
            Reference energy value at zero field.

        Raises
        ------
        ValueError
            If no field value is close to zero.
        """
        zero_field_mask = np.isclose(self.field, 0.0)
        if np.any(zero_field_mask):
            e_ref = self.energy[zero_field_mask][0]
        else:
            raise ValueError("No field value is close to zero.")
        mask = (self.field != 0.0) & (self.field <= self.threshold + 1e-8)
        reduced_field = self.field[mask]
        reduced_energy = self.energy[mask]
        sort_indices = np.argsort(reduced_field)
        sorted_field = reduced_field[sort_indices]
        sorted_energy = reduced_energy[sort_indices]
        return sorted_field, sorted_energy, e_ref

    def solve_least_squares(self, A, b):
        """
        Solve the least squares problem and return the solution, residual, and normalization factors.

        Parameters
        ----------
        A : ndarray
            The fitting matrix.
        b : ndarray
            The energy difference vector.

        Returns
        -------
        x_svd : ndarray
            Solution to the least squares problem.
        residual : float
            The residual sum of squares.
        norms : ndarray
            Normalization factors for the columns of the fitting matrix.
        """
        norms = np.linalg.norm(A, axis=0)
        A_normalized = A / norms
        x_svd_normalized, residuals, _, _ = lstsq(A_normalized, b, cond=self.rcond)

        if A.shape[0] > A.shape[1]:
            assert np.ndim(b) == 1
            assert np.isscalar(residuals)
            residual = residuals
        else:
            assert A.shape[1] == b.shape[0]
            # an empty array is returned
            assert len(residuals) == 0
            residual = 0.0
        x_svd = x_svd_normalized / norms

        ## Old implementation
        # self.rcond = None
        # norms = 1
        # A_normalized = A / norms
        # x_svd_normalized, residuals, _, _ = np.linalg.lstsq(
        #     A_normalized, b, rcond=self.rcond
        # )
        # if len(residuals) == 1:
        #     residual = residuals[0]
        # elif len(residuals) == 0:
        #     residual = 0
        # x_svd = x_svd_normalized / norms
        return x_svd, residual, norms

    def calculate_covariance_matrix(self, A, residual, degrees_of_freedom):
        """
        Calculate the covariance matrix for the least squares solution.

        Parameters
        ----------
        A : ndarray
            The fitting matrix.
        residual : float
            The residual sum of squares.
        degrees_of_freedom : int
            Degrees of freedom for the fitting.

        Returns
        -------
        ndarray
            The covariance matrix.
        """
        if degrees_of_freedom == 0:
            cov_matrix = np.zeros((A.shape[1], A.shape[1]))
        else:
            MSE = residual / degrees_of_freedom
            cov_matrix = np.linalg.inv(np.dot(A.T, A)) * MSE
        return cov_matrix

    def check_results(self, x_svd, errors):
        """
        Check the results of the least squares fitting for validity.

        Parameters
        ----------
        x_svd : ndarray
            Solution to the least squares problem.
        errors : ndarray
            Standard errors of the solution.

        Returns
        -------
        bool
            True if the results are valid, False otherwise.
        """
        return check_results(x_svd, errors, self.nb_coeffs)

    def fit(self):
        """
        Perform the fitting process and return the fitted coefficients and their standard errors.

        Returns
        -------
        x_svd : ndarray
            Fitted coefficients.
        standard_errors : ndarray
            Standard errors of the fitted coefficients.

        Raises
        ------
        RuntimeError
            If an unexpected condition occurs during fitting.
        """
        sorted_field, sorted_energy, e_ref = self.preprocess_data()
        A = self.calc_mat(sorted_field)
        b = sorted_energy - e_ref
        x_svd, residual, _ = self.solve_least_squares(A, b)
        degrees_of_freedom = len(sorted_energy) - len(x_svd)
        cov_matrix = self.calculate_covariance_matrix(A, residual, degrees_of_freedom)
        standard_errors = np.sqrt(np.diag(cov_matrix))

        if self.check_hyper_polar:
            if self.check_results(x_svd, standard_errors):
                return x_svd, standard_errors
            elif self.nb_coeffs > 1:
                x_svd, standard_errors = FiniteFieldLstsqSolver(
                    self.energy,
                    self.field,
                    1,
                    self.rcond,
                    self.threshold,
                    self.calc_type,
                    self.check_hyper_polar,
                ).fit()
                tmp_svd = np.zeros((self.nb_coeffs,))
                tmp_errors = np.zeros((self.nb_coeffs,))
                tmp_svd[0] = x_svd[0]
                tmp_errors[0] = standard_errors[0]
                return tmp_svd, tmp_errors
            else:
                raise RuntimeError("You found a bug!")
        else:
            return x_svd, standard_errors


class PolarizabilityCalculator:
    """
    A class to calculate dipole polarizability and quadrupole polarizability using finite difference
    method.

    Parameters
    ----------
    calc_type : str, optional
        The type of calculation to perform. Must be one of "dipole" or "quadrupole".
        Default is "dipole".
    nb_coeffs : int, optional
        The number of coefficients to use for the calculation. Must be 2 or 3.
        Default is 2.

    Attributes
    ----------
    calc_type : str
        The type of calculation to perform. Must be one of "dipole" or "quadrupole".
    nb_coeffs : int
        The number of coefficients to use for the calculation.
    """

    def __init__(self, calc_type="dipole", nb_coeffs=2):
        """
        Initialize the PolarizabilityCalculator instance.

        Parameters
        ----------
        calc_type : str, optional
            The type of calculation to perform. Must be one of "dipole" or "quadrupole".
            Default is "dipole".
        nb_coeffs : int, optional
            The number of coefficients to use for the calculation. Must be 2 or 3.
            Default is 2.
        """
        _ct = calc_type.lower()
        if _ct in ["d", "dipole"]:
            _ct = "dipole"
        elif _ct in ["q", "quad", "quadrupole"]:
            _ct = "quadrupole"
        else:
            raise ValueError("Wrong calculation type!")
        self.calc_type = _ct

        if nb_coeffs not in (1, 2, 3):
            raise ValueError("nb_coeffs must be 1, 2 or 3!")
        self.nb_coeffs = nb_coeffs

    def get_svd_from_array(
        self, energy, field, rcond=1e-10, threshold=None, check_hyper_polar=False
    ):
        return FiniteFieldLstsqSolver(
            energy,
            field,
            self.nb_coeffs,
            rcond,
            threshold,
            self.calc_type,
            check_hyper_polar,
        ).fit()


def get_polarizability(
    dirname: str = "./",
    calc_dir_patters=None,
    deepth: int = 0,
    verbos=False,
    threshold=None,
) -> dict:
    """Get polarizability from a directory

    A calculation directory may contain some information as followings:
        (1) no more calculation directory (calc_dir_list) but several output files (curr_dir_output_list)
        (2) several output files (curr_dir_output_list) with other calculation directories (calc_dir_list)
        (3) several output files (curr_dir_output_list) with subdirectories (sub_dir_list) --> call self again when deepth > 0
        (4) only calculations directories (calc_dir_list)
        (5) calculations directories (calc_dir_list) and (sub_dir_list)
        (6) only subdirectories (sub_dir_list)
        (7) all (calc_dir_list), (curr_dir_output_list) and (sub_dir_list)

    Args:
        dirname: dirname
        suffix:

    Returns:
        None
    """
    # deepth = 0: only calc current output files
    if calc_dir_patters is None:
        calc_dir_patters = ["*dyall*"]

    if isinstance(calc_dir_patters, (str, list)):
        if isinstance(calc_dir_patters, str):
            calc_dir_patters = [calc_dir_patters]
    else:
        raise TypeError("calc_dir_patters can only be <str> or <list>")

    deepth = 0 if deepth < 0 else deepth
    do_sub_dir = deepth > 0

    # for debug
    if verbos:
        print(f"cd {dirname}")

    all_res = {}
    curr_dir_output_lis = []
    calc_dir_output_lis = []
    with cd(dirname):
        # Step 1. deal with all current output files
        # ------------------------------
        for f in glob.glob("*.out"):
            obj = Output(filename=f)
            if obj.is_ok:
                curr_dir_output_lis.append(obj)

        # Step 2. deal with all calc_dir
        # ------------------------------
        clc_dirs = []
        for ptn in calc_dir_patters:
            clc_dirs.extend([d for d in glob.glob("*" + ptn + "*") if os.path.isdir(d)])
        # if has calc_dir then I would try
        for clc_d in clc_dirs:
            with cd(clc_d):
                outs = glob.glob("*.out")
                if len(outs) > 1:
                    warnings.warn(
                        "there are more than two output file in this directory {}, and we take "
                        "the first one {}".format(clc_d, outs[0])
                    )
                elif len(outs) == 1:
                    obj = Output(filename=outs[0])
                    if obj.is_ok:
                        calc_dir_output_lis.append(obj)
                else:
                    warnings.warn(f"there is no output file in this directory {clc_d}")

        all_res["curr_dir"] = {}
        if is_valid(curr_dir_output_lis):
            e_dict = get_polarizability_from_output_list(
                dirname,
                curr_dir_output_lis,
                tag="curr_dir",
                verbos=verbos,
                threshold=threshold,
            )
            all_res["curr_dir"] = e_dict

        all_res["calc_dir"] = {}
        if is_valid(calc_dir_output_lis):
            e_dict = get_polarizability_from_output_list(
                dirname,
                calc_dir_output_lis,
                tag="clc_dir",
                verbos=verbos,
                threshold=threshold,
            )
            all_res["calc_dir"] = e_dict

        # Step 3. deal with all sub_dir
        # ------------------------------
        all_res["sub_dir"] = {}
        if do_sub_dir:
            sub_dirs = [
                d for d in glob.glob("*") if os.path.isdir(d) and d not in clc_dirs
            ]

            for sd in sub_dirs:
                res = get_polarizability(
                    os.path.join(dirname, sd),
                    calc_dir_patters,
                    deepth - 1,
                    verbos=verbos,
                    threshold=threshold,
                )
                if sd not in all_res["sub_dir"] and len(res) > 0:
                    all_res["sub_dir"][os.path.basename(sd)] = res
    return all_res


def do_one_basis(output_lis, threshold=None):
    if not len(output_lis):
        warnings.warn("there is no valid output file here")
        return

    # check the calc_type from the first output file
    ct = output_lis[0].inp.calc_type
    calc_orbit = output_lis[0].calc_orbit
    pc = PolarizabilityCalculator(calc_type=ct)

    # maybe we just use a dict to restore all information
    energies = {}
    # cause for different type calculations, the electric fields are the same
    fields = []

    for o in output_lis:
        if o.inp.calc_type != ct:
            warnings.warn(
                "we found a error output whose calc_type {} is "
                "different with the first one {}".format(o.inp.calc_type, ct)
            )
            continue

        # CC or CI
        fields.append(float(o.inp.electric_field))
        if o.inp.calc_method == "CC":
            for k, v in o.energy_settings.items():
                if k not in energies.keys():
                    energies[k] = []
                energies[k].append(v)

        elif o.inp.calc_method == "CI":
            # check all roots converged
            for i in o.energy_settings["ci_converged"]:
                if not np.all(np.asarray(i)):
                    raise RuntimeError("Not all roots are converged!")

            for k, v in o.energy_settings["ci_e"].items():
                if k not in energies.keys():
                    energies[k] = []
                energies[k].append(v)

            if "scf_e" not in energies:
                energies["scf_e"] = []
            energies["scf_e"].append(o.energy_settings["scf_e"])
        else:
            # raise NotImplementedError
            continue

    # remove items where the length of energy is not equal to the one of field
    del_k_lis = []
    for k, v in energies.items():
        if len(v) != len(fields):
            warnings.warn(
                'The length of energy set "{}" {} does '
                "not equal the length of field {}".format(k, len(v), len(fields))
            )
            del_k_lis.append(k)
    for k in del_k_lis:
        energies.pop(k)

    # {'scf': alpha_scf, 'mp2':alpha_scf, 'ccsd':alpha_ccsd}
    res = {}
    res_error = {}
    for k, energy in energies.items():
        polar_key = k.strip("_e")
        _polar, _error = pc.get_svd_from_array(energy, fields, threshold=threshold)
        res[polar_key] = _polar
        res_error[polar_key] = _error
    res_all = {
        "energy": {"energies": energies, "fields": fields},
        "polar": res,
        "polar_error": res_error,
    }
    return res_all


def get_polarizability_from_output_list(
    dirname, output_lis, tag=None, verbos=True, threshold=None
):
    all_basis_res = {}
    for o in output_lis:
        task_type, orbit = o.task_type, o.calc_orbit
        obt_info = (
            get_orbital_info(orbit["occ"], orbit["vir"]) if len(orbit) else "null"
        )
        k = get_keyword(o.mol.molecule.atomic_info.symbol, task_type, obt_info)
        if k not in all_basis_res:
            all_basis_res[k] = []
        all_basis_res[k].append(o)

    e_res = {}
    for k, v in all_basis_res.items():
        one_res = do_one_basis(v, threshold=threshold)
        if one_res:
            e_res[k] = one_res
    return e_res


def is_valid(output_lis, verbos=False):
    """Check whether a output_lis is valid

    Args:
        output_lis: a list of Output obj

    Returns:
        True or False
    """

    if len(output_lis) < 3:
        if verbos:
            warnings.warn("the nubmer of output objects is less than 3")
        return False

    # check all there file if they are all 'CC' or 'CI' calculations
    task_record = {}
    for o in output_lis:
        if o.inp.calc_method in ["CC", "CI"]:
            if not o.task_type in task_record:
                task_record[o.task_type] = 1
            else:
                task_record[o.task_type] += 1
    for v in task_record.values():
        if v >= 3:
            return True
    else:
        if verbos:
            print(task_record)
            warnings.warn(
                "the maximum of output objects with the same type is less than 3"
            )
        return False
