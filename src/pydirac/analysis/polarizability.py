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
    "PolarizabilityCalculator",
    "get_polarizability",
    "get_polarizability_from_output_list",
    "do_one_basis",
    "is_valid",
]


def calc_mat(fields, nb_coeffs, calc_type):
    assert calc_type in ["dipole", "quadrupole"]
    if calc_type == "dipole":
        f_list = []
        for i in range(nb_coeffs):
            order = 2 * i + 2
            f_list.append(-1 / factorial(order) * fields**order)
    elif calc_type == "quadrupole":
        assert 1 <= nb_coeffs <= 3
        if nb_coeffs == 1:
            f_list = [1 / 2 * fields]
        elif nb_coeffs == 2:
            f_list = [1 / 2 * fields, -1 / 8 * fields**2]
        elif nb_coeffs == 3:
            f_list = [1 / 2 * fields, -1 / 8 * fields**2, 1 / 24 * fields**3]
        else:
            raise RuntimeError(
                f"nb_coeffs={nb_coeffs} should be less than 4 and larger than 0"
            )
    else:
        raise NotImplementedError()

    A = np.column_stack(f_list)
    return A


def check_results(x_svd, errors, nb_coeffs, atol=0.1):
    if nb_coeffs > 1:
        if np.isclose(x_svd[1], 0.0):
            assert np.allclose(x_svd[1:], 0.0) and np.allclose(errors[1:], 0.0)
            return True

        # The hyperpolarizability of atom should be positive.
        if x_svd[1] < 0:
            return False

        # The errors are too large
        if errors[1] / x_svd[1] < atol:
            return False
    return True


def fit_ff_eq(
    energy,
    field,
    calc_mat,
    nb_coeffs,
    rcond,
    threshold,
    calc_type,
    check_hyper_polar,
):
    """
    Calculate the SVD-based solution for the least-squares problem.

    Parameters
    ----------
    energy : array_like of float
        The energy values.
    field : array_like of float
        The electric field strengths.

    Returns
    -------
    array_like of float
        The coefficients
    """
    assert calc_type in ["dipole", "quadrupole"]
    field = np.asarray(field)
    energy = np.asarray(energy)

    # Reference energy for zero field
    zero_field_mask = np.isclose(field, 0.0)
    if np.any(zero_field_mask):
        e_ref = energy[zero_field_mask][0]
    else:
        raise ValueError("No field value is close to zero.")

    # Create a mask to filter out values where the field is 0 or larger than the threshold
    mask = (field != 0.0) & (field <= threshold + 1e-8)

    # Apply the mask to both field and energy arrays
    reduced_field = field[mask]
    reduced_energy = energy[mask]

    # Sort the field array and use the same indices to sort the energy array
    sort_indices = np.argsort(reduced_field)
    sorted_field = reduced_field[sort_indices]
    sorted_energy = reduced_energy[sort_indices]

    # Define the matrix
    A = calc_mat(sorted_field, nb_coeffs, calc_type)
    b = sorted_energy - e_ref

    # Normalize the columns of A
    norms = np.linalg.norm(A, axis=0)
    A_normalized = A / norms

    # Solve the least squares problem using SVD on the normalized matrix
    x_svd_normalized, residuals, rank, s = lstsq(A_normalized, b, cond=rcond)

    # Adjust the solution to account for the normalization
    x_svd = x_svd_normalized / norms

    # Calculate the residual sum of squares
    residual_sum_of_squares = np.sum(residuals**2)

    # Check for degrees of freedom
    degrees_of_freedom = len(reduced_energy) - len(x_svd)
    if degrees_of_freedom == 0:
        # This implies an exact fit; no error to estimate.
        MSE = 0
        # Covariance matrix is not defined
        cov_matrix = np.zeros((len(x_svd), len(x_svd)))
    else:
        # Calculate the covariance matrix of the coefficients
        # Assuming homoscedasticity (constant variance of residuals)
        MSE = residual_sum_of_squares / degrees_of_freedom  # Mean Squared Error
        cov_matrix = np.linalg.inv(np.dot(A.T, A)) * MSE

    # Standard errors are the square roots of the diagonal elements of the covariance matrix
    standard_errors = np.sqrt(np.diag(cov_matrix))

    if check_hyper_polar:
        if check_results(x_svd, standard_errors, nb_coeffs):
            return x_svd, standard_errors
        elif nb_coeffs > 1:
            x_svd, standard_errors = fit_ff_eq(
                energy, field, calc_mat, 1, rcond, threshold, calc_type
            )
            tmp_svd = np.zeros((nb_coeffs,))
            tmp_errors = np.zeros((nb_coeffs,))
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

    def get_coeff(self, f):
        """
        Return the coefficients for a given electric field strength.

        Parameters
        ----------
        f : float
            The electric field strength.

        Returns
        -------
        array_like of float
            The coefficients.
        """
        coeffs = {
            "dipole": {
                1: [-np.power(f, 2) / 2.0],
                2: [-np.power(f, 2) / 2.0, -np.power(f, 4) / 24.0],
                3: [
                    -np.power(f, 2) / 2.0,
                    -np.power(f, 4) / 24.0,
                    -np.power(f, 6) / 720.0,
                ],
            },
            "quadrupole": {
                1: [np.power(f, 1) / 2.0],
                2: [np.power(f, 1) / 2.0, -np.power(f, 2) / 8.0],
                3: [np.power(f, 1) / 2.0, -np.power(f, 2) / 8.0, np.power(f, 3) / 24.0],
            },
        }

        if self.calc_type not in coeffs:
            raise TypeError(
                f'calc_type {self.calc_type} is not valid, please use "dipole" or "quadrupole".'
            )

        if self.nb_coeffs not in coeffs[self.calc_type]:
            raise ValueError(
                f"nb_coeffs {self.nb_coeffs} is not valid for calc_type {self.calc_type}."
            )

        return np.array(coeffs[self.calc_type][self.nb_coeffs])

    def get_A(self, field):
        """
        Create the matrix A for the given electric field strengths.

        Parameters
        ----------
        field : array_like of float
            The electric field strengths.

        Returns
        -------
        array_like of float, shape (n_samples, nb_coeffs)
            The matrix A for the given electric field strengths.
        """
        A = np.zeros((len(field), self.nb_coeffs))
        for i, v in enumerate(field):
            A[i, :] = self.get_coeff(v)
        return A

    def get_svd_from_array_old(self, energy, field, rcond=1e-8, threshold=None):
        """
        Calculate the SVD-based solution for the least-squares problem.

        Parameters
        ----------
        energy : array_like of float
            The energy values.
        field : array_like of float
            The electric field strengths.

        Returns
        -------
        array_like of float
            The coefficients
        """
        field = np.asarray(field)
        energy = np.asarray(energy)

        # Reference energy for zero field
        zero_field_mask = np.isclose(field, 0.0)
        if np.any(zero_field_mask):
            e_ref = energy[zero_field_mask][0]
        else:
            raise ValueError("No field value is close to zero.")

        # Create a mask to filter out values where the field is 0 or larger than the threshold
        mask = (field != 0.0) & (field <= threshold + 1e-8)

        # Apply the mask to both field and energy arrays
        field = field[mask]
        energy = energy[mask]

        # Sort the field array and use the same indices to sort the energy array
        sort_indices = np.argsort(field)
        field = field[sort_indices]
        energy = energy[sort_indices]

        # Define the matrix
        A = self.get_A(field)
        b = energy - e_ref

        # Normalize the columns of A
        norms = np.linalg.norm(A, axis=0)
        A_normalized = A / norms

        # Solve the least squares problem using SVD on the normalized matrix
        x_svd_normalized, residuals, rank, s = lstsq(A_normalized, b, cond=rcond)

        # Adjust the solution to account for the normalization
        x_svd = x_svd_normalized / norms

        # Check condition to fit only the first coefficient
        if self.nb_coeffs > 1 and (x_svd[1] < 0 or x_svd[1] > 1e7):
            # Fit only the first coefficient
            A_normalized = A_normalized[:, [0]]  # Reduce A_normalized to single column
            x_svd_normalized, residuals, rank, s = lstsq(A_normalized, b, cond=rcond)
            x_svd = np.zeros(self.nb_coeffs)
            x_svd[0] = x_svd_normalized[0] / norms[0]  # Adjust the first coefficient

            # Update the residuals for the single coefficient fit
            residuals = b - np.dot(A[:, [0]], x_svd[:1])
        else:
            # Manually calculate residuals if not provided by scipy.linalg.lstsq
            if np.ndim(residuals) == 0:  # Scalar residual case
                residuals = b - np.dot(A, x_svd)

        # Calculate the residual sum of squares
        residual_sum_of_squares = np.sum(residuals**2)

        # Check for degrees of freedom
        degrees_of_freedom = len(energy) - len(x_svd)
        if degrees_of_freedom == 0:
            # This implies an exact fit; no error to estimate.
            MSE = 0
            # Covariance matrix is not defined
            cov_matrix = np.zeros((len(x_svd), len(x_svd)))
        else:
            # Calculate the covariance matrix of the coefficients
            # Assuming homoscedasticity (constant variance of residuals)
            MSE = residual_sum_of_squares / degrees_of_freedom  # Mean Squared Error
            cov_matrix = np.linalg.inv(np.dot(A.T, A)) * MSE

        # Standard errors are the square roots of the diagonal elements of the covariance matrix
        standard_errors = np.sqrt(np.diag(cov_matrix))
        return x_svd, standard_errors

    def get_svd_from_array(
        self, energy, field, rcond=None, threshold=None, check_hyper_polar=False
    ):
        return fit_ff_eq(
            energy=energy,
            field=field,
            calc_mat=calc_mat,
            nb_coeffs=self.nb_coeffs,
            rcond=rcond,
            threshold=threshold,
            calc_type=self.calc_type,
            check_hyper_polar=check_hyper_polar,
        )


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
