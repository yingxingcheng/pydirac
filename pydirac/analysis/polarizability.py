#!/usr/bin/env python
import os
import numpy as np


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

        if nb_coeffs not in (2, 3):
            raise ValueError("nb_coeffs must be 2 or 3!")
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
        if self.nb_coeffs == 2:
            if self.calc_type == "dipole":
                return np.array([-np.power(f, 2) / 2.0, -np.power(f, 4) / 24.0])
            elif self.calc_type == "quadrupole":
                return np.array([np.power(f, 1) / 2.0, -np.power(f, 2) / 8.0])
            else:
                raise TypeError(
                    'calc_type {0} is not valid, please use "dipole" '
                    'or "quadrupole"'.format(self.calc_type)
                )

        elif self.nb_coeffs == 3:
            if self.calc_type == "dipole":
                return np.array(
                    [-np.power(f, 2) / 2.0, -np.power(f, 4) / 24.0, -np.power(f, 3) / 144.0]
                )
            elif self.calc_type == "quadrupole":
                return np.array(
                    [-np.power(f, 1) / 2.0, -np.power(f, 2) / 8.0, -np.power(f, 3) / 24.0]
                )
            else:
                raise TypeError(
                    'calc_type {0} is not valid, please use "dipole" '
                    'or "quadrupole"'.format(self.calc_type)
                )
        else:
            raise ValueError(
                f"calc_type {self.calc_type} is not valid, please use 'dipole' or 'quadrupole'."
            )

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

    def get_svd_from_array(self, energy, field):
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
        e_ref = 0.0
        for e, f in zip(energy, field):
            if np.isclose(f, 0.0):
                e_ref = e

        # define a matrix
        A = self.get_A(np.asarray(field))
        b = np.asarray(energy) - e_ref

        x_svd = np.linalg.lstsq(A, b, rcond=1e-8)[0]
        return x_svd

    def get_res_svd(self, filename):
        """
        Calculate the SVD-based solution for the least-squares problem using the data from a file.

        Parameters
        ----------
        filename : str
            The name of the file containing the electric field strengths and energy values.

        Returns
        -------
        array_like of float
            The coefficients.
        """
        with open(filename, "r") as f:
            context = f.readlines()

        if len(context) < 3:
            raise ValueError(f"Data points in {filename} is less than 3")

        filed_energy = [tuple(c.split()) for c in context]
        energy, field = [], []
        e_ref = 0.0
        for f, e in filed_energy:
            f = np.float(f)
            if np.isclose(f, 0.00):
                e_ref = np.float(e)
            energy.append(np.float(e))
            field.append(f)

        x_svd = self.get_svd_from_array(energy, field)
        return x_svd


def get_dipole_polarizability_from_cc(opt, option, file_list, parser):
    pc = PolarizabilityCalculator(calc_type="dipole", nb_coeffs=2)

    print(" Begin ".center(80, "*"))
    methods, res = [], []
    for f in file_list:
        filename = os.path.basename(f)
        calc_type = filename.split("_")[1]
        if "(" in calc_type:
            calc_type = "CCSD_T"
        x_svd = pc.get_res_svd(f)
        methods.append("Res of {0}".format(f))
        res.append((x_svd[0]))

        with open("polar_" + calc_type, "w") as f:
            f.write(" ".join([str(i) for i in x_svd[0:2]]))

    print(" ".join(methods))
    print(" ".join(["{0:.3f}".format(r) for r in res]))


def get_quadrupole_polarizability_from_cc(opt, opotion, file_list, parser):
    pc = PolarizabilityCalculator(calc_type="quadrupole", nb_coeffs=2)

    print(" Begin ".center(80, "*"))
    for f in file_list:
        x_svd = pc.get_res_svd(f)
        x_svd = np.asarray(x_svd)
        x_svd[1] = x_svd[1] / 4.0
        print(
            "{0}: \n\tquadru momentum: {1:.3f} \n\tquadru polarizability: "
            "{2:.3f} (a.u.)".format(f, x_svd[0], x_svd[1])
        )


def get_quadrupole_polarizability_from_ci(opt, option, file_list, parser):
    pc = PolarizabilityCalculator(calc_type="quadrupole", nb_coeffs=2)

    methods, res = [], []
    total_momentum = []
    total_polarizability = []
    print(" Begin ".center(80, "*"))
    latex_res = []
    sym_idx = 1
    root_idx = "Root 1"
    for f in file_list:
        calc_type = f.split(".")[0]
        x_svd = pc.get_res_svd(f)
        if x_svd is None:
            continue

        x_svd[1] = x_svd[1] / 4.0
        print(
            "{0}: \n\tquadru momentum: {1:.3f} \n\tquadru "
            "polarizability: {2:.3f} (a.u.)".format(calc_type, x_svd[0], x_svd[1])
        )

        words = calc_type.split("_")
        if len(words) > 2:
            sym_idx = int(words[-2])
            root_idx = "Root " + words[-1]
        else:
            if "scf" in words[1]:
                root_idx = "DHF"
            else:
                raise IndexError("Invalid index: {0}".format(words[1]))

        if sym_idx == 1:
            latex_str = "{0} & {1:.3f} & {2:.3f} & & \\\\ \\hline".format(
                root_idx, x_svd[0], x_svd[1]
            )
        else:
            latex_str = "{0} & & & {1:.3f} & {2:.3f} \\\\ \\hline".format(
                root_idx, x_svd[0], x_svd[1]
            )
        latex_res.append(latex_str)

        if "scf" not in calc_type:
            total_momentum.append(x_svd[0])
            total_polarizability.append(x_svd[1])

        with open("quadru_" + calc_type, "w") as f:
            f.write(" ".join([str(i) for i in x_svd[0:2]]))

    print(" ")
    print("The final resutls: ")
    for m, r in zip(methods, [str(r) for r in res]):
        print("{0} : {1}".format(m, r))

    if pc.calc_type == "quadrupole" and len(total_momentum) and len(total_polarizability):
        ave_m = sum(total_momentum) / len(total_momentum)
        ave_p = sum(total_polarizability) / len(total_polarizability)
        print("Average total quadrupole momentum is: {0:.3f}".format(ave_m))
        print("Average total quadrupole polarizability is: {0:.3f}".format(ave_p))
        print(" ")
        print(
            "Note: the different roots are not all degenerate states, \n"
            "and ground states and excited states are both included here.\n "
            "For lighter atoms, only this method is valid to reproduce \n"
            "correct quadrupole polarizability. This is because an \n"
            "average-over-state SCF used in DIRAC, and for lighter atoms, \n"
            "these states are all degenerated for non-relativistic \n"
            "calculations but not in spin-orbit coupling calculations"
        )
        if sym_idx == 1:
            latex_str = "Ave. & {0:.3f} & {1:.3f} & & \\\\ \\hline".format(ave_m, ave_p)
        else:
            latex_str = "Ave. & & & {0:.3f} & {1:.3f} \\\\ \\hline".format(ave_m, ave_p)
        latex_res.append(latex_str)
        print(" End ".center(80, "*"))
        print(" ")
    else:
        print(" End ".center(80, "*"))

    with open("ref_{0}E.tex".format(sym_idx), "w") as f:
        f.write("\n".join(latex_res))
