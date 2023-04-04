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
        for f, e in filed_energy:
            f = float(f)
            energy.append(float(e))
            field.append(f)

        x_svd = self.get_svd_from_array(energy, field)
        return x_svd
