#!/usr/bin/env python

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

import os
import numpy as np
from scipy.linalg import svd
np.set_printoptions(suppress=True)

#DIPOLE = False
DIPOLE = True

class PolarizabilityCalculator():
    """
    A class to calculate dipole polarizability and quadrupole polarizability
    """

    def __init__(self, calc_type='dipole', nb_coeffs=2):
        self.calc_type = calc_type
        self.nb_coeffs = nb_coeffs

    def get_coeff(self, f):
        if self.nb_coeffs == 2:
            if self.calc_type == 'dipole':
                return np.array([-np.power(f, 2)/2., -np.power(f, 4)/24.])
            elif self.calc_type == 'quadrupole':
                return np.array([np.power(f, 1)/2., -np.power(f, 2)/8.])
            else:
                raise TypeError(
                    'calc_type {0} is not valid, please use "dipole" '
                    'or "quadrupole"'.format(self.calc_type))

        elif self.nb_coeffs == 3:
            if self.calc_type == 'dipole':
                return np.array([-np.power(f, 2)/2., -np.power(f, 4)/24.
                                 -np.power(f, 3)/144.])
            elif self.calc_type == 'quadrupole':
                return np.array([-np.power(f, 1)/2., -np.power(f, 2)/8.,
                                 -np.power(f, 3)/24.])
            else:
                raise TypeError(
                    'calc_type {0} is not valid, please use "dipole" '
                    'or "quadrupole"'.format(self.calc_type))
        else:
            raise TypeError('np_coeffs {0} is not valid, please use 2 '
                            'or 3!'.format(self.nb_coeffs))

    def get_A(self, f):
        A = np.zeros((len(f), self.nb_coeffs))
        for i, v in enumerate(f):
            A[i, :] = self.get_coeff(v)
        return A

    def get_res_svd(self, filename):
        with open(filename, 'r') as f:
            context = f.readlines()

        if len(context) < 3:
            raise ValueError('Data points in {0} is less than 3'.format(filename))

        filed_energy = [tuple(c.split()) for c in context]
        energy, field = [], []
        e_ref = 0.0
        for f, e in filed_energy:
            f = np.float(f)
            if np.isclose(f, 0.00):
                e_ref = np.float(e)
            energy.append(np.float(e))
            field.append(f)

        # define a matrix
        A = self.get_A(np.asarray(field))
        # SVD
        U, sigma, VT = svd(A)
        # Make a matrix Sigma of the correct size:
        Sigma = np.zeros(A.shape)
        Sigma[:A.shape[1], :A.shape[1]] = np.diag(sigma)

        Sigma_prinv = np.zeros(A.shape).T
        Sigma_prinv[:self.nb_coeffs, :self.nb_coeffs] = \
            np.diag(1/sigma[:self.nb_coeffs])

        # Now compute the SVD-based solution for the least-squares problem
        b = np.asarray(energy) - e_ref
        x_svd = VT.T.dot(Sigma_prinv).dot(U.T).dot(b)
        return x_svd


def get_dipole_polarizability_from_cc(opt, option, file_list, parser):
    pc = PolarizabilityCalculator(calc_type='dipole', nb_coeffs=2)

    print(" Begin ".center(80, '*'))
    methods, res = [], []
    for f in file_list:
        filename = os.path.basename(f)
        calc_type = filename.split('_')[1]
        if '(' in calc_type:
            calc_type = 'CCSD_T'
        x_svd = pc.get_res_svd(f)
        methods.append('Res of {0}'.format(f))
        res.append((x_svd[0]))

        with open('polar_' + calc_type, 'w') as f:
            f.write(' '.join([str(i) for i in x_svd[0:2]]))

    print(' '.join(methods))
    print(' '.join(['{0:.3f}'.format(r) for r in res]))


def get_quadrupole_polarizability_from_cc(opt, opotion, file_list, parser):
    pc = PolarizabilityCalculator(calc_type='quadrupole', nb_coeffs=2)

    print(" Begin ".center(80, '*'))
    for f in file_list:
        x_svd = pc.get_res_svd(f)
        x_svd = np.asarray(x_svd)
        x_svd[1] = x_svd[1]/4.
        print('{0}: \n\tquadru momentum: {1:.3f} \n\tquadru polarizability: '
              '{2:.3f} (a.u.)'.format(f, x_svd[0], x_svd[1]))


def get_quadrupole_polarizability_from_ci(opt, option, file_list, parser):
    pc = PolarizabilityCalculator(calc_type='quadrupole', nb_coeffs=2)

    methods, res = [], []
    total_momentum = []
    total_polarizability = []
    print(' Begin '.center(80, '*'))
    latex_res = []
    sym_idx = 1
    root_idx = 'Root 1'
    for f in file_list:
        calc_type = f.split('.')[0]
        x_svd = pc.get_res_svd(f)
        if x_svd is None:
            continue

        x_svd[1] = x_svd[1]/4.
        print('{0}: \n\tquadru momentum: {1:.3f} \n\tquadru '
              'polarizability: {2:.3f} (a.u.)'.
              format(calc_type, x_svd[0], x_svd[1]))

        words = calc_type.split('_')
        if len(words) > 2:
            sym_idx = int(words[-2])
            root_idx = 'Root ' + words[-1]
        else:
            if 'scf' in words[1]:
                root_idx = 'DHF'
            else:
                raise IndexError('Invalid index: {0}'.format(words[1]))

        if sym_idx == 1:
            latex_str = '{0} & {1:.3f} & {2:.3f} & & \\\\ \\hline'.\
                format(root_idx, x_svd[0], x_svd[1])
        else:
            latex_str = '{0} & & & {1:.3f} & {2:.3f} \\\\ \\hline'.\
                format(root_idx, x_svd[0], x_svd[1])
        latex_res.append(latex_str)

        if 'scf' not in calc_type:
            total_momentum.append(x_svd[0])
            total_polarizability.append(x_svd[1])

        with open('quadru_' + calc_type, 'w') as f:
            f.write(' '.join([str(i) for i in x_svd[0:2]]))

    print(' ')
    print('The final resutls: ')
    for m, r in zip(methods, [str(r) for r in res]):
        print('{0} : {1}'.format(m, r))

    if pc.calc_type == 'quadrupole' and len(total_momentum) \
            and len(total_polarizability):
        ave_m = sum(total_momentum)/ len(total_momentum)
        ave_p = sum(total_polarizability) / \
                             len(total_polarizability)
        print("Average total quadrupole momentum is: {0:.3f}".format(ave_m))
        print("Average total quadrupole polarizability is: {0:.3f}".format(ave_p))
        print(" ")
        print("Note: the different roots are not all degenerate states, \n"
              "and ground states and excited states are both included here.\n "
              "For lighter atoms, only this method is valid to reproduce \n"
              "correct quadrupole polarizability. This is because an \n"
              "average-over-state SCF used in DIRAC, and for lighter atoms, \n"
              "these states are all degenerated for non-relativistic \n"
              "calculations but not in spin-orbit coupling calculations")
        if sym_idx == 1:
            latex_str = 'Ave. & {0:.3f} & {1:.3f} & & \\\\ \\hline'.format(ave_m, ave_p)
        else:
            latex_str = 'Ave. & & & {0:.3f} & {1:.3f} \\\\ \\hline'.format(ave_m, ave_p)
        latex_res.append(latex_str)
        print(" End ".center(80, '*'))
        print(" ")
    else:
        print(" End ".center(80, '*'))

    with open('ref_{0}E.tex'.format(sym_idx), 'w') as f:
        f.write('\n'.join(latex_res))


if __name__ == '__main__':
    pass
