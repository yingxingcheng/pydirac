#!/usr/bin/env python

import numpy as np
from scipy.linalg import svd
np.set_printoptions(suppress=True)

#DIPOLE = False
DIPOLE = True

def get_res_svd(file_name, NB_COEFFS=3):

    def get_coeff(f):
        if NB_COEFFS == 2:
            if DIPOLE:
                return np.array([-np.power(f, 2)/2, -np.power(f, 4)/24])
            else:
                # qudratuple polarizability
                return np.array([-n.power(f, 1)/2, -np.power(f, 2)/8])
        else:
            if DIPOLE:
                return np.array([-np.power(f, 2)/2, -np.power(f, 4)/24,
                                 -np.power(f, 6)/144])
            else:
                return np.array([-np.power(f, 1)/2, -np.power(f, 2)/8,
                                 -np.power(f, 3)/24])

    def get_A(f):
        # f = np.linspace(+0.019, 0.00, 20)
        # print(f)
        A = np.zeros((len(f), NB_COEFFS))
        for i, v in enumerate(f):
            A[i, :] = get_coeff(v)
        return A


    with open(file_name, 'r') as f:
        context = f.readlines()

    field_energy = [tuple(c.split()) for c in context]
    if len(field_energy) < 3:
        print('Wrong: Data point is not enough!')
    elif len(field_energy) == 3:
        NB_COEFFS = 2
    else:
        NB_COEFFS = 3

    energy, field = [], []
    e_ref = 0.0
    for f, e in field_energy:
        f =  np.float(f)
        if np.isclose(f, 0.00):
            e_ref = np.float(e)
        energy.append(np.float(e))
        field.append(f)

    # define a matrix
    A = get_A(np.asarray(field))

    # SVD
    U, sigma, VT = svd(A)

    # Make a matrix Sigma of the correct size:
    Sigma = np.zeros(A.shape)
    Sigma[:A.shape[1], :A.shape[1]] = np.diag(sigma)

    # check that we've attually factorized A:
    #print((U.dot(Sigma).dot(VT) - A).round(4))

    Sigma_pinv = np.zeros(A.shape).T
    Sigma_pinv[:NB_COEFFS, :NB_COEFFS] = np.diag(1/sigma[:NB_COEFFS])
    # print(Sigma_pinv.round(3))

    # Now compute the SVD-based solution for the least-squares problem
    b = np.asarray(energy) - e_ref
    x_svd = VT.T.dot(Sigma_pinv).dot(U.T).dot(b)
    return x_svd


if __name__ == '__main__':
    import sys

    argv = sys.argv[1:]
    methods, res = [], []
    for f in argv:
        x_svd = get_res_svd(f)
        methods.append('Res of {0}: '.format(f))
        res.append((x_svd[0]))

    print(' '.join(methods))
    print(' '.join([str(r) for r in res]))
