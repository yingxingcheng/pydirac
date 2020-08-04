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

"""
TODO:
This script is used to collect all results from different type calculations,
e.g., CC or CI with different Hamiltonians and different calculation type, i.e.,
dipole and quadrupole. Using this results, one can generate atomic information
table and write it in a tex file.
"""

import os
import shutil
import glob
import numpy as np
from pathlib import Path
import warnings
from monty.os import cd
from pydirac.io.outputs import Output
from pydirac.analysis.polarizability import PolarizabilityCalculator


def get_polarizabiltiy(dirname: str = './',
                       sub_dir_tag: str ='*dyall*',
                       deepth: int =0) -> None:
    """Get polarizability from a directory

    Args:
        dirname: dirname
        suffix:

    Returns:
        None
    """

    output_lis = []
    with cd(dirname):
        # here, maybe 'JOB_DONE' is not the file which specify the end of
        # a calculation, e.g., basis_set_JOB_DONE
        #
        # but, for q_mrci type calculation, JOB_DONE exists and the output
        # is a SCF calculation normally.
        find_dir = False
        if 'JOB_DONE' in glob.glob('*'):
            for f in glob.glob('*.out'):
                output_lis.append(Output(filename=f))

            # check all there file if they are all 'CC' or 'CI' calculations
            for o in output_lis:
                if not o.calc_method in ['CC' or 'CI']:
                    find_dir = False
                    output_lis = []
                    break
                else:
                    find_dir = True

        find_sub_dir = False
        if not find_dir:
            # check *dyall* type directory which contains individual calculation
            # with a different electric field
            #
            # but, maybe there exist two basis-set type calculations
            for sub_dir in glob.glob('*' + sub_dir_tag + '*'):
                sub_dir_fullname = os.path.join(dirname, sub_dir)
                if os.path.isdir(sub_dir_fullname):
                    with cd(sub_dir):
                        all_out = glob.glob('*.out')
                        if len(all_out) > 1:
                            warnings.warn('there are more than two output file '
                                          'in this directory {0}, and we take '
                                          'the first one {1}'.format(sub_dir,
                                                                     all_out[0]))
                        elif len(all_out) == 1:
                            output_lis.append(Output(filename=all_out[0]))
                        else:
                            warnings.warn('there is no output file in this'
                                          ' directory {0}'.format(sub_dir))

        if len(output_lis):
            find_sub_dir = True

        if not find_sub_dir:
            warnings.warn('current directory does not have any calc directory '
                          'whose name contain "dyall" and useful output file '
                          'please check it')
            return

    def _deal_with_single_basis(output_lis):
        # here, we are do the one basis_set
        if not len(output_lis):
            warnings.warn('there is no valid output file here')
            return

        # for o in output_lis:
        #     print(o.task_type)
        # return

        # check the calc_type from the first output file
        calc_type = output_lis[0].calc_type
        if calc_type not in 'QD':
            warnings.warn('Unknown type calculation!')
            return

        if calc_type == 'Q':
            pc = PolarizabilityCalculator(calc_type='quadrupole')
        else:
            pc = PolarizabilityCalculator(calc_type='dipole', )

        # maybe we just use a dict to restore all information
        energies = {}
        # cause for different type calculations, the electric fields are the same
        fields = []

        for o in output_lis:
            if o.calc_type != calc_type:
                warnings.warn('we found a error output whose calc_type {0} is '
                              'different with the first one {1}'.format(
                    o.calc_type, calc_type))
                continue

            # TODO: different calc formula for CC and CI
            # check CC or CI
            if o.calc_method not in ['CC', 'CI']:
                warnings.warn('This calculation {0} does not belong any '
                              'method "CC" or "CI".'.format(o.calc_method))
                continue

            # CC or CI
            fields.append(float(o.electric_field))
            if o.calc_method == 'CC':
                for k, v in  o.energy_settings.items():
                    if k not in energies.keys():
                        energies[k] = []
                    energies[k].append(v)

            else:
                # check all roots converged
                if not np.alltrue(np.asarray(o.energy_settings['ci_converged'])):
                    warnings.warn('not all roots are converged, please use data carfully')

                for k, v in o.energy_settings['ci_e'].items():
                    if k not in energies.keys():
                        energies[k] = []
                    energies[k].append(v)

        res = {}
        for k, v in energies.items():
            res[k] = pc.get_svd_from_array(v, fields)
        return res

    # extract all basis_type info
    all_basis_res = {}
    for o in output_lis:
        k = o.task_type
        if not k in all_basis_res.keys():
            all_basis_res[k] = []
        all_basis_res[k].append(o)

    e_res = {}
    for k, v in all_basis_res.items():
        e_res[k] = _deal_with_single_basis(v)

    for k, v in e_res.items():
        print('Table: results of {0} \n\tfrom {1}'.format(k, dirname))
        print('='*80)
        print('  {0:<20s} {1:<20s} {2:<20s}'.format('calc_method', 'momentum', 'polarizability'))
        print('-'*80)
        for i_k, i_v in v.items():
            print('  {0:<20s} {1:<20.3f} {2:<20.3f}'.format(i_k, i_v[0], i_v[1]))
        print('='*80)
        print()


    if deepth > 0:
        calc_from_sub_dir(dirname,sub_dir_tag, deepth-1)


def calc_from_sub_dir(dirname, exclude, deepth):
    sub_dirs = [ sub_path for sub_path in glob.glob('{0}/*'.format(dirname))
                 if os.path.isdir(sub_path) and exclude not in sub_path]

    for sub_dir in sub_dirs:
        sub_dir_fullname = os.path.join(dirname, sub_dir)
        get_polarizabiltiy(sub_dir_fullname, sub_dir_tag=exclude, deepth=deepth)



def get_atomDB(args):
    print(args)
    dir_fullname = os.path.abspath(args.dirname)
    get_polarizabiltiy(dir_fullname, sub_dir_tag=args.sub_dir_tag,
                       deepth=args.deepth )


