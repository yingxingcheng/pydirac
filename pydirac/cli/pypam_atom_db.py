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
import glob
import numpy as np
import warnings
from monty.os import cd
from pydirac.io.outputs import Output
from pydirac.analysis.polarizability import PolarizabilityCalculator


def get_polarizability(dirname: str = './',
                       sub_dir_tag: str ='*dyall*',
                       deepth: int =0) -> None:
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
    if deepth < 0: deepth = 0
    do_curr_dir = False
    do_clc_dir = False
    do_sub_dir = deepth > 0

    # for debug
    print('cd {0}'.format(dirname))

    curr_dir_output_lis = []
    calc_dir_output_lis = []
    with cd(dirname):
        # Step 1. deal with all current output files
        # ------------------------------
        for f in glob.glob('*.out'):
            obj = Output(filename=f)
            if obj.is_ok:
                curr_dir_output_lis.append(Output(filename=f))

        if is_valid(curr_dir_output_lis):
            do_curr_dir = True

        # Step 2. deal with all calc_dir
        # ------------------------------
        clc_dirs = [d for d in glob.glob('*' + sub_dir_tag + '*') if
                      os.path.isdir(d)]
        # if has calc_dir then I would try
        for clc_d in clc_dirs:
            with cd(clc_d):
                outs = glob.glob('*.out')
                if len(outs) > 1:
                    warnings.warn('there are more than two output file '
                                  'in this directory {0}, and we take '
                                  'the first one {1}'.format(clc_d,
                                                             outs[0]))
                elif len(outs) == 1:
                    obj = Output(filename=outs[0])
                    if obj.is_ok:
                        calc_dir_output_lis.append(obj)
                else:
                    warnings.warn('there is no output file in this'
                                  ' directory {0}'.format(clc_d))

        # check if all output obj in calc_dir_output_lis is valid and
        # do we need to go on
        if is_valid(calc_dir_output_lis):
            do_clc_dir = True

        if do_curr_dir:
            get_polarizability_from_output_list(dirname, curr_dir_output_lis, tag='curr_dir')
        if do_clc_dir:
            get_polarizability_from_output_list(dirname, calc_dir_output_lis, tag='clc_dir')

        # Step 3. deal with all sub_dir
        if do_sub_dir:
            sub_dirs = [d for d in glob.glob('*') if os.path.isdir(d)
                        and d not in clc_dirs]

            for sd in sub_dirs:
                get_polarizability(os.path.join(dirname, sd), sub_dir_tag, deepth-1)


def get_polarizability_from_output_list(dirname, output_lis, tag=None):
    def _deal_with_single_basis(output_lis):
        # here, we are do the one basis_set
        if not len(output_lis):
            warnings.warn('there is no valid output file here')
            return

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

        del_k_lis = []
        for k, v in energies.items():
            if len(v) != len(fields):
                warnings.warn('The length of energy set "{0}" {1} does '
                              'not equal the length of field {2}'.format(k, len(v), len(fields)))
                del_k_lis.append(k)
        for k in del_k_lis:
            energies.pop(k)

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
        single_res = _deal_with_single_basis(v)
        if single_res:
            e_res[k] = single_res

    for k, v in e_res.items():
        if tag:
            print('Table: results of {0} from \n "{1}" ({2})'.format(k, dirname, str(tag)))
        else:
            print('Table: results of {0} from \n "{1}"'.format(k, dirname))
        print('='*80)
        if k.startswith('D'):
            print(
                '  {0:<20s} {1:<20s} {2:<20s}'.format('method', 'polarizability',
                                                      'hyper-polarizability'))
        else:
            print(
                '  {0:<20s} {1:<20s} {2:<20s}'.format('method', 'momentum',
                                                      'polarizability'))

        print('-'*80)
        for i_k, i_v in v.items():
            print('  {0:<20s} {1:<20.3f} {2:<20.3f}'.format(i_k, i_v[0], i_v[1]))
        print('='*80)
        print()


def is_valid(output_lis):
    """Check whether a output_lis is valid

    Args:
        output_lis: a list of Output obj

    Returns:
        True or False
    """

    if len(output_lis) < 3:
        return False

    # check all there file if they are all 'CC' or 'CI' calculations
    task_record = {}
    for o in output_lis:
        if not o.calc_method in ['CC', 'CI']:
            continue
        else:
            if not o.task_type in task_record:
                task_record[o.task_type] = 0
            else:
                task_record[o.task_type] += 1
    for v in task_record.values():
        if v >= 3:
            return True
    else:
        return False


def get_atomDB(args):
    dir_fullname = os.path.abspath(args.dirname)
    get_polarizability(dir_fullname, sub_dir_tag=args.sub_dir_tag,
                       deepth=args.deepth )


