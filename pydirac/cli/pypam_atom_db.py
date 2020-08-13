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
                       calc_dir_patters=None,
                       deepth: int =0, verbos=True) -> dict:
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
        calc_dir_patters = ['*dyall*']

    if isinstance(calc_dir_patters, (str, list)):
        if isinstance(calc_dir_patters, str):
            calc_dir_patters = [calc_dir_patters]
    else:
        raise TypeError('calc_dir_patters can only be <str> or <list>')

    if deepth < 0: deepth = 0
    do_curr_dir = False
    do_clc_dir = False
    do_sub_dir = deepth > 0

    # for debug
    if verbos:
        print('cd {0}'.format(dirname))

    all_res = {}
    curr_dir_output_lis = []
    calc_dir_output_lis = []
    with cd(dirname):
        # Step 1. deal with all current output files
        # ------------------------------
        for f in glob.glob('*.out'):
            obj = Output(filename=f)
            if obj.is_ok:
                curr_dir_output_lis.append(obj)

        if is_valid(curr_dir_output_lis):
            do_curr_dir = True

        # Step 2. deal with all calc_dir
        # ------------------------------
        clc_dirs = []
        for ptn in calc_dir_patters:
            clc_dirs.extend([d for d in glob.glob('*' + ptn + '*') if
                          os.path.isdir(d)])
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

        all_res['curr_dir'] = {}
        if do_curr_dir:
            e_dict = get_polarizability_from_output_list(dirname, curr_dir_output_lis, tag='curr_dir')
            all_res['curr_dir'] = e_dict

        all_res['calc_dir'] = {}
        if do_clc_dir:
            e_dict = get_polarizability_from_output_list(dirname, calc_dir_output_lis, tag='clc_dir')
            all_res['calc_dir'] = e_dict

        # Step 3. deal with all sub_dir
        all_res['sub_dir'] = {}
        if do_sub_dir:
            sub_dirs = [d for d in glob.glob('*') if os.path.isdir(d)
                        and d not in clc_dirs]

            for sd in sub_dirs:
                res = get_polarizability(os.path.join(dirname, sd), calc_dir_patters, deepth-1)
                if sd not in all_res['sub_dir'] and len(res) > 0:
                    all_res['sub_dir'][os.path.basename(sd)] = res
    return all_res



def get_polarizability_from_output_list(dirname, output_lis, tag=None):
    def _deal_with_single_basis(output_lis):
        # here, we are do the one basis_set
        if not len(output_lis):
            warnings.warn('there is no valid output file here')
            return

        # check the calc_type from the first output file
        calc_type = output_lis[0].inp.calc_type
        calc_orbit = output_lis[0].calc_orbit
        precision = 'null'
        if len(calc_orbit):
            precision = '(core ' + str(calc_orbit['occ']) + ')[vir ' + str(calc_orbit['vir']) + ']'
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
            if o.inp.calc_type != calc_type:
                warnings.warn('we found a error output whose calc_type {0} is '
                              'different with the first one {1}'.format(
                    o.inp.calc_type, calc_type))
                continue

            # TODO: different calc formula for CC and CI
            # check CC or CI
            if o.inp.calc_method not in ['CC', 'CI']:
                warnings.warn('This calculation {0} does not belong any '
                              'method "CC" or "CI".'.format(o.inp.calc_method))
                continue

            # CC or CI
            fields.append(float(o.inp.electric_field))
            if o.inp.calc_method == 'CC':
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
        return res, precision

    # extract all basis_type info
    all_basis_res = {}
    for o in output_lis:
        k = o.task_type
        if not k in all_basis_res.keys():
            all_basis_res[k] = []
        all_basis_res[k].append(o)

    e_res = {}
    for k, v in all_basis_res.items():
        single_res, precision = _deal_with_single_basis(v)
        if single_res:
            key = k+'@' + precision
            e_res[key] = single_res

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
            # In DIRAC, the correct quadrupole polarizability should be divided 4
            for _k in v.keys():
                v[_k] = v[_k]/4.

        print('-'*80)
        for i_k, i_v in v.items():
            print('  {0:<20s} {1:<20.3f} {2:<20.3f}'.format(i_k.strip('_e'), i_v[0], i_v[1]))
        print('='*80)
        print()
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
            warnings.warn('the nubmer of output objects is less than 3')
        return False

    # check all there file if they are all 'CC' or 'CI' calculations
    task_record = {}
    for o in output_lis:
        if not o.inp.calc_method in ['CC', 'CI']:
            continue
        else:
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
            warnings.warn('the maximum of output objects with the same type '
                          'is less than 3')
        return False


def get_atomDB(args):
    if isinstance(args.dir_list, str):
        dirname_lis = [args.dir_list]
    else:
        dirname_lis = args.dir_list

    for d in dirname_lis:
        dir_fullname = os.path.abspath(d)
        get_polarizability(dir_fullname, calc_dir_patters=args.patterns,
                           deepth=args.deepth)


# TODO: Rb_mrci, d-aug-dyall.cv3z directory
