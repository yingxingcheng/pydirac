import os
import glob
import numpy as np
import warnings
from monty.os import cd
from pydirac.io.outputs import Output
from pydirac.analysis.polarizability import PolarizabilityCalculator
from pydirac.analysis.constant import *
from pydirac.core.settings import Settings

__all__ = [
    "get_polarizability",
    "get_polarizability_from_output_list",
    "do_one_basis",
    "is_valid",
    "get_atomDB",
]

"""
This script is used to collect all results from different type calculations,
e.g., CC or CI with different Hamiltonians and different calculation type, i.e.,
dipole and quadrupole. Using this results, one can generate atomic information
table and write it in a tex file.
"""


def get_polarizability(
    dirname: str = "./", calc_dir_patters=None, deepth: int = 0, verbos=True
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
        print("cd {0}".format(dirname))

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
                        "there are more than two output file in this directory {0}, and we take "
                        "the first one {1}".format(clc_d, outs[0])
                    )
                elif len(outs) == 1:
                    obj = Output(filename=outs[0])
                    if obj.is_ok:
                        calc_dir_output_lis.append(obj)
                else:
                    warnings.warn("there is no output file in this directory {0}".format(clc_d))

        all_res["curr_dir"] = {}
        if is_valid(curr_dir_output_lis):
            e_dict = get_polarizability_from_output_list(
                dirname, curr_dir_output_lis, tag="curr_dir", verbos=verbos
            )
            all_res["curr_dir"] = e_dict

        all_res["calc_dir"] = {}
        if is_valid(calc_dir_output_lis):
            e_dict = get_polarizability_from_output_list(
                dirname, calc_dir_output_lis, tag="clc_dir", verbos=verbos
            )
            all_res["calc_dir"] = e_dict

        # Step 3. deal with all sub_dir
        # ------------------------------
        all_res["sub_dir"] = {}
        if do_sub_dir:
            sub_dirs = [d for d in glob.glob("*") if os.path.isdir(d) and d not in clc_dirs]

            for sd in sub_dirs:
                res = get_polarizability(
                    os.path.join(dirname, sd), calc_dir_patters, deepth - 1, verbos=verbos
                )
                if sd not in all_res["sub_dir"] and len(res) > 0:
                    all_res["sub_dir"][os.path.basename(sd)] = res
    return all_res


def do_one_basis(output_lis):
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
                "we found a error output whose calc_type {0} is "
                "different with the first one {1}".format(o.inp.calc_type, ct)
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
                if not np.alltrue(np.asarray(i)):
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
                'The length of energy set "{0}" {1} does '
                "not equal the length of field {2}".format(k, len(v), len(fields))
            )
            del_k_lis.append(k)
    for k in del_k_lis:
        energies.pop(k)

    # {'scf': alpha_scf, 'mp2':alpha_scf, 'ccsd':alpha_ccsd}
    res = {}
    for k, energy in energies.items():
        polar_key = k.strip("_e")
        res[polar_key] = pc.get_svd_from_array(energy, fields)
    res_all = {"energy": {"energies": energies, "fields": fields}, "polar": res}
    return res_all


def get_polarizability_from_output_list(dirname, output_lis, tag=None, verbos=True):
    all_basis_res = {}
    for o in output_lis:
        task_type, orbit = o.task_type, o.calc_orbit
        obt_info = get_orbital_info(orbit["occ"], orbit["vir"]) if len(orbit) else "null"
        k = get_keyword(o.mol.molecule.atomic_info.symbol, task_type, obt_info)
        if k not in all_basis_res:
            all_basis_res[k] = []
        all_basis_res[k].append(o)

    e_res = {}
    # print("[Debug] get_polarizability_from_output_list():")
    for k, v in all_basis_res.items():
        one_res = do_one_basis(v)
        if one_res:
            # print(k, "==>", one_res)
            e_res[k] = one_res

    if verbos:
        for k, v in e_res.items():
            if tag:
                print('Table: results of {0} from \n "{1}" ({2})'.format(k, dirname, str(tag)))
            else:
                print('Table: results of {0} from \n "{1}"'.format(k, dirname))
            print("=" * 80)
            # if k.startswith('D'):
            if "@D" in k:
                print(
                    "  {0:<20s} {1:<20s} {2:<20s}".format(
                        "method", "polarizability", "hyper-polarizability"
                    )
                )
            else:
                print("  {0:<20s} {1:<20s} {2:<20s}".format("method", "momentum", "polarizability"))
                # In DIRAC, the correct quadrupole polarizability should be divided 4
                for _k in v.keys():
                    v[_k] = v[_k] / 4.0

            print("-" * 80)
            for i_k, i_v in v.items():
                print("  {0:<20s} {1:<20.3f} {2:<20.3f}".format(i_k.strip("_e"), i_v[0], i_v[1]))
            print("=" * 80)
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
            warnings.warn("the maximum of output objects with the same type is less than 3")
        return False


def get_atomDB(args):
    if isinstance(args.dir_list, str):
        dirname_lis = [args.dir_list]
    else:
        dirname_lis = args.dir_list

    for d in dirname_lis:
        dir_fullname = os.path.abspath(d)
        get_polarizability(dir_fullname, calc_dir_patters=args.patterns, deepth=args.deepth)
