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

import argparse
import sys
from pydirac.cli.pypam_atom_db import get_atomDB
from pydirac.core.basis import get_custom_basis_from_ele


def get_inp():
    pass


def main():
    pypam_parser = argparse.ArgumentParser(
        description="pypam is a script tool to create, "
        "run, analyse a DIRAC calculation "
        "or a buntch of calculations"
    )

    # input generation script
    subparsers = pypam_parser.add_subparsers()
    parser_input = subparsers.add_parser(
        "input", help="Tools for creating inputs for a DIRAC calculation"
    )
    parser_input.add_argument("-c", "--calc_method", help="calculation method")
    parser_input.add_argument(
        "-d", "--dirname", help="where to generate all input files for a DIRAC calculation"
    )
    parser_input.set_defaults(func=get_inp)

    # atomic database script
    parser_atomdb = subparsers.add_parser(
        "atomdb", help="Tools for compute polarizability " "from a calculation directory"
    )
    parser_atomdb.add_argument(
        "dir_list", nargs="+", type=str, help="A list of directories to compute polarizability"
    )
    parser_atomdb.add_argument(
        "-p",
        "--patterns",
        nargs="+",
        type=str,
        help="to specify which kind of directory should be considered",
    )
    parser_atomdb.add_argument(
        "--deepth",
        type=int,
        default=0,
        help="how deep of directory with respect to current directory specified "
        "by '--dirname' to find the calculation information",
    )
    # parser_atomdb.add_argument(
    #     '--verbos', nargs=1, type=int, help="print level", default=1)

    parser_atomdb.set_defaults(func=get_atomDB)
    parser_atomdb.set_defaults(patterns=["dyall"])
    parser_atomdb.set_defaults(dir_list=["./"])

    # atomic basis set
    parser_basis = subparsers.add_parser(
        "basis",
        help="Tools for get explicit basis set based on mol file " "in which default basis is used",
    )
    parser_basis.add_argument("element_type", type=str, help="Element type")
    parser_basis.add_argument(
        "basis_type",
        type=str,
        help="Basis type. At present, this script only supports "
        "dyall.acv4z, dyall.acv3z, dyall.cv3z, and dyall.cv4z",
    )
    parser_basis.add_argument(
        "-f",
        "--filename_out",
        nargs="+",
        type=str,
        help="A list of directories to compute polarizability",
    )
    parser_basis.set_defaults(func=get_explicit_bs)

    try:
        import argcomplete

        argcomplete.autocomplete(pypam_parser)
    except ImportError:
        # argcompolete not present
        pass

    args = pypam_parser.parse_args()
    # args = pypam_parser.parse_args(['input', '-d', './'])
    # args = pypam_parser.parse_args(['atomdb', './', '-p', 'faegri', 'dyall'])
    # args = pypam_parser.parse_args(['atomdb', './'])
    # print(args)

    try:
        getattr(args, "func")
    except AttributeError:
        pypam_parser.print_help()
        sys.exit(-1)
    return args.func(args)


def get_explicit_bs(args):
    fout = args.filename_out
    get_custom_basis_from_ele(args.element_type, args.basis_type, fout)


if __name__ == "__main__":
    main()
