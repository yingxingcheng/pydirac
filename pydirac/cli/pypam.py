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
from pydirac.cli.pypam_input import get_inp
from pydirac.cli.pypam_atom_db import get_atomDB


def main():
    pypam_parser = argparse.ArgumentParser(description="pypam is a script tool to create, "
                                                 "run, analyse a DIRAC calculation "
                                                 "or a buntch of calculations")

    # input generation script
    subparsers = pypam_parser.add_subparsers()
    parser_input = subparsers.add_parser(
        "input", help="Tools for creating inputs for a DIRAC calculation")
    parser_input.add_argument('-c', '--calc_method', help="calculation method")
    parser_input.add_argument(
        '-d', '--dirname',
        help="where to generate all input files for a DIRAC calculation")
    parser_input.set_defaults(func=get_inp)

    # atomic database script
    parser_atomdb = subparsers.add_parser(
        "atomdb", help="Tools for compute polarizability "
                       "from a calculation directory")
    parser_atomdb.add_argument('-d', '--dirname', required=True)
    parser_atomdb.add_argument(
        '-t', '--sub_dir_tag', help="to specify which kind of "
                                    "directory should be considered",
        default='dyall')
    parser_atomdb.add_argument(
        '--deepth', type=int, default=0,
        help="how deep of directory with respect to current directory specified "
             "by '--dirname' to find the calculation information")

    parser_atomdb.set_defaults(func=get_atomDB)

    try:
        import argcomplete
        argcomplete.autocomplete(pypam_parser)
    except ImportError:
        # argcompolete not present
        pass

    args = pypam_parser.parse_args()
    #args = pypam_parser.parse_args(['input', '-d', './'])
    #args = pypam_parser.parse_args(['atomdb', '-d', './'])

    try:
        getattr(args, "func")
    except AttributeError:
        pypam_parser.print_help()
        sys.exit(-1)
    return args.func(args)


if __name__ == '__main__':
    main()
