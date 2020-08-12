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
import warnings
from fnmatch import fnmatch
import yaml


__author__ = "Yingxing Cheng"
__email__ = "yingxing.cheng@ugent.be"
__maintainer__ = "Yingxing Cheng"
__maintainer_email__ = "yingxing.cheng@ugent.be"
__version__ = "1.0.0"

SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pmgrc.yaml")


def _load_pmg_settings():
    try:
        with open(SETTINGS_FILE, "rt") as f:
            d = yaml.safe_load(f)
    except IOError:
        # If there are any errors, default to using environment variables
        # if present.
        d = {}
        # for k, v in os.environ.items():
        #     if k.startswith("PMG_"):
        #         d[k] = v
        #     elif k in ["VASP_PSP_DIR", "MAPI_KEY", "DEFAULT_FUNCTIONAL"]:
        #         d["PMG_" + k] = v
    d = d or {}
    return dict(d)


SETTINGS = _load_pmg_settings()


def loadfn(fname):
    """
    Convenience method to perform quick loading of data from a filename. The
    type of object returned depends the file type.

    Args:
        fname (string): A filename.

    Returns:
        Note that fname is matched using unix-style, i.e., fnmatch.
        (Structure) if *POSCAR*/*CONTCAR*/*.cif
        (Vasprun) *vasprun*
        (obj) if *json* (passthrough to monty.serialization.loadfn)
    """
    # if (fnmatch(fname, "*POSCAR*") or fnmatch(fname, "*CONTCAR*") or ".cif" in fname.lower()) or \
    #         fnmatch(fname, "*.vasp"):
    #     return Mol.from_file(fname)
    # if fnmatch(fname, "*vasprun*"):
    #     from pymatgen.io.vasp import Vasprun
    #     return Vasprun(fname)
    # if fnmatch(fname, "*.json*"):
    #     from monty.serialization import loadfn
    #     return loadfn(fname)
    # raise ValueError("Unable to determine how to process %s." % fname)
    pass
