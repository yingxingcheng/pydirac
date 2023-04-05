import os

import yaml

__author__ = "YingXing Cheng"
__email__ = "yingxing.cheng@ugent.be"
__maintainer__ = "YingXing Cheng"
__maintainer_email__ = "yingxing.cheng@ugent.be"
__version__ = "1.0.0"

SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pmgrc.yaml")


def _load_pmg_settings():
    try:
        with open(SETTINGS_FILE) as f:
            d = yaml.safe_load(f)
    except OSError:
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

from pydirac.analysis import *
from pydirac.cli import *
from pydirac.core import *
from pydirac.io import *
