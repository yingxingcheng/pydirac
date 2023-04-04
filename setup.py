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
from setuptools import setup, find_packages

module_dir = os.path.dirname(os.path.abspath(__file__))


def main():
    setup(
        name="pydirac",
        version="2023-04-04",
        description="DIRAC interface",
        long_description=open(os.path.join(module_dir, "README.md")).read(),
        url="https://github.com/yingxingcheng/pydirac",
        author="YingXing Cheng",
        author_email="yingxing.cheng@ugent.be",
        license="GNU",
        packages=find_packages(),
        include_package_data=True,
        package_data={
            "pydirac": ["data/**/*.*"],
            "pydirac.tests": ["data/**/*.*", "data/**/**/*.*", "data/**/**/**/*.*"],
        },
        zip_safe=False,
        install_requires=[],
        extras_require={},
        python_requires=">3.6",
        entry_points={"console_scripts": ["pypam = pydirac.cli.pypam:main"]},
    )


if __name__ == "__main__":
    main()
