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

def package_files(directory, extensions):
    """Walk package directory to make sure we include all relevant files in package."""
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            if any([filename.endswith(ext) for ext in extensions]):
                paths.append(os.path.join('..', path, filename))
    return paths
json_yaml_csv_files = package_files('pydirac', ['yaml', 'json', 'csv'])

if __name__ == "__main__":
    print(module_dir)
    setup(
        name='pydirac',
        version='1.0.0',
        description='DIRAC interface',
        long_description=open(os.path.join(module_dir, 'README.md')).read(),
        url='https://github/yxcheng/pydirac',
        author='Yingxing Cheng',
        author_email='Yingxing.Cheng@ugent.be',
        license='GNU',
        packages=find_packages(),
        package_data={'pydirac': json_yaml_csv_files},
        zip_safe=False,
        install_requires=[],
        extras_require={},
        classifiers=['Programming Language :: Python :: 2.7',
                     "Programming Language :: Python :: 3",
                     "Programming Language :: Python :: 3.6",
                     'Development Status :: 5 - Production/Stable',
                     'Intended Audience :: Science/Research',
                     'Intended Audience :: System Administrators',
                     'Intended Audience :: Information Technology',
                     'Operating System :: OS Independent',
                     'Topic :: Other/Nonlisted Topic',
                     'Topic :: Scientific/Engineering'],
        #test_suite='nose.collector',
        #tests_require=['nose'],
        #scripts=[os.path.join('scripts', f) for f in os.listdir(os.path.join(module_dir, 'scripts'))]
        entry_points={
            'console_scripts': [
                'pypam = pydirac.cli.pypam:main',
            ]
        }
    )
