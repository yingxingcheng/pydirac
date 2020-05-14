#!/usr/bin/env python

import os

from setuptools import setup, find_packages

module_dir = os.path.dirname(os.path.abspath(__file__))

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
        license='modified BSD',
        packages=find_packages(),
        #package_data={'pygace': ['examples/*']},
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
    )
