#!/usr/bin/env python

import os
import sys
import subprocess
import re


def get_symbol_and_charge(fname='atom.mol'):
    with open(fname, 'r') as f:
        context = f.read()

    pattern = r'^\s+(\d+)\.\s+(\d+\.?\d?)\s+'
    re_obj = re.compile(pattern)
    atoms = re.findall(pattern, context, re.MULTILINE)
    print(atoms)
    return atoms


# @ Total CCSD(T) energy :                     -15.138231245354447
def get_energy(fname, method='CCSD(T)'):
    """
    method: CCSD(T), MP2, SCF
    """
    with open(fname, 'r') as f:
        context = f.read()

    if '(' and ')' in method:
        method = method.replace('(', '\(').replace(')', '\)')
    pattern = r'^@.*{method} energy[\s:]*(-?[\d\.]+)'.format(
        **{'method': method})

    # re_obj = re.compile(pattern)
    energy = re.findall(pattern, context, re.MULTILINE)
    print(energy)
    return energy



if __name__ == '__main__':
    pass
