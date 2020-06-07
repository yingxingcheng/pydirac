#!/usr/bin/env python

"""
Function: read and write DIRAC input files based on parameters specified by
user.
Author: Yingxing Cheng
Email: Yingxing.Cheng@ugent.be
Date: 10/21/2019
"""

from collections import OrderedDict
import re
from os.path import join as opj
from mendeleev import element
from pydirac.core.settings import Settings
from pydirac.utility.config import *
from pydirac.input.mole import get_mole_file
from pydirac.input.basic import input_from_calctype

__all__ = ['DiracJob', 'Inpobj', 'Molobj', 'OldDiracJob', 'JobType']

dir_path = os.path.dirname(os.path.abspath(__file__))


class Inpobj(object):
    """
    Input object for DIRAC inp file
    """
    first_level = {
        '**DIRAC': {'.TITLE': None, '.WAVE' : True, '.4INDEX': True, '.ANALYZE': True},
        '**ANALYZE': {'.MULPOP': True, },
        '**MULPOP' : '1..oo',
        '**HAMILTONIAN' : {'.X2C': True, '.NOSPIN': True, '.OPERATOR': {
            'ZDIPLEN': True, 'COMFACTOR': 0.0000}},
        '**INTEGRALS': {'*READINP': '.UNCONTRACT'},
        '**WAVE FUNCTIONS': {'.SCF': True, '.RELCCSD': True, '.RESOLVE': True}
        }

    second_level = {
    }


    def __init__(self, keyword_dict):
        self.origin_keyword_dict = keyword_dict
        self.keywords_list = list(self.origin_keyword_dict.items())

    def write_to_file(self, filename):
        """
        write context to file named by filename
        """
        # inp_str = []
        # for k1st, v1st in self.keyword_dict.items():
        #     inp_str.append(k1st)
        #     if len(v1st):
        #         inp_str.extend(v1st)

        inp_str = []
        for k, v in self.keywords_list:
            inp_str.append(k)
            if len(v):
                inp_str.extend(v)
        
        with open(filename, 'w') as f:
            f.write('\n'.join(inp_str))
        #print('\n'.join(inp_str))

    def add_keywords(self, keywords):
        """
        Add keywords which will be wrote to file
        """
        if not isinstance(keywords, tuple):
            raise TypeError('Wrong type for keywords, it should be tuple')
        
        # for k, v in self.keywords_list:
        #     for i in range(len(v)):
        #         v[i] = v[i].strip()

        if keywords[0] in ['**RELCC', '*KRCICALC', '*SCF']:
            for k,v in self.keywords_list:
                if k.startswith('**WAVE'):
                    # if keywords[0] in [_keywords.strip() for _keywords in v]:
                    #     break

                    str_added  = '#'
                    if keywords[0] == '*KRCICALC':
                        str_added = '.KR CI'
                    elif keywords[0] == '**RELCC':
                        str_added = '.RELCCSD'
                    else:
                        pass

                    for idx, kw in enumerate(v):
                        if kw.startswith('*'):
                            v.insert(idx, str_added)
                            break
                    break
                else:
                    # add '**WAVE'
                    pass

        for k, v in self.keywords_list:
            if k.startswith('**HAMI'):
                # change COMFACTOR
                for idx, value in enumerate(v):
                    if value.strip() in ['COMFACTOR', 'ZDIPLEN']:
                        v[idx] = ' '+value.strip()

                for idx, value in enumerate(v):
                    if value.strip() == 'COMFACTOR':
                        v[idx+1] = ' zff'
                        break
                break

        print('check resolved')
        for k, v in self.keywords_list:
            if k.startswith('**WAVE'):
                # check resolved
                if '.RESOLVE' in [_keywords.strip() for _keywords in v]:
                    break

                for idx, value in enumerate(v):
                    if value.startswith('*'):
                        v.insert(idx,'.RESOLVE')
                        break
                break
        
        print('check VECPOP')
        for k, v in self.keywords_list:
            if k.startswith('**ANAL'):
                # check .VECPOP
                for idx, value in enumerate(v):
                    if value.strip() == '.VECPOP':
                        v[idx+1] = '1..oo'
                        break
                break

        self.keywords_list.insert(-1, keywords)

    @classmethod
    def from_file(cls, filename):
        """
        create Inpobj from filename which can be inp file or output file of DIRAC
        """
        with open(filename, 'r') as f:
            lines = f.readlines()

        #start_line = re.compile(r"^Contents of the input file\s+$")
        #end_line = re.compile(r"^Contents of the molecule file\s+$")

        start_line = re.compile(r"^\*\*DIRAC\s*$")
        end_line = re.compile(r"^\*+END .*$")

        for i, l in enumerate(lines):
            if start_line.match(l):
                # print(l)
                break
        #print(i)
        lines = lines[i:]
        inp_lines = []
        for l in lines:
            if end_line.match(l):
                inp_lines.append(l)
                break

            if len(l.strip()):
                inp_lines.append(l)

        keyword_dict = OrderedDict()
        cur_1st_level = ''
        #cur_2nd_level = ''
        for l in inp_lines:
            #l = l.strip() # yxcheng 02/15/2020
            l = l.strip('\n')
            if l.startswith('**'):
                cur_1st_level = l
                if l not in keyword_dict.keys():
                    keyword_dict[l] = []
            # elif l.startswith('*'):
            #     cur_2nd_level = l
            #     if l not in keyword_dict[cur_1st_level].keys():
            #         keyword_dict[cur_1st_level][l] = []
            elif l.startswith('*END'):
                if l not in keyword_dict.keys():
                    keyword_dict[l] = []
            else:
               #keyword_dict[cur_1st_level][cur_2nd_level].append(l)
               keyword_dict[cur_1st_level].append(l)
        return Inpobj(keyword_dict)


class Molobj(object):

    def __init__(self):
        pass

    def write_to_file(self, filename):
        """write Molobj to a file named filename
        """
        pass

    @classmethod
    def from_file(cls, filename):
        pass


class JobType:
    d_DOSSSS_SCF = 0
    q_DOSSSS_SCF = 1
    d_DOSSSS_RELCCSD = 2
    q_DOSSSS_RELCCSD = 3
    d_X2C_SCF =  4
    q_X2C_SCF = 5
    d_X2C_NOSPIN_SCF = 6
    q_X2C_NOSPIN_SCF = 7
    d_X2C_NOSPIN_RELCCSD = 8
    q_X2C_NOSPIN_RELCCSD = 9



class OldDiracJob:

    @staticmethod
    def get_all_input(atom_type='He', basis_type='dyall.v2z',
                      basis_choice='BASIS',
                      field_value=None, create_script=True, scratch_dir=None,
                      suffix=None, calc_type=None):
        """Create input file according atom type and basis specified by users.
        """
        if scratch_dir is None:
            if suffix:
                scratch_dir = atom_type + suffix
            else:
                scratch_dir = atom_type
        elif suffix:
                basename = os.path.basename(scratch_dir)
                dirname = os.path.dirname(scratch_dir)
                scratch_dir = os.path.join(dirname, basename + "_"+ suffix)
        else:
            pass

        if not os.path.exists(scratch_dir):
            os.makedirs(scratch_dir)

        atom_info = element(atom_type)
        atom_type = atom_info.symbol

        fname = atom_type + '_' + basis_type + '.mol'
        fname = os.path.join(scratch_dir, fname)
        get_mole_file(atom_info=atom_info, filename_out=fname,
                      basis_type=basis_type, basis_choice=basis_choice)

        # create input file of DIRAC
        calc_type = calc_type or 0
        calc_funcs = ['d_DOSSSS_SCF', 'q_DOSSSS_SCF',
                      'd_DOSSSS_RELCCSD', 'q_DOSSSS_RELCCSD',
                      'd_X2C_SCF', 'q_X2C_SCF',
                      'd_X2C_NOSPIN_SCF', 'q_X2C_NOSPIN_SCF',
                      'd_X2C_NOSPIN_RELCCSD', 'd_X2C_NOSPIN_RELCCSD',
                      ]

        inp_fname = atom_type + '_' + basis_type + '.inp'
        inp_fname = os.path.join(scratch_dir, inp_fname)
        input_from_calctype(atom_info, calc_funcs[calc_type], inp_fname)

        # # create script to execute with different field values.
        # if create_script:
        #     field_str = []
        #     if field_value is None:
        #         field_str.append('+0.0000')
        #     else:
        #         for f in field_value:
        #             if f > 0:
        #                 field_str.append('+{0:.4f}'.format(f))
        #             else:
        #                 field_str.append('{0:.4f}'.format(f))

        #     field_str = ' '.join(field_str)
        #     # TODO: support multiple methods (post-HF)
        #     reference_str = '' if 0.00 in field_value else '+0.0000'

        #     calc_method = 'CCSD(T)'
        #     res_fname = '_'.join(['res', atom_type,
        #                           calc_method.replace('(', '\(').
        #                          replace(')', '\)'),
        #                           basis_type]) + '.dat'
        #     script_tp = get_script(atom_type, field_str, res_fname, calc_method,
        #                            field_value, dir_path, reference_str,
        #                            basis_type)

        #     spt_fname = 'run_{0}_{1}.sh'.format(atom_type, basis_type)
        #     spt_fname = os.path.join(scratch_dir, spt_fname)
        #     with open(spt_fname, 'w') as f:
        #         f.write(script_tp)

        #     # os.chmod(scratch_dir, 777)
        #     os.chmod(spt_fname, 0o777)


class DiracJob:
    """A class representing a single computational job with DIRAC."""

    _top = ['dirac']

    def __init__(self, settings):
        self.settings = settings


    def get_input(self):
        """Transform all contents of ``input`` branch of ``settings``
        into string with blocks, subblocks, keys and values.

        On the highest level alphabetic order of iteration is modified:
        keys occuring in class attribute ``_top`` are printed first.
        See :ref:`dirac-input` for details.
        """
        is_empty = lambda x: isinstance(x, Settings) and len(x) == 0

        def parse_key(key, value):
            ret = '.' + key.upper() + '\n'
            if not (value is True or is_empty(value)):
                if isinstance(value, list):
                    for i in value:
                        ret += str(i) + '\n'
                else:
                    ret += str(value) + '\n'
            return ret

        def parse_block(block):
            enabler = '_en'
            ret = '**' + block.upper() + '\n'
            s = self.settings.input[block]
            for k,v in s.items():
                if not isinstance(v, Settings) or is_empty(v):
                    ret += parse_key(k, v)
            for k,v in s.items():
                if isinstance(v, Settings) and enabler in v:
                    ret += parse_key(k, v[enabler])
            for k,v in s.items():
                if isinstance(v, Settings) and len(v) > 0:
                    ret += '*' + k.upper() + '\n'
                    for kk,vv in v.items():
                        if kk != enabler:
                            ret += parse_key(kk, vv)
            return ret

        inp = ''
        for block in self._top:
            if block in self.settings.input:
                inp += parse_block(block)
        for block in self.settings.input:
            if block not in self._top:
                inp += parse_block(block)
        inp += '*END OF INPUT\n'
        return inp


    def get_runscript(self):
        """Generate a runscript. Returned string is a ``pam`` call followed
        by option flags generated based on ``self.settings.runscript.pam``
        contents. See :ref:`dirac-runscript` for details."""
        r = self.settings.runscript.pam
        ret = 'pam'
        for k,v in r.items():
            ret += ' --%s'%k
            if v is not True:
                if isinstance(v, list):
                    ret += '="%s"' % ' '.join(v)
                else:
                    ret += '='+str(v)

        if self.settings.runscript.stdout_redirect:
            ret += ' >'+self._filename('out')
        ret += '\n\n'
        return ret

    def _get_ready(self):
        """Before generating runscript and input with parent method
        :meth:`SingleJob._get_ready<scm.plams.core.basejob.SingleJob._get_ready>`
        add proper ``mol`` and ``inp`` entries to ``self.settings.runscript.pam``.
        If already present there, ``mol`` will not be added.
        """
        s = self.settings.runscript.pam
        if 'mol' not in s:
            s.mol = self.name+'.xyz'
            with open(opj(self.path, self.name+'.xyz'), 'w') as f:
                f.write(str(len(self.molecule)) + '\n\n')
                for atom in self.molecule:
                    suffix = 'b={block}' if hasattr(atom,'block') else ''
                    f.write(atom.str(suffix=suffix)+'\n')
        s.inp = self._filename('inp')


if __name__ == '__main__':
    import sys

    argv = sys.argv[1:]

    relcc = ('**RELCC', ['.NEL_F1', '3 0 0 0 2 0 0 0', '.ENERGY',
                         '.PRINT', '1', '*CCENER',
                       '.MAXIT', '60', '.NTOL', '10'])
    # inactivate = (closed_elec - elec_in_gas)//2
    # gas1 = '10 12 / 6' # for elements from third row
    # gas2 = '{0} {1} / 3'.format(open_elec+10, open_elec+12) # for p open-shell orbit
    # gas3 = '{0} {1} / {2}'.format(open_elec+12, open_elec+12, nb_virtual)
    # krci = ('*KRCICALC',['.CI PROGRAM', 'LUCIAREL','.INACTIVE',inactivate,
    #                      '.GAS SHELLS',3, gas1,gas2,gas3,'.CIROOTS',ciroot1,'.CIROOTS',ciroot2,'.MAX CI','120', '.MXCIVE','60','.ANALYZ','.RSTRCI','rstr','.CHECKP'])
    for f in argv:
        inp = Inpobj.from_file(f)
        inp.add_keywords(relcc)
        # inp.add_keywords(krci)
        inp.write_to_file('tmp_B.inp')
