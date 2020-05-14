#!/usr/bin/env python

"""
Function: read and write DIRAC input files based on parameters specified by user.
Author: Yingxing Cheng
Email: Yingxing.Cheng@ugent.be
Date: 10/21/2019
"""

from collections import OrderedDict
import re

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


if __name__ == '__main__':
    import sys

    argv = sys.argv[1:]

    relcc = ('**RELCC', ['.NEL_F1', '3 0 0 0 2 0 0 0', '.ENERGY', '.PRINT', '1', '*CCENER',
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
