import os

from os.path import join as opj

from pydirac.settings import Settings

__all__ = ['DiracJob']



class DiracJob:
    """A class representing a single computational job with DIRAC."""

    _result_type = DiracResults
    _top = ['dirac']
    _filenames = {'inp':'$JN.inp', 'run':'$JN.run', 'out':'$JN.out', 'err': '$JN.err'}


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
        SingleJob._get_ready(self)

