#!/usr/bin/env python 

from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pickle
import os
import pandas as pd
from collections import OrderedDict

from pydirac.cli.pypam_atom_db import *
from pydirac.core.settings import Settings
from pydirac.core.periodic_table import Element
import warnings
warnings.filterwarnings("ignore")

#mpl.rcParams['font.family'] = 'Avenir'
plt.rc('font', family='serif')
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2

# red
red = '#E82229'
# blue
blue = '#2E378E'


class AbstractDataProcessor():
    best_dict = {
        'CC-NR': {
            'Li': ['CC-NR-2C@s-aug-ANO-RCC@(core 3)[vir 229]' ,],
            'Na': ['CC-NR-2C@s-aug-ANO-RCC@(core 11)[vir 287]',],
            'K' : ['CC-NR-2C@s-aug-ANO-RCC@(core 19)[vir 311]',],
            'Rb': ['CC-NR-2C@s-aug-ANO-RCC@(core 27)[vir 247]',],
            'Cs': ['CC-NR-2C@s-aug-ANO-RCC@(core 27)[vir 259]',],
            'Fr': ['CC-NR-2C@s-aug-ANO-RCC@(core 41)[vir 217]',],

            'Be': ['CC-NR-2C@d-aug-dyall.cv4z@(core 4)[vir 266]' , ],
            'Mg': ['CC-NR-2C@d-aug-dyall.cv4z@(core 10)[vir 314]', ],
            'Ca': ['CC-NR-2C@d-aug-dyall.cv4z@(core 20)[vir 510]', ],
            'Sr': ['CC-NR-2C@s-aug-dyall.cv4z@(core 28)[vir 346]', ],
            'Ba': ['CC-NR-2C@s-aug-dyall.cv4z@(core 28)[vir 398]', ],
            'Ra': ['CC-NR-2C@s-aug-dyall.cv4z@(core 42)[vir 378]', ],

            'Cu': ['CC-NR-2C@s-aug-dyall.cv4z@(core 29)[vir 299]', ],
            'Ag': ['CC-NR-2C@s-aug-dyall.cv4z@(core 29)[vir 355]', ],
            'Au': ['CC-NR-2C@s-aug-dyall.cv4z@(core 43)[vir 381]', ],

            'Zn': ['CC-NR-2C@s-aug-dyall.cv4z@(core 20)[vir 276]', ],
            'Cd': ['CC-NR-2C@s-aug-dyall.cv4z@(core 30)[vir 344]', ],
            'Hg': ['CC-NR-2C@s-aug-dyall.cv4z@(core 44)[vir 362]', ],
            'Cn': ['CC-NR-2C@s-aug-dyall.cv4z@(core 44)[vir 434]', ],

            'B' : ['CC-NR-2C@dyall.acv4z@(core 5)[vir 251]' , ],
            'Al': ['CC-NR-2C@dyall.acv4z@(core 13)[vir 363]', ],
            'Ga': ['CC-NR-2C@dyall.acv4z@(core 31)[vir 481]', ],
            'In': ['CC-NR-2C@dyall.acv4z@(core 49)[vir 553]', ],
            'Tl': ['CC-NR-2C@dyall.acv4z@(core 45)[vir 309]', ],
            'Nh': ['CC-NR-2C@dyall.acv4z@(core 45)[vir 305]', ],

            'C' : ['CC-NR-2C@dyall.acv4z@(core 6)[vir 250]' , ],
            'Si': ['CC-NR-2C@dyall.acv4z@(core 14)[vir 362]', ],
            'Ge': ['CC-NR-2C@dyall.acv4z@(core 32)[vir 480]', ],
            'Sn': ['CC-NR-2C@dyall.acv4z@(core 50)[vir 552]', ],
            'Pb': ['CC-NR-2C@dyall.acv4z@(core 46)[vir 308]', ],
            'Fl': ['CC-NR-2C@dyall.acv4z@(core 52)[vir 302]', ],

            'N' : ['CC-NR-2C@dyall.acv4z@(core 7)[vir 249]' , ],
            'P' : ['CC-NR-2C@dyall.acv4z@(core 15)[vir 361]', ],
            'As': ['CC-NR-2C@dyall.acv4z@(core 33)[vir 479]', ],
            'Sb': ['CC-NR-2C@dyall.acv4z@(core 51)[vir 551]', ],
            'Bi': ['CC-NR-2C@dyall.acv4z@(core 47)[vir 307]', ],
            'Mc': ['CC-NR-2C@dyall.acv4z@(core 47)[vir 303]', ],

            'O' : ['CC-NR-2C@dyall.acv4z@(core 8)[vir 248]' , ],
            'S' : ['CC-NR-2C@dyall.acv4z@(core 16)[vir 360]', ],
            'Se': ['CC-NR-2C@dyall.acv4z@(core 34)[vir 478]', ],
            'Te': ['CC-NR-2C@dyall.acv4z@(core 52)[vir 550]', ],
            'Po': ['CC-NR-2C@dyall.acv4z@(core 48)[vir 296]', ],
            'Lv': ['CC-NR-2C@dyall.acv4z@(core 48)[vir 294]', ],

            'F' : ['CC-NR-2C@dyall.acv4z@(core 9)[vir 247]' , ],
            'Cl': ['CC-NR-2C@dyall.acv4z@(core 17)[vir 359]', ],
            'Br': ['CC-NR-2C@dyall.acv4z@(core 35)[vir 477]', ],
            'I' : ['CC-NR-2C@dyall.acv4z@(core 53)[vir 549]', ],
            'At': ['CC-NR-2C@dyall.acv4z@(core 49)[vir 295]', ],
            'Ts': ['CC-NR-2C@dyall.acv4z@(core 49)[vir 293]', ],

            'He': ['CC-NR-2C@dyall.acv4z@(core 2)[vir 104]' , ],
            'Ne': ['CC-NR-2C@dyall.acv4z@(core 10)[vir 246]', ],
            'Ar': ['CC-NR-2C@dyall.acv4z@(core 18)[vir 358]', ],
            'Kr': ['CC-NR-2C@dyall.acv4z@(core 36)[vir 476]', ],
            'Xe': ['CC-NR-2C@dyall.acv4z@(core 26)[vir 286]', ],
            'Rn': ['CC-NR-2C@dyall.acv4z@(core 50)[vir 294]', ],
            'Og': ['CC-NR-2C@dyall.acv4z@(core 50)[vir 278]', ],
        },

        'CC-SR': {
            'Li': ['CC-SR-2C@s-aug-ANO-RCC@(core 3)[vir 229]' , 'CC-SR-2C@s-aug-ANO-RCC@(core 3)[vir 229]' ,],
            'Na': ['CC-SR-2C@s-aug-ANO-RCC@(core 11)[vir 287]', 'CC-SR-2C@s-aug-ANO-RCC@(core 11)[vir 287]',],
            'K' : ['CC-SR-2C@s-aug-ANO-RCC@(core 19)[vir 311]', 'CC-SR-2C@s-aug-ANO-RCC@(core 19)[vir 311]',],
            'Rb': ['CC-SR-2C@s-aug-ANO-RCC@(core 27)[vir 247]', 'CC-SR-2C@s-aug-ANO-RCC@(core 27)[vir 247]',],
            'Cs': ['CC-SR-2C@s-aug-ANO-RCC@(core 27)[vir 259]', 'CC-SR-2C@s-aug-ANO-RCC@(core 27)[vir 259]',],
            'Fr': ['CC-SR-2C@s-aug-ANO-RCC@(core 41)[vir 217]', 'CC-SR-2C@s-aug-ANO-RCC@(core 41)[vir 217]',],

            'Be': ['CC-SR-2C@d-aug-dyall.cv4z@(core 4)[vir 266]' , 'CC-SR-2C@d-aug-dyall.cv4z@(core 4)[vir 266]' ,],
            'Mg': ['CC-SR-2C@d-aug-dyall.cv4z@(core 10)[vir 314]', 'CC-SR-2C@d-aug-dyall.cv4z@(core 10)[vir 314]',],
            'Ca': ['CC-SR-2C@d-aug-dyall.cv4z@(core 20)[vir 510]', 'CC-SR-2C@d-aug-dyall.cv4z@(core 20)[vir 510]',], # (core 20)[vir 346]
            'Sr': ['CC-SR-2C@s-aug-dyall.cv4z@(core 28)[vir 346]', 'CC-SR-2C@s-aug-dyall.cv4z@(core 28)[vir 346]',],
            'Ba': ['CC-SR-2C@s-aug-dyall.cv4z@(core 28)[vir 398]', 'CC-SR-2C@s-aug-dyall.cv4z@(core 28)[vir 398]',],
            'Ra': ['CC-SR-2C@s-aug-dyall.cv4z@(core 42)[vir 378]', 'CC-SR-2C@s-aug-dyall.cv4z@(core 42)[vir 378]',], # (core 42)[vir 380]

            'Cu': ['CC-SR-2C@s-aug-dyall.cv4z@(core 29)[vir 299]', 'CC-SR-2C@s-aug-dyall.cv4z@(core 19)[vir 271]',],
            'Ag': ['CC-SR-2C@s-aug-dyall.cv4z@(core 29)[vir 355]', 'CC-SR-2C@s-aug-dyall.cv4z@(core 29)[vir 355]',],
            'Au': ['CC-SR-2C@s-aug-dyall.cv4z@(core 43)[vir 381]', 'CC-SR-2C@s-aug-dyall.cv4z@(core 33)[vir 381]',],

            'Zn': ['CC-SR-2C@s-aug-dyall.cv4z@(core 20)[vir 276]', 'CC-SR-2C@s-aug-dyall.cv4z@(core 20)[vir 276]',],
            'Cd': ['CC-SR-2C@s-aug-dyall.cv4z@(core 30)[vir 344]', 'CC-SR-2C@s-aug-dyall.cv4z@(core 30)[vir 344]',],
            'Hg': ['CC-SR-2C@s-aug-dyall.cv4z@(core 44)[vir 362]', 'CC-SR-2C@s-aug-dyall.cv4z@(core 44)[vir 362]',],
            'Cn': ['CC-SR-2C@s-aug-dyall.cv4z@(core 44)[vir 434]', 'CC-SR-2C@s-aug-dyall.cv4z@(core 44)[vir 434]',], # (core 48)[vir 434]

            'B' : ['CC-SR-2C@dyall.acv4z@(core 5)[vir 251]' , 'CC-SR-2C@dyall.acv4z@(core 5)[vir 251]' ],
            'Al': ['CC-SR-2C@dyall.acv4z@(core 13)[vir 363]', 'CC-SR-2C@dyall.acv4z@(core 13)[vir 363]'],
            'Ga': ['CC-SR-2C@dyall.acv4z@(core 31)[vir 481]', 'CC-SR-2C@dyall.acv4z@(core 13)[vir 203]'],
            'In': ['CC-SR-2C@dyall.acv4z@(core 49)[vir 553]', 'CC-SR-2C@dyall.acv4z@(core 31)[vir 309]'],
            'Tl': ['CC-SR-2C@dyall.acv4z@(core 45)[vir 309]', 'CC-SR-2C@dyall.acv4z@(core 35)[vir 309]'],
            'Nh': ['CC-SR-2C@dyall.acv4z@(core 45)[vir 305]', 'CC-SR-2C@dyall.acv4z@(core 35)[vir 305]'],

            'C' : ['CC-SR-2C@dyall.acv4z@(core 6)[vir 250]' , 'CC-SR-2C@dyall.acv4z@(core 6)[vir 250]' ,],
            'Si': ['CC-SR-2C@dyall.acv4z@(core 14)[vir 362]', 'CC-SR-2C@dyall.acv4z@(core 14)[vir 362]',],
            'Ge': ['CC-SR-2C@dyall.acv4z@(core 32)[vir 480]', 'CC-SR-2C@dyall.acv4z@(core 22)[vir 246]',],
            'Sn': ['CC-SR-2C@dyall.acv4z@(core 50)[vir 552]', 'CC-SR-2C@dyall.acv4z@(core 32)[vir 302]',],
            'Pb': ['CC-SR-2C@dyall.acv4z@(core 46)[vir 308]', 'CC-SR-2C@dyall.acv4z@(core 46)[vir 308]',],
            'Fl': ['CC-SR-2C@dyall.acv4z@(core 52)[vir 302]', 'CC-SR-2C@dyall.acv4z@(core 52)[vir 302]',], # (core 50)[vir 302]

            'N' : ['CC-SR-2C@dyall.acv4z@(core 7)[vir 249]' , 'CC-SR-2C@dyall.acv4z@(core 7)[vir 159]' ,],
            'P' : ['CC-SR-2C@dyall.acv4z@(core 15)[vir 361]', 'CC-SR-2C@dyall.acv4z@(core 15)[vir 185]',],
            'As': ['CC-SR-2C@dyall.acv4z@(core 33)[vir 479]', 'CC-SR-2C@dyall.acv4z@(core 15)[vir 185]',],
            'Sb': ['CC-SR-2C@dyall.acv4z@(core 51)[vir 551]', 'CC-SR-2C@dyall.acv4z@(core 15)[vir 227]',],
            'Bi': ['CC-SR-2C@dyall.acv4z@(core 47)[vir 307]', 'CC-SR-2C@dyall.acv4z@(core 15)[vir 265]',],
            'Mc': ['CC-SR-2C@dyall.acv4z@(core 47)[vir 303]', 'CC-SR-2C@dyall.acv4z@(core 15)[vir 265]',],

            'O' : ['CC-SR-2C@dyall.acv4z@(core 8)[vir 248]' , 'CC-SR-2C@dyall.acv4z@(core 8)[vir 186]' ,],
            'S' : ['CC-SR-2C@dyall.acv4z@(core 16)[vir 360]', 'CC-SR-2C@dyall.acv4z@(core 16)[vir 184]',],
            'Se': ['CC-SR-2C@dyall.acv4z@(core 34)[vir 478]', 'CC-SR-2C@dyall.acv4z@(core 16)[vir 184]',],
            'Te': ['CC-SR-2C@dyall.acv4z@(core 52)[vir 550]', 'CC-SR-2C@dyall.acv4z@(core 16)[vir 216]',],
            'Po': ['CC-SR-2C@dyall.acv4z@(core 48)[vir 296]', 'CC-SR-2C@dyall.acv4z@(core 16)[vir 262]',],
            'Lv': ['CC-SR-2C@dyall.acv4z@(core 48)[vir 294]', 'CC-SR-2C@dyall.acv4z@(core 16)[vir 262]',],

            'F' : ['CC-SR-2C@dyall.acv4z@(core 9)[vir 247]' , 'CC-SR-2C@dyall.acv4z@(core 9)[vir 183]' ,], # (core 9)[vir 115]
            'Cl': ['CC-SR-2C@dyall.acv4z@(core 17)[vir 359]', 'CC-SR-2C@dyall.acv4z@(core 17)[vir 169]',],
            'Br': ['CC-SR-2C@dyall.acv4z@(core 35)[vir 477]', 'CC-SR-2C@dyall.acv4z@(core 17)[vir 183]',],
            'I' : ['CC-SR-2C@dyall.acv4z@(core 53)[vir 549]', 'CC-SR-2C@dyall.acv4z@(core 17)[vir 215]',],
            'At': ['CC-SR-2C@dyall.acv4z@(core 49)[vir 295]', 'CC-SR-2C@dyall.acv4z@(core 17)[vir 261]',],
            'Ts': ['CC-SR-2C@dyall.acv4z@(core 49)[vir 293]', 'CC-SR-2C@dyall.acv4z@(core 17)[vir 247]',],

            'He': ['CC-SR-2C@dyall.acv4z@(core 2)[vir 104]' , 'CC-SR-2C@dyall.acv4z@(core 2)[vir 104]' ,],
            'Ne': ['CC-SR-2C@dyall.acv4z@(core 10)[vir 246]', 'CC-SR-2C@dyall.acv4z@(core 10)[vir 246]',],
            'Ar': ['CC-SR-2C@dyall.acv4z@(core 18)[vir 358]', 'CC-SR-2C@dyall.acv4z@(core 18)[vir 358]',],
            'Kr': ['CC-SR-2C@dyall.acv4z@(core 36)[vir 476]', 'CC-SR-2C@dyall.acv4z@(core 36)[vir 476]',],
            'Xe': ['CC-SR-2C@dyall.acv4z@(core 26)[vir 286]', 'CC-SR-2C@dyall.acv4z@(core 26)[vir 286]',], # (core 50)[vir 286]
            'Rn': ['CC-SR-2C@dyall.acv4z@(core 50)[vir 294]', 'CC-SR-2C@dyall.acv4z@(core 50)[vir 294]',],
            'Og': ['CC-SR-2C@dyall.acv4z@(core 50)[vir 278]', 'CC-SR-2C@dyall.acv4z@(core 50)[vir 278]',],
        },
        'CC-SO': {
            'Li': ['CC-SO-4C@s-aug-ANO-RCC@(core 3)[vir 229]' , ],
            'Na': ['CC-SO-4C@s-aug-ANO-RCC@(core 11)[vir 287]', ],
            'K' : ['CC-SO-4C@s-aug-ANO-RCC@(core 19)[vir 311]', ],
            'Rb': ['CC-SO-4C@s-aug-ANO-RCC@(core 27)[vir 247]', ],
            'Cs': ['CC-SO-4C@s-aug-ANO-RCC@(core 27)[vir 259]', ],
            'Fr': ['CC-SO-4C@s-aug-ANO-RCC@(core 41)[vir 217]', ],

            'Be': ['CC-SO-4C@d-aug-dyall.cv4z@(core 4)[vir 266]' , ],
            'Mg': ['CC-SO-4C@d-aug-dyall.cv4z@(core 10)[vir 314]', ],
            'Ca': ['CC-SO-4C@d-aug-dyall.cv4z@(core 20)[vir 346]', ],
            'Sr': ['CC-SO-2C@s-aug-dyall.cv4z@(core 28)[vir 346]', ],
            'Ba': ['CC-SO-2C@s-aug-dyall.cv4z@(core 28)[vir 398]', ],
            'Ra': ['CC-SO-2C@s-aug-dyall.cv4z@(core 42)[vir 380]', ],

            'Cu': ['CC-SO-4C@s-aug-dyall.cv4z@(core 19)[vir 271]', ],
            'Ag': ['CC-SO-4C@s-aug-dyall.cv4z@(core 29)[vir 355]', ],
            'Au': ['CC-SO-4C@s-aug-dyall.cv4z@(core 33)[vir 381]', ],

            'Zn': ['CC-SO-2C@s-aug-dyall.cv4z@(core 20)[vir 276]', ],
            'Cd': ['CC-SO-2C@s-aug-dyall.cv4z@(core 30)[vir 344]', ],
            'Hg': ['CC-SO-2C@s-aug-dyall.cv4z@(core 44)[vir 362]', ],
            'Cn': ['CC-SO-2C@s-aug-dyall.cv4z@(core 48)[vir 434]', ],

            'Al': ['CC-SO-4C@dyall.acv4z@(core 13)[vir 363]', ],
            'Ga': ['CC-SO-4C@dyall.acv4z@(core 13)[vir 203]', ],
            'In': ['CC-SO-2C@dyall.acv4z@(core 31)[vir 309]', ],
            'Tl': ['CC-SO-4C@dyall.acv4z@(core 35)[vir 309]', ],
            'Nh': ['CC-SO-4C@dyall.acv4z@(core 35)[vir 305]', ],

            'C' : ['CC-SO-4C@dyall.acv4z@(core 6)[vir 250]' , ],
            'Si': ['CC-SO-2C@dyall.acv4z@(core 14)[vir 362]', ],
            'Ge': ['CC-SO-2C@dyall.acv4z@(core 22)[vir 246]', ],
            'Sn': ['CC-SO-2C@dyall.acv4z@(core 32)[vir 302]', ],
            'Pb': ['CC-SO-2C@dyall.acv4z@(core 46)[vir 308]', ],
            'Fl': ['CC-SO-2C@dyall.acv4z@(core 50)[vir 302]', ],

            'He': ['CC-SO-4C@dyall.acv4z@(core 2)[vir 104]' , ],
            'Ne': ['CC-SO-4C@dyall.acv4z@(core 10)[vir 246]', ],
            'Ar': ['CC-SO-4C@dyall.acv4z@(core 18)[vir 358]', ],
            'Kr': ['CC-SO-2C@dyall.acv4z@(core 36)[vir 476]', ],
            'Xe': ['CC-SO-4C@dyall.acv4z@(core 50)[vir 286]', ],
            'Rn': ['CC-SO-4C@dyall.acv4z@(core 50)[vir 294]', ],
            'Og': ['CC-SO-4C@dyall.acv4z@(core 50)[vir 278]', ],
        },
        'CI-SO': {
            'Li': ['CI-SO-4C@s-aug-ANO-RCC@(core 3)[vir 209]' , ],
            'Na': ['CI-SO-4C@s-aug-ANO-RCC@(core 9)[vir 245]' , ],
            'K' : ['CI-SO-4C@s-aug-ANO-RCC@(core 19)[vir 199]', ],
            'Rb': ['CI-SO-4C@s-aug-ANO-RCC@(core 19)[vir 231]', ],
            'Cs': ['CI-SO-4C@s-aug-ANO-RCC@(core 19)[vir 241]', ],
            'Fr': ['CI-SO-4C@s-aug-ANO-RCC@(core 19)[vir 187]', ],

            'Be': ['CI-SO-4C@dyall.cv4z@(core 4)[vir 132]' , ],
            'Mg': ['CI-SO-4C@dyall.cv4z@(core 10)[vir 182]', ],
            'Ca': ['CI-SO-4C@dyall.cv4z@(core 20)[vir 228]', ],
            'Sr': ['CI-SO-4C@dyall.cv4z@(core 20)[vir 214]', ],
            'Ba': ['CI-SO-4C@dyall.cv4z@(core 20)[vir 252]', ],
            'Ra': ['CI-SO-4C@dyall.cv4z@(core 20)[vir 234]', ],

            'Cu': ['CI-SO-4C@dyall.cv4z@(core 11)[vir 187]', ],
            'Ag': ['CI-SO-4C@dyall.cv4z@(core 11)[vir 223]', ],
            'Au': ['CI-SO-4C@dyall.cv4z@(core 11)[vir 241]', ],

            'Zn': ['CI-SO-4C@dyall.cv4z@(core 12)[vir 154]', ],
            'Cd': ['CI-SO-4C@dyall.cv4z@(core 12)[vir 200]', ],
            'Hg': ['CI-SO-4C@dyall.cv4z@(core 12)[vir 240]', ],
            'Cn': ['CI-SO-4C@dyall.cv4z@(core 12)[vir 264]', ],

            'B' : ['CI-SO-4C@dyall.acv4z@(core 5)[vir 251]' , ],
            'Al': ['CI-SO-4C@dyall.acv4z@(core 13)[vir 223]', ],
            'Ga': ['CI-SO-4C@dyall.acv4z@(core 13)[vir 203]', ],
            'In': ['CI-SO-4C@dyall.acv4z@(core 13)[vir 237]', ],
            'Tl': ['CI-SO-4C@dyall.acv4z@(core 13)[vir 277]', ],
            'Nh': ['CI-SO-4C@dyall.acv4z@(core 13)[vir 271]', ],

            'C' : ['CI-SO-4C@dyall.acv4z@(core 6)[vir 162]' , ],
            'Si': ['CI-SO-4C@dyall.acv4z@(core 14)[vir 204]', ],
            'Ge': ['CI-SO-4C@dyall.acv4z@(core 14)[vir 186]', ],
            'Sn': ['CI-SO-4C@dyall.acv4z@(core 14)[vir 230]', ],
            'Pb': ['CI-SO-4C@dyall.acv4z@(core 14)[vir 276]', ],
            'Fl': ['CI-SO-4C@dyall.acv4z@(core 14)[vir 270]', ],

            'N' : ['CI-SO-4C@dyall.acv4z@(core 7)[vir 159]' , ],
            'P' : ['CI-SO-4C@dyall.acv4z@(core 15)[vir 185]', ],
            'As': ['CI-SO-4C@dyall.acv4z@(core 15)[vir 185]', ],
            'Sb': ['CI-SO-4C@dyall.acv4z@(core 15)[vir 227]', ],
            'Bi': ['CI-SO-4C@dyall.acv4z@(core 15)[vir 265]', ],
            'Mc': ['CI-SO-4C@dyall.acv4z@(core 15)[vir 265]', ],

            'O' : ['CI-SO-4C@dyall.acv4z@(core 8)[vir 116]' , ],
            'S' : ['CI-SO-4C@dyall.acv4z@(core 16)[vir 184]', ],
            'Se': ['CI-SO-4C@dyall.acv4z@(core 16)[vir 184]', ],
            'Te': ['CI-SO-4C@dyall.acv4z@(core 16)[vir 216]', ],
            'Po': ['CI-SO-4C@dyall.acv4z@(core 16)[vir 262]', ],
            'Lv': ['CI-SO-4C@dyall.acv4z@(core 16)[vir 262]', ],

            'F' : ['CI-SO-4C@dyall.acv4z@(core 9)[vir 115]' , ],
            'Cl': ['CI-SO-4C@dyall.acv4z@(core 17)[vir 169]', ],
            'Br': ['CI-SO-4C@dyall.acv4z@(core 17)[vir 183]', ],
            'I' : ['CI-SO-4C@dyall.acv4z@(core 17)[vir 215]', ],
            'At': ['CI-SO-4C@dyall.acv4z@(core 17)[vir 261]', ],
            'Ts': ['CI-SO-4C@dyall.acv4z@(core 17)[vir 247]', ],

            'He': ['CI-SO-4C@dyall.acv4z@(core 2)[vir 60]'  , ],
            'Ne': ['CI-SO-4C@dyall.acv4z@(core 10)[vir 108]', ],
            'Ar': ['CI-SO-4C@dyall.acv4z@(core 18)[vir 156]', ],
            'Kr': ['CI-SO-4C@dyall.acv4z@(core 18)[vir 182]', ],
            'Xe': ['CI-SO-4C@dyall.acv4z@(core 18)[vir 214]', ],
            'Rn': ['CI-SO-4C@dyall.acv4z@(core 18)[vir 246]', ],
            'Og': ['CI-SO-4C@dyall.acv4z@(core 18)[vir 246]', ],
        },
    }

    def __init__(self):
        pass

    def __str__(self):
        raise NotImplementedError

    @staticmethod
    def get_new_ct(old_ct):
        type_dict = {
            'CC-NR-2C': '2C-NR-CC',
            'CC-SR-2C': '2C-SR-CC',
            'CC-SO-2C': '2C-DC-CC',
            'CC-SO-4C': '4C-DC-CC',
            'CI-SO-4C': '4C-DC-CI',
        }

        new_ct = old_ct
        for k, v in type_dict.items():
            if k in old_ct:
                new_ct = old_ct.replace(k, v)

        return new_ct

    @staticmethod
    def is_the_best(calc_type, ele_symbol):
        for key in ['CC-NR', 'CC-SR', 'CC-SO', 'CI-SO']:
            if key in calc_type:
                if ele_symbol in AbstractDataProcessor.best_dict[key] and \
                        calc_type in AbstractDataProcessor.best_dict[key][
                    ele_symbol]:
                    # AbstractDataProcessor.best_dict[key][ele_symbol][0] == calc_type:
                    return True
                else:
                    return False
        else:
            return False

    @staticmethod
    def two_SRCC(calc_type, ele_symbol):
        for key in ['CC-NR', 'CC-SR', 'CC-SO', 'CI-SO']:
            if key in calc_type:
                if ele_symbol in AbstractDataProcessor.best_dict[key] and \
                        calc_type in AbstractDataProcessor.best_dict[key][ele_symbol]:
                    return True
                else:
                    return False
        else:
            return False


class SRCCDataProcessor(AbstractDataProcessor):
    def __init__(self, dirname, has_header=False, only_best=False):
        super(SRCCDataProcessor, self).__init__()

        self.dirname = dirname
        self.has_header = has_header
        self.only_best = only_best

        self.res_ml0 = None
        self.res_ml1 = None

        for p in Path(self.dirname).glob('Ml*'):
            if p.is_dir() and p.name == 'Ml0':
                self.res_ml0 = get_polarizability(
                    str(p), ['dyall', 'ANO-RCC', 'faegri'], verbos=False)
            elif p.is_dir() and p.name == 'Ml1':
                self.res_ml1 = get_polarizability(
                    str(p), ['dyall', 'ANO-RCC', 'faegri'],verbos=False)

        self._extract_info()


    def _extract_info(self):
        self.all_res = Settings()
        ele_type = None
        self.is_quadru = False

        if self.res_ml1:
            # we assume that res_ml0 is not empty
            assert self.res_ml0 is not None

            f_ml0 = Settings(self.res_ml0).flatten()
            f_ml1 = Settings(self.res_ml1).flatten()

            ave = Settings()
            for k in f_ml0:
                if k in f_ml1:
                    if len(f_ml0[k]) and len(f_ml0[k]):
                        ave[k] = (np.asarray(f_ml0[k]) + np.asarray(
                            f_ml1[k]) * 2) / 3.
                    else:
                        ave[k] = Settings()

            for k, v in self.res_ml0.items():
                if k in ['curr_dir', 'calc_dir'] and len(v) > 0:
                    for kk, vv in v.items():
                        # kk is calc_type, e.g., Li@D-CC-2C-SR@ANO-RCC
                        # then, we collect data for different calc_type, i.e., kk
                        ele_type = kk.split('@')[0].strip()
                        if '@Q' in kk:
                            self.is_quadru = True
                        self.all_res[kk] = Settings()
                        self.all_res[kk].ml0 = vv
                        self.all_res[kk].ml1 = self.res_ml1[k][kk]
                        self.all_res[kk].ave = ave.unflatten()[k][kk]

        elif self.res_ml0:
            for k, v in self.res_ml0.items():
                if k in ['curr_dir', 'calc_dir'] and len(v) > 0:
                    for kk, vv in v.items():
                        # kk is calc_type, e.g., Li@D-CC-2C-SR@ANO-RCC
                        # then, we collect data for different calc_type, i.e., kk
                        ele_type = kk.split('@')[0].strip()
                        if '@Q' in kk:
                            self.is_quadru = True
                        self.all_res[kk] = Settings()
                        self.all_res[kk].ml0 = vv
                        self.all_res[kk].ml1 = Settings()
                        self.all_res[kk].ave = vv
        else:
            raise RuntimeError('res_ml0 is None!')

        self.element = Element(ele_type)

        if self.element.group in [1, 11]:
            self.gs_term = r'^2S'
        elif self.element.group in [2, 12, 18]:
            self.gs_term = r'^1S'
        elif self.element.group in [13]:
            self.gs_term = r'^2P'
        elif self.element.group in [14, 16]:
            self.gs_term = r'^3P'
        elif self.element.group in [15]:
            self.gs_term = r'^4S'
        elif self.element.group in [17]:
            self.gs_term = r'^2P'
        else:
            raise RuntimeError('Unsupport group {0}'.format(self.element.group))

        if self.is_quadru:
            self.dq = 'quadrupole'
            self.dq_tag = 'Q'
        else:
            self.dq = 'dipole'
            self.dq_tag = 'D'


    def __str__(self):
        header, body, footer = self.get_tex()
        return '\n'.join(header + body + footer)

    def get_tex(self):
        long_table_header = r"""
\begin{{longtable}}{{ccccccc}}
\caption{{The results of {DQ} polarizability of {ele_type}}} \\
\hline 
\multicolumn{{1}}{{c}}{{\textbf{{Z}}}} & 
\multicolumn{{1}}{{c}}{{\textbf{{Atom}}}} & 
\multicolumn{{1}}{{c}}{{\textbf{{State}}}} & 
\multicolumn{{1}}{{c}}{{\textbf{{$\alpha_{dq_tag}$ (a.u.)}}}} & 
\multicolumn{{1}}{{c}}{{\textbf{{$\delta$ (\%)}}}} & 
\multicolumn{{1}}{{c}}{{\textbf{{Method}}}} & 
\multicolumn{{1}}{{c}}{{\textbf{{Comments}}}} \\ \hline
\endfirsthead

\multicolumn{{7}}{{l}}
{{{{\bfseries \tablename\ \thetable{{}}. continued.}}}} \\ \hline 
\multicolumn{{1}}{{c}}{{\textbf{{Z}}}} & 
\multicolumn{{1}}{{c}}{{\textbf{{Atom}}}} & 
\multicolumn{{1}}{{c}}{{\textbf{{State}}}} & 
\multicolumn{{1}}{{c}}{{\textbf{{$\alpha_D$ (a.u.)}}}} & 
\multicolumn{{1}}{{c}}{{\textbf{{$\delta$}}}} & 
\multicolumn{{1}}{{c}}{{\textbf{{Method}}}} & 
\multicolumn{{1}}{{c}}{{\textbf{{Comments}}}} \\ \hline
\endhead
\hline \multicolumn{{7}}{{r}}{{{{(continued)}}}} \\
\endfoot
\endlastfoot
"""

        header = r"""\begin{{table}}[H]
    \centering
    \caption{{The results of {DQ} polarizability of {ele_type}}}
        \begin{{tabular}}{{ccccccc}}
        \hline\hline
        Z & Atom & State & $\alpha_{dq_tag}$ (a.u.) & $\delta$ (\%) & Method & Comments   \\   \hline
"""

        # body = r"""         {ele_Z} & {ele_symbol} & \multicolumn{{5}}{{r}}{{{calc_type}}} \\
        #      &      & ${gs_term}, M_L=0$        & {ml0_scf:.2f}  & {pe_scf_ml0:.2f} & DHF     &            \\
        body = r"""         {has_line}{ele_Z} & {ele_symbol}      & ${gs_term}, M_L=0$        & {ml0_scf:.2f}     & {pe_scf_ml0:.2f} & DHF     & {calc_type}\\
          &      & ${gs_term}, M_L=\pm1$     & {ml1_scf:.2f}    & {pe_scf_ml1:.2f} &        &            \\
          &      & ${gs_term}$    & {ave_scf:.2f}    & {pe_scf:.2f} &       & $\Delta \alpha = {delta_scf:.2f}$             \\
          &      & ${gs_term}, M_L=0$        & {ml0_mp2:.2f}   & {pe_mp2_ml0:.2f} & MP2     &             \\
          &      & ${gs_term}, M_L=\pm1$     & {ml1_mp2:.2f}   & {pe_mp2_ml1:.2f} &         &            \\
          &      & ${gs_term}$    & {ave_mp2:.2f}    & {pe_mp2:.2f} &         & $\Delta \alpha = {delta_mp2:.2f}$             \\
          &      & ${gs_term}, M_L=0$        & {ml0_ccsd:.2f}  & {pe_ccsd_ml0:.2f} & CCSD    &             \\
          &      & ${gs_term}, M_L=\pm1$     & {ml1_ccsd:.2f}  & {pe_ccsd_ml1:.2f} &         &            \\
          &      & ${gs_term}$    & {ave_ccsd:.2f}    & {pe_ccsd:.2f} &         & $\Delta \alpha = {delta_ccsd:.2f}$              \\
          &      & ${gs_term}, M_L=0$        & {ml0_ccsd_t:.2f}  &   & CCSD(T) &              \\
          &      & ${gs_term}, M_L=\pm1$     & {ml1_ccsd_t:.2f}  &   &         &             \\
          &      & ${gs_term}$    & {ave_ccsd_t:.2f}    &   &         & $\Delta \alpha = {delta_ccsd_t:.2f}$            \\
"""
        # body_q = r"""         {ele_Z} & {ele_symbol} & \multicolumn{{5}}{{r}}{{{calc_type}}} \\
        #      &      & ${gs_term}, M_L=0$        & {ml0_scf:.2f}  & {pe_scf_ml0:.2f} & DHF     &            \\
        body_q = r"""         {has_line}{ele_Z} & {ele_symbol} & ${gs_term}, M_L=0$        & {ml0_scf:.2f}  & {pe_scf_ml0:.2f} & DHF     & {calc_type}           \\
          &      & ${gs_term}, M_L=\pm1$     & {ml1_scf:.2f}  & {pe_scf_ml1:.2f} &         &            \\
          &      & ${gs_term}     & {ave_scf:.2f}    & {pe_scf:.2f} &         &  $\Delta \alpha = {delta_scf:.2f}$             \\
          &      & ${gs_term}, M_L=0$        & {ml0_ccsd:.2f}    &  & CCSD    &             \\
          &      & ${gs_term}, M_L=\pm1$     & {ml1_ccsd:.2f}    &  &         &             \\
          &      & ${gs_term}     & {ave_ccsd:.2f}    &  &         &  $\Delta \alpha = {delta_ccsd:.2f}$            \\
"""

        # body_only_ml0 = r"""         {ele_Z} & {ele_symbol}  & \multicolumn{{5}}{{r}}{{{calc_type}}}  \\
        #      &      & ${gs_term}, M_L=0$        & {ml0_scf:.2f}  & {pe_scf_ml0:.2f} & DHF     & \\
        body_only_ml0 = r"""         {has_line}{ele_Z} & {ele_symbol}  & ${gs_term}, M_L=0$        & {ml0_scf:.2f}  & {pe_scf_ml0:.2f} & DHF     & {calc_type}\\
          &      &                & {ml0_mp2:.2f}    & {pe_mp2_ml0:.2f} & MP2     &            \\
          &      &                & {ml0_ccsd:.2f}    & {pe_ccsd_ml0:.2f} & CCSD    &            \\
          &      &                & {ml0_ccsd_t:.2f}  &   & CCSD(T) &             \\
"""
        # body_only_ml0_q = r"""        {ele_Z} & {ele_symbol}  & \multicolumn{{4}}{{r}}{{{calc_type}}} \\
        #      &      & ${gs_term}, M_L=0$ &        {ml0_scf:.2f}     & {pe_scf_ml0:.2f} & DHF     &  \\
        body_only_ml0_q = r"""        {has_line}{ele_Z} & {ele_symbol}  & ${gs_term}, M_L=0$ &        {ml0_scf:.2f}     & {pe_scf_ml0:.2f} & DHF     & {calc_type}  \\
          &      &                & {ml0_ccsd:.2f}    & CCSD    &             \\
"""

        footer = r"""    \hline\hline
    \end{{tabular}}
    \label{{tab:{DQ}_{ele_type}}}
\end{{table}}
"""
        long_table_footer = r"""    \hline
    \label{{tab:{DQ}_{ele_type}}}
\end{{longtable}}
"""
        # return CCsr_tex(self.res_ml0, self.res_ml1,
        #                 has_header=self.has_header,
        #                 only_best=self.only_best)
        ele = self.element
        ele_Z = ele.Z
        ele_symbol = ele.symbol

        if self.is_quadru:
            dq = 'quadrupole'
            dq_tag = 'Q'
        else:
            dq = 'dipole'
            dq_tag = 'D'

        if self.has_header:
            # template = header.format(**{'DQ':dq, 'dq_tag':dq_tag, 'ele_type': ele_symbol})
            header_tex = long_table_header.format(**{'DQ':dq, 'dq_tag':dq_tag, 'ele_type': ele_symbol})
            has_line = ''
        else:
            header_tex = ''
            has_line = '\cline{3-7}\n         '

        body_tex = ''
        if self.res_ml1:
            for ct, v in sorted(self.all_res.items(), reverse=True):
                if not self.is_quadru:
                    if not self.has_header:
                        ele_Z = ''
                        ele_symbol = '   '

                    tmp_ct = ct.replace(ele.symbol + '@D-', '')

                    if self.only_best and not SRCCDataProcessor.is_the_best(tmp_ct, ele.symbol):
                        continue

                    body_tex += body.format(**{
                        'ml0_scf':    round(v.ml0.scf_e[0], 2),
                        'ml0_mp2':    round(v.ml0.mp2_e[0], 2),
                        'ml0_ccsd':   round(v.ml0.ccsd_e[0], 2),
                        'ml0_ccsd_t': round(v.ml0.ccsd_p_T_e[0], 2),
                        'ml1_scf':    round(v.ml1.scf_e[0], 2),
                        'ml1_mp2':    round(v.ml1.mp2_e[0], 2),
                        'ml1_ccsd':   round(v.ml1.ccsd_e[0], 2),
                        'ml1_ccsd_t': round(v.ml1.ccsd_p_T_e[0], 2),
                        'ave_scf':    round(v.ave.scf_e[0], 2),
                        'ave_mp2':    round(v.ave.mp2_e[0], 2),
                        'ave_ccsd':   round(v.ave.ccsd_e[0], 2),
                        'ave_ccsd_t': round(v.ave.ccsd_p_T_e[0], 2),
                        'calc_type': SRCCDataProcessor.get_new_ct(tmp_ct),
                        'ele_Z' : ele_Z,
                        'ele_symbol': ele_symbol,
                        'gs_term': self.gs_term,
                        'delta_mp2':   round( v.ml1.mp2_e[0] - v.ml0.mp2_e[0], 2),
                        'delta_scf':   round( v.ml1.scf_e[0] - v.ml0.scf_e[0], 2),
                        'delta_ccsd':  round( v.ml1.ccsd_e[0] - v.ml0.ccsd_e[0], 2),
                        'delta_ccsd_t':round( v.ml1.ccsd_p_T_e[0] - v.ml0.ccsd_p_T_e[0], 2),
                        'pe_scf_ml0':  round((v.ml0.scf_e[0] - v.ml0.ccsd_p_T_e[0])/v.ml0.ccsd_p_T_e[0] * 100, 2),
                        'pe_scf_ml1':  round((v.ml1.scf_e[0] - v.ml1.ccsd_p_T_e[0])/v.ml1.ccsd_p_T_e[0] * 100, 2),
                        'pe_scf':      round((v.ave.scf_e[0] - v.ave.ccsd_p_T_e[0])/v.ave.ccsd_p_T_e[0] * 100, 2),
                        'pe_mp2_ml0':  round((v.ml0.mp2_e[0] - v.ml0.ccsd_p_T_e[0])/v.ml0.ccsd_p_T_e[0] * 100, 2),
                        'pe_mp2_ml1':  round((v.ml1.mp2_e[0] - v.ml1.ccsd_p_T_e[0])/v.ml1.ccsd_p_T_e[0] * 100, 2),
                        'pe_mp2':      round((v.ave.mp2_e[0] - v.ave.ccsd_p_T_e[0])/v.ave.ccsd_p_T_e[0] * 100, 2),
                        'pe_ccsd_ml0': round((v.ml0.ccsd_e[0] - v.ml0.ccsd_p_T_e[0])/v.ml0.ccsd_p_T_e[0] * 100, 2),
                        'pe_ccsd_ml1': round((v.ml1.ccsd_e[0] - v.ml1.ccsd_p_T_e[0])/v.ml1.ccsd_p_T_e[0] * 100, 2),
                        'pe_ccsd':     round((v.ave.ccsd_e[0] - v.ave.ccsd_p_T_e[0])/v.ave.ccsd_p_T_e[0] * 100, 2),
                        'has_line':    has_line
                    })
                    ele_Z = ''
                    ele_symbol = '   '
                    has_line = '\cline{3-7}\n         '
                else:
                    if not self.has_header:
                        ele_Z = ''
                        ele_symbol = '   '

                    tmp_ct = ct.replace(ele.symbol + '@Q-', '')
                    if self.only_best and not SRCCDataProcessor.is_the_best(tmp_ct, ele.symbol):
                        continue

                    body_tex += body_q.format(**{
                        'ml0_scf':  round(v.ml0.scf_e[1]/4., 2),
                        'ml0_ccsd': round(v.ml0.ccsd_e[1]/4., 2),
                        'ml1_scf':  round(v.ml1.scf_e[1]/4., 2),
                        'ml1_ccsd': round(v.ml1.ccsd_e[1]/4., 2),
                        'ave_scf':  round(v.ave.scf_e[1]/4., 2),
                        'ave_ccsd': round(v.ave.ccsd_e[1]/4., 2),
                        'calc_type': SRCCDataProcessor.get_new_ct(tmp_ct),
                        'ele_Z': ele_Z,
                        'ele_symbol': ele_symbol,
                        'gs_term': self.gs_term,
                        'delta_scf':  round((v.ml1.scf_e[0] - v.ml0.scf_e[0])/4., 2),
                        'delta_ccsd': round((v.ml1.ccsd_e[0] - v.ml0.ccsd_e[0])/4., 2),
                        'pe_scf_ml0': round((v.ml0.scf_e[0] - v.ml0.ccsd_e[0])/v.ml0.ccsd_e[0] * 100, 2),
                        'pe_scf_ml1': round((v.ml1.scf_e[0] - v.ml1.ccsd_e[0])/v.ml1.ccsd_e[0] * 100, 2),
                        'pe_scf':     round((v.ave.scf_e[0] - v.ave.ccsd_e[0])/v.ave.ccsd_e[0] * 100, 2),
                        'has_line':    has_line
                    })
                    ele_Z = ''
                    ele_symbol = '   '


        elif self.res_ml0:
            for ct, v in sorted(self.all_res.items(), reverse=True):
                if not self.is_quadru:
                    if not self.has_header:
                        ele_Z = ''
                        ele_symbol = '   '

                    tmp_ct = ct.replace(ele.symbol + '@D-', '')
                    if self.only_best and not SRCCDataProcessor.is_the_best(tmp_ct, ele.symbol):
                        continue

                    body_tex += body_only_ml0.format(**{
                        'ml0_scf':    round(v.ml0.scf_e[0], 2),
                        'ml0_mp2':    round(v.ml0.mp2_e[0], 2),
                        'ml0_ccsd':   round(v.ml0.ccsd_e[0], 2),
                        'ml0_ccsd_t': round(v.ml0.ccsd_p_T_e[0], 2),
                        'ave_scf':    round(v.ave.scf_e[0], 2),
                        'ave_mp2':    round(v.ave.mp2_e[0], 2),
                        'ave_ccsd':   round(v.ave.ccsd_e[0], 2),
                        'ave_ccsd_t': round(v.ave.ccsd_p_T_e[0], 2),
                        'calc_type': SRCCDataProcessor.get_new_ct(tmp_ct),
                        'ele_Z' : ele_Z,
                        'ele_symbol': ele_symbol,
                        'gs_term': self.gs_term,
                        'pe_scf_ml0':  round((v.ml0.scf_e[0] - v.ml0.ccsd_p_T_e[0])/v.ml0.ccsd_p_T_e[0] * 100, 2),
                        'pe_mp2_ml0':  round((v.ml0.mp2_e[0] - v.ml0.ccsd_p_T_e[0])/v.ml0.ccsd_p_T_e[0] * 100, 2),
                        'pe_ccsd_ml0': round((v.ml0.ccsd_e[0] - v.ml0.ccsd_p_T_e[0])/v.ml0.ccsd_p_T_e[0] * 100, 2),
                        'has_line':    has_line
                    }
                                                     )
                    ele_Z = ''
                    ele_symbol = '   '
                    has_line = '\cline{3-7}\n         '
                else:
                    if not self.has_header:
                        ele_Z = ''
                        ele_symbol = '   '

                    tmp_ct = ct.replace(ele.symbol + '@Q-', '')
                    if self.only_best and not SRCCDataProcessor.is_the_best(tmp_ct, ele.symbol):
                        continue

                    body_tex += body_only_ml0_q.format(**{
                        'ml0_scf':  round(v.ml0.scf_e[1]/4., 2),
                        'ml0_ccsd': round(v.ml0.ccsd_e[1]/4., 2),
                        'ave_scf':  round(v.ave.scf_e[1]/4., 2),
                        'ave_ccsd': round(v.ave.ccsd_e[1]/4., 2),
                        'calc_type': SRCCDataProcessor.get_new_ct(tmp_ct),
                        'ele_Z' : ele_Z,
                        'ele_symbol': ele_symbol,
                        'gs_term': self.gs_term,
                        'pe_scf_ml0': round((v.ml0.scf_e[0] - v.ml0.ccsd_e[0])/v.ml0.ccsd_e[0] * 100, 2),
                        'has_line':    has_line,
                    })
                    ele_Z = ''
                    ele_symbol = '   '

        #template += footer.format(**{'DQ':dq})
        #print(template)
        # return template, footer.format(**{'DQ':dq, 'ele_type':ele.symbol})
        # return template, long_table_footer.format(**{'DQ':dq, 'ele_type':ele.symbol})
        footer_tex = long_table_footer.format(**{'DQ':dq, 'ele_type':ele.symbol})
        return header_tex, body_tex, footer_tex

    def get_gs(self):
        ele_sybmol = self.element.symbol

        if ele_sybmol in SRCCDataProcessor.best_dict['CC-SR']:
            gs_ct_nr = SRCCDataProcessor.best_dict['CC-SR'][ele_sybmol][0]
            gs_ct_so = SRCCDataProcessor.best_dict['CC-SR'][ele_sybmol][1]
            if self.is_quadru:
                key_nr = ele_sybmol + '@Q-' + gs_ct_nr
                key_so = ele_sybmol + '@Q-' + gs_ct_so
            else:
                key_nr = ele_sybmol + '@D-' + gs_ct_nr
                key_so = ele_sybmol + '@D-' + gs_ct_so
            gs_res_nr = self.all_res[key_nr]
            gs_res_so = self.all_res[key_so]
        else:
            gs_res_nr = None
            gs_res_so = None

        return gs_res_nr, gs_res_so


class NRCCDataProcessor(SRCCDataProcessor):

    def __init__(self, dirname, **kwargs):
        super(NRCCDataProcessor, self).__init__(dirname, **kwargs)

    def get_gs(self):
        round_nb = 2
        ele_sybmol = self.element.symbol
        gs_res_nr = Settings()

        if ele_sybmol in NRCCDataProcessor.best_dict['CC-NR']:
            gs_ct_nr = NRCCDataProcessor.best_dict['CC-NR'][ele_sybmol][0]
            if self.is_quadru:
                key_nr = ele_sybmol + '@Q-' + gs_ct_nr
                pos = 1
                scale = 4
            else:
                key_nr = ele_sybmol + '@D-' + gs_ct_nr
                pos = 0
                scale = 1
            gs_res_nr = self.all_res[key_nr]

        else:
            gs_alpha = np.nan
        return gs_res_nr


class SOCCDataProcessor(AbstractDataProcessor):
    def __init__(self, dirname, only_best=False):
        super(SOCCDataProcessor, self).__init__()

        self.dirname = dirname
        self.only_best = only_best

        self.raw_res = get_polarizability(self.dirname,
                                 calc_dir_patters = ['dyall', 'ANO-RCC', 'faegri'],
                                 verbos=False, deepth=0)
        self._extract_info()

    def _extract_info(self):
        self.all_res = Settings()
        ele_type = None
        self.is_quadru = False

        if len(self.raw_res) < 1:
            return

        for k, v in self.raw_res.items():
            if k in ['curr_dir', 'calc_dir'] and len(v) > 0:
                for kk, vv in v.items():
                    # kk is calc_type, e.g., Li@D-CC-2C-SR@ANO-RCC
                    # then, we collect data for different calc_type, i.e., kk
                    ele_type = kk.split('@')[0].strip()
                    if '@Q' in kk:
                        self.is_quadru = True
                    self.all_res[kk] = Settings(vv)


        self.element = Element(ele_type)

        if self.element.group in [1, 11]:
            self.gs_term = r'$^2S_{1/2}$'
        elif self.element.group in [2, 12, 18]:
            self.gs_term = r'$^1S_0$'
        elif self.element.group in [13]:
            self.gs_term = r'$^2P_{1/2}$'
        elif self.element.group in [14]:
            self.gs_term = r'$^3P_0$'
        else:
            raise RuntimeError('Unsupport group {0}'.format(self.element.group))

    def __str__(self):
        return self.get_tex()

    def get_tex(self):
        # return CCso_tex(self.raw_res)
        body = r"""          \cline{{3-7}}
          &      & {gs_term}    & {mj_scf:.2f}    & {pe_scf:.2f} & DHF     & {calc_type} \\
          &      &                & {mj_mp2:.2f}    & {pe_mp2:.2f} & MP2     &             \\
          &      &                & {mj_ccsd:.2f}    & {pe_ccsd:.2f} & CCSD    &             \\
          &      &                & {mj_ccsd_t:.2f}    &  & CCSD(T) &             \\
"""

        body_d_noT = r"""          \cline{{3-7}}
          &      & {gs_term}    & {mj_scf:.2f}    & {pe_scf:.2f} & DHF     & {calc_type} \\
          &             &               & {mj_mp2:.2f}   & {pe_mp2:.2f} &  MP2        &             \\
          &             &               & {mj_ccsd:.2f}   &  &  CCSD       &             \\
"""

        body_q = r"""          \cline{{3-7}}
          &      & {gs_term}    & {mj_scf:.2f}    & {pe_scf:.2f}  & DHF     & {calc_type} \\
          &      &                  & {mj_ccsd:.2f}    &  &  CCSD    &             \\
"""

        if self.is_quadru:
            pos = 1
            scale = 4
        else:
            pos = 0
            scale = 1

        body_lis = []
        for ct, v in sorted(self.all_res.items(), reverse=True):
            if not self.is_quadru:
                tmp_ct = ct.replace(self.element.symbol + '@D-', '')
                if self.only_best and not SOCCDataProcessor.is_the_best(tmp_ct, self.element.symbol):
                    continue

                if 'ccsd_p_T_e' not in v:
                    body_lis.append(body_d_noT.format(**{
                        'mj_scf': round(v.scf_e[pos], 2),
                        'mj_mp2': round(v.mp2_e[pos], 2),
                        'mj_ccsd': round(v.ccsd_e[pos], 2),
                        'calc_type': SOCCDataProcessor.get_new_ct(tmp_ct),
                        'gs_term': self.gs_term,
                        'pe_scf': round(
                            (v.scf_e[pos] - v.ccsd_e[pos]) / v.ccsd_e[
                                pos] * 100, 2),
                        'pe_mp2': round(
                            (v.mp2_e[pos] - v.ccsd_e[pos]) / v.ccsd_e[
                                pos] * 100, 2),
                    }))
                else:
                    body_lis.append(body.format(**{
                        'mj_scf': round(v.scf_e[pos], 2),
                        'mj_mp2': round(v.mp2_e[pos], 2),
                        'mj_ccsd': round(v.ccsd_e[pos], 2),
                        'mj_ccsd_t': round(v.ccsd_p_T_e[pos], 2),
                        'calc_type': SOCCDataProcessor.get_new_ct(tmp_ct),
                        'gs_term': self.gs_term,
                        'pe_scf': round(
                            (v.scf_e[pos] - v.ccsd_p_T_e[pos]) / v.ccsd_p_T_e[
                                pos] * 100, 2),
                        'pe_mp2': round(
                            (v.mp2_e[pos] - v.ccsd_p_T_e[pos]) / v.ccsd_p_T_e[
                                pos] * 100, 2),
                        'pe_ccsd': round(
                            (v.ccsd_e[pos] - v.ccsd_p_T_e[pos]) / v.ccsd_p_T_e[
                                pos] * 100, 2),
                    }))
            else:
                tmp_ct = ct.replace(self.element.symbol + '@Q-', '')
                if self.only_best and not SOCIDataProcessor.is_the_best(tmp_ct, self.element.symbol):
                    continue

                body_lis.append(body_q.format(**{
                    'mj_scf': round(v.scf_e[pos] / scale, 2),
                    'mj_ccsd': round(v.ccsd_e[pos] / scale, 2),
                    'calc_type': SOCCDataProcessor.get_new_ct(tmp_ct),
                    'gs_term': self.gs_term,
                    'pe_scf': round(
                        (v.scf_e[pos] - v.ccsd_e[pos]) / v.ccsd_e[pos] * 100,
                        2),
                }))

        body_tex = '\n'.join(body_lis)
        return body_tex

    def get_gs(self):
        ele_sybmol = self.element.symbol

        if ele_sybmol in SOCCDataProcessor.best_dict['CC-SO']:
            gs_ct = SOCCDataProcessor.best_dict['CC-SO'][ele_sybmol][0]
            if self.is_quadru:
                key = ele_sybmol + '@Q-' + gs_ct
            else:
                key = ele_sybmol + '@D-' + gs_ct
            gs_res = self.all_res[key]

        else:
            gs_res = Settings()
        return gs_res


class SOCIDataProcessor(AbstractDataProcessor):
    def __init__(self, dirname, only_best=False):
        super(SOCIDataProcessor, self).__init__()

        self.dirname = dirname
        self.only_best = only_best

        self.raw_res = get_polarizability(
            self.dirname, ['dyall', 'ANO-RCC', 'faegri'], 0, verbos=False)

        self._extract_info()

    def _extract_info(self):
        if len(self.raw_res) < 1:
            return

        s = Settings(self.raw_res).flatten()

        for k in s.keys():
            if len(k) == 3:
                ele = k[1].split('@')[0]
                break
        else:
            ele = ''
        self.element = Element(ele)

        calc_types = []
        for k in s.keys():
            if len(k) == 3:
                calc_types.append(k[1])
        calc_types = set(calc_types)

        # new_s =  { ('I@D-CI-4C-SO@dyall.acv4z@(core 17)[vir 215]', 'sym_3_root_1'): [1111.11, 2222.222], ... }
        tmp_all_res = Settings()
        for sk in s.keys():
            # here, s is a flatten result dict
            # sk is ('curr_dir', 'I@D-CI-4C-SO@dyall.acv4z@(core 17)[vir 215]', 'sym_3_root_1')
            if len(sk) == 3:
                tmp_all_res[sk[1:]] = s[sk]

        if '@Q' in ''.join(list(tmp_all_res.keys())[0]):
            self.is_quadru = True
        else:
            self.is_quadru = False
        self.all_res = tmp_all_res.unflatten()

    def __str__(self):
        return self.get_tex()

    def get_tex(self):
        # return CIso_tex(self.raw_res, only_best=self.only_best)
        if self.is_quadru:
            pos = 1
            dq_tag = 'Q'
            scale = 4
        else:
            pos = 0
            dq_tag = 'D'
            scale = 1

        body_lis = []
        for k, v in sorted(self.all_res.items(), reverse=True):
            # k is calc_type, e.g., I@D-CI-4C-SO@dyall.acv4z@(core 17)[vir 215]
            # if not k in tex_str:
            #     tex_str[k] = Settings()

            body = ''
            ct = k.replace(self.element.symbol + '@' + dq_tag + '-', '')
            if self.only_best and not SOCIDataProcessor.is_the_best(ct, self.element.symbol):
                continue

            ct = SOCIDataProcessor.get_new_ct(ct)

            if self.element.group in [1, 11]:
                # 2^S_{1/2}
                body = r"""          \cline{{3-7}}
                  &      & $^2S_{{1/2}}$    & {gs:.2f}    &  & MRCISD  & {calc_type}             \\
""".format(**{
                    'gs': round(v['sym_3_root_1'][pos] / scale, 2),
                    'calc_type': ct,
                })
            elif self.element.group in [13]:
                # B, Al, Ga, In, Tl, Nh
                # 2^P_{1/2}
                body = r"""          \cline{{3-7}}
          &      & $^2P_{{1/2}}$               & {gs:.2f}    &  & MRCISD  & {calc_type}            \\
          &      & $^2P_{{3/2}}$, $M_J=\pm1/2$ & {32_12:.2f}    &  &  &             \\
          &      & $^2P_{{3/2}}$, $M_J=\pm3/2$ & {32_32:.2f}    &  &  &             \\
          &      & $^2P_{{3/2}}$               & {32_ave:.2f}    &  &  &  $\bar{{\alpha}}_{{J=3/2}} - \alpha_{{J=1/2}} = {diff:.2f}$           \\
""".format(**{
                    'gs': round(v['sym_3_root_1'][pos] / scale, 2),
                    '32_12': round(v['sym_3_root_2'][pos] / scale, 2),
                    '32_32': round(v['sym_3_root_3'][pos] / scale, 2),
                    '32_ave': round((v['sym_3_root_2'][pos] + v['sym_3_root_3'][
                        pos]) / 2. / scale, 2),
                    'calc_type': ct,
                    'diff': round((v['sym_3_root_2'][pos] + v['sym_3_root_3'][
                        pos]) / 2. / scale - v['sym_3_root_1'][pos] / scale, 2),
                })
            elif self.element.group in [2, 12, 18]:
                # 1^S_0
                body = r"""          \cline{{3-7}}
          &      & $^1S_0$        & {gs:.2f}    &  & MRCISD  & {calc_type}    \\
""".format(**{
                    'gs': round(v['sym_1_root_1'][pos] / scale, 2),
                    'calc_type': ct,
                })
            elif self.element.group in [14]:
                # 3^P_0
                body = r"""          \cline{{3-7}}
          &      & $^3P_0$        & {gs:.2f}    &  & MRCISD  & {calc_type}    \\
""".format(**{
                    'gs': round(v['sym_1_root_1'][pos] / scale, 2),
                    'calc_type': ct,
                })
            elif self.element.group in [15]:
                # N, P, As, Sb, Bi, Mc
                # 4^S_{3/2}
                body = r"""          \cline{{3-7}}
          &      & $^4S_{{3/2}}$, $M_J=\pm1/2$ & {32_12:.2f}    &  & MRCISD  &  {calc_type}           \\
          &      & $^4S_{{3/2}}$, $M_J=\pm3/2$ & {32_32:.2f}    &  &  &             \\
          &      & $^4S_{{3/2}}$               & {32_ave:.2f}    &  &  &             \\
""".format(**{
                    '32_12': round(v['sym_3_root_1'][pos] / scale, 2),
                    '32_32': round(v['sym_3_root_2'][pos] / scale, 2),
                    '32_ave': round((v['sym_3_root_1'][pos] + v['sym_3_root_2'][
                        pos]) / 2. / scale, 2),
                    'calc_type': ct
                })
            elif self.element.group in [16]:
                # 3^P_2
                body = r"""          \cline{{3-7}}
          &      & $^3P_2$, $M_J=0$    & {2_0:.2f}    &  & MRCISD  &  {calc_type}           \\
          &      & $^3P_2$, $M_J=\pm1$ & {2_1:.2f}   &  &       &             \\
          &      & $^3P_2$, $M_J=\pm2$ & {2_2:.2f}   &  &       &             \\
          &      & $^3P_2$             & {2_ave:.2f}     &  &       &             \\
""".format(**{
                    '2_0': round(v['sym_1_root_1'][pos] / scale, 2),
                    '2_1': round((v['sym_1_root_2'][pos] + v['sym_1_root_3'][
                       pos]) / 2. / scale, 2),
                    '2_2': round((v['sym_2_root_1'][pos] + v['sym_2_root_2'][
                       pos]) / 2. / scale, 2),
                    '2_ave': round(
                        (v['sym_1_root_1'][pos] + v['sym_1_root_2'][pos] +
                         v['sym_1_root_3'][pos] + v['sym_2_root_1'][pos] +
                         v['sym_2_root_2'][pos]) / 5. / scale, 2),
                    'calc_type': ct,
                })

            elif self.element.group in [17]:
                # 2^P_{3/2}
                body = r"""          \cline{{3-7}}
          &      & $^2P_{{3/2}}$, $M_J=\pm1/2$ & {32_12:.2f}    &  & MRCISD  & {calc_type}            \\
          &      & $^2P_{{3/2}}$, $M_J=\pm3/2$ & {32_32:.2f}   &  &       &             \\
          &      & $^2P_{{3/2}}$               & {32_ave:.2f}      &  &       &             \\
          &      & $^2P_{{1/2}}$               & {12_12:.2f}      &  &       & $\bar{{\alpha}}_{{J=3/2}} - \alpha_{{J=1/2}} = {diff:.2f}$             \\
""".format(**{
                    '32_12': round(v['sym_3_root_2'][pos] / scale, 2),
                    '32_32': round(v['sym_3_root_1'][pos] / scale, 2),
                    '32_ave': round((v['sym_3_root_1'][pos] + v['sym_3_root_2'][
                        pos]) / 2. / scale, 2),
                    '12_12': round(v['sym_3_root_3'][pos] / scale, 2),
                    'calc_type': ct,
                    'diff': round((v['sym_3_root_1'][pos] + v['sym_3_root_2'][
                        pos]) / 2. / scale - v['sym_3_root_3'][pos] / scale, 2),
                })
            else:
                raise RuntimeError('Not support group {0}'.format(self.element.group))
            body_lis.append(body)

        body_tex = '\n'.join(body_lis)
        return body_tex

    def get_gs(self):
        round_nb = 2
        ele_sybmol = self.element.symbol

        if ele_sybmol in SOCIDataProcessor.best_dict['CI-SO']:
            gs_ct = SOCIDataProcessor.best_dict['CI-SO'][ele_sybmol][0]
            if self.is_quadru:
                key = ele_sybmol + '@Q-' + gs_ct
                pos = 1
                scale = 4
            else:
                key = ele_sybmol + '@D-' + gs_ct
                pos = 0
                scale = 1
            gs_res = self.all_res[key]

            if self.element.group in [1, 11]:
                gs_alpha = round(gs_res['sym_3_root_1'][pos] / scale, round_nb)
            elif self.element.group in [13]:
                gs_alpha = round(gs_res['sym_3_root_1'][pos] / scale, round_nb)
            elif self.element.group in [2, 12, 18]:
                gs_alpha = round(gs_res['sym_1_root_1'][pos] / scale, round_nb)
            elif self.element.group in [14]:
                gs_alpha = round(gs_res['sym_1_root_1'][pos] / scale, round_nb)
            elif self.element.group in [15]:
                gs_alpha = round((gs_res['sym_3_root_1'][pos] +
                                  gs_res['sym_3_root_2'][pos]) / 2. / scale, round_nb)
            elif self.element.group in [16]:
                gs_alpha = round(
                    (gs_res['sym_1_root_1'][pos] + gs_res['sym_1_root_2'][pos] +
                     gs_res['sym_1_root_3'][pos] + gs_res['sym_2_root_1'][pos] +
                     gs_res['sym_2_root_2'][pos]) / 5. / scale, round_nb)
            elif self.element.group in [17]:
                gs_alpha = round((gs_res['sym_3_root_1'][pos] +
                                  gs_res['sym_3_root_2'][pos]) / 2. / scale, round_nb)
            else:
                raise RuntimeError('Not support group {0}'.format(self.element.group))
        else:
            gs_alpha = np.nan
            gs_res = Settings()

        gs_res['gs'] = gs_alpha
        return gs_res


class Ploter:
    def __init__(self, data, save_path='.', fname='group_{}.pdf'):
        self.data = data
        self.group = self.data['group']
        self.atom_list = data['atom_list']
        self.save_path = os.path.abspath(save_path)
        self.fname = fname

        self.fig = plt.figure(figsize=(10.24, 8))
        self.ax1 = plt.subplot2grid((2, 2), (0, 0))
        self.ax2 = plt.subplot2grid((2, 2), (0, 1))
        self.ax3 = plt.subplot2grid((2, 2), (1, 0))
        self.ax4 = plt.subplot2grid((2, 2), (1, 1))
        self.axs = [self.ax1, self.ax2, self.ax3, self.ax4]


    def plot(self):
        self.set_tick()
        self.set_locator_size()

        self.plot_polarizability()
        self.plot_correlation_contribution()
        self.plot_tot_rel()
        self.plot_resolved_rel()

        self.set_xylim()

        # global plot setup
        self.fig.tight_layout()  # otherwise the right y-label is slightly clipped
        filename = os.path.join(self.save_path, self.fname.format(self.group))
        plt.savefig(filename, bbox_inches='tight')
        plt.show()

    def plot_polarizability(self):
        # -----------------------------------
        # NR-, SR-, DC-polarizabilities
        # -----------------------------------
        if self.group in [13, 14, 16, 17]:
            # plot for Ml0 and Ml1 of alpha_0
            alpha_0_ml0, alpha_0_ml1, alpha_1_ml0, alpha_1_ml1 = \
                self.get_data_for_ml(ml0_callback=lambda x: x.ml0.ccsd_p_T_e[0],
                                ml1_callback=lambda x: x.ml1.ccsd_p_T_e[0])

            self.ax1.plot(self.atom_list, alpha_0_ml0, '--', color=red, alpha=0.8,
                     label=r'$\alpha^{(0)}_{M_L=0}$', markersize=7)
            self.ax1.plot(self.atom_list, alpha_0_ml1, '-.', color=red, alpha=0.8,
                     label=r'$\alpha^{(0)}_{M_L=\pm1}$', markersize=7)
            self.ax1.plot(self.atom_list, self.data['alpha_0'], '-D', color=red,
                     label=r'$\alpha^{(0)}$', markersize=7)

            self.ax1.plot(self.atom_list, alpha_1_ml0, '--', color=blue, alpha=0.8,
                     label=r'$\alpha^{(1)}_{M_L=0}$', markersize=7)
            self.ax1.plot(self.atom_list, alpha_1_ml1, '-.', color=blue, alpha=0.8,
                     label=r'$\alpha^{(1)}_{M_L=\pm1}$', markersize=7)
            self.ax1.plot(self.atom_list, self.data['alpha_1'], '-s', color=blue,
                     label=r'$\alpha^{(1)}$', markersize=7)

            self.ax1.plot(self.atom_list, self.data['alpha_2'], '-o', color='k',
                     label=r'$\alpha^{(2)}$', markersize=7,
                     markerfacecolor='w', markeredgewidth=1,
                     markeredgecolor='k')

        else:
            self.ax1.plot(self.atom_list, self.data['alpha_0'], '-D', color=red,
                     label=r'$\alpha^{(0)}$', markersize=7)
            self.ax1.plot(self.atom_list, self.data['alpha_1'], '-s', color=blue,
                     label=r'$\alpha^{(1)}$', markersize=7)
            # ax1.plot(atom_list, data['alpha_2'], 'x', color='k', label=r'$\alpha^{(2)}$', markersize=10)
            self.ax1.plot(self.atom_list, self.data['alpha_2'], '--o', color='k',
                     label=r'$\alpha^{(2)}$', markersize=7,
                     markerfacecolor='w', markeredgewidth=1,
                     markeredgecolor='k')

        self.set_xylabel(self.ax1)
        self.set_legend(self.ax1, 1)

    def plot_correlation_contribution(self):
        # -----------------------------------
        # Electron-correlation contribution
        # -----------------------------------
        self.ax2.plot(self.atom_list, self.data['corr_0'], '-D', color=red,
                 label=r'$\alpha_{corr}^{(0)}$', markersize=7)
        self.ax2.plot(self.atom_list, self.data['corr_1'], '-s', color=blue,
                 label=r'$\alpha_{corr}^{(1)}$', markersize=7)
        if not self.group in [15, 16, 17]:
            self.ax2.plot(self.atom_list, self.data['corr_2'], '--o', color='k',
                     label=r'$\alpha_{corr}^{(2)}$', markersize=7,
                     markerfacecolor='w', markeredgewidth=1,
                     markeredgecolor='k')

        self.set_xylabel(self.ax2)
        self.set_legend(self.ax2, 2)

    def plot_tot_rel(self):
        # ---------------------------
        # Relativistic contribution
        # ---------------------------
        atom_Z_list = [Element(e).Z for e in self.atom_list]
        self.ax3.plot(atom_Z_list, self.data['delta_tot_1'], '-D', color=red,
                 label=r'$\Delta \alpha_{tot}^{(1)}$', markersize=7, )
        self.ax3.plot(atom_Z_list, self.data['delta_tot_so'], '-s', color=blue,
                 label=r'$\Delta \alpha_{tot}^{(SO)}$', markersize=7, )
        self.ax3.plot(atom_Z_list, np.asarray(self.data['delta_tot_1']) +
                 np.asarray(self.data['delta_tot_so']), '--o', color='k',
                 label=r'$\Delta \alpha_{tot}^{(2)}$', markersize=7,
                 markerfacecolor='w', markeredgewidth=1, markeredgecolor='k',
                 )

        self.set_xylabel(self.ax3, use_Z=True)
        self.set_legend(self.ax3, 3)

    def plot_resolved_rel(self):
        # ---------------------------------------------
        # Component-resolved relativistic contribution
        # ---------------------------------------------
        self.ax4.plot(self.atom_list, self.data['delta_orb_1'], '-D', color=red,
                 label=r'$\Delta \alpha_{orb}^{(1)}$', markersize=7)
        self.ax4.plot(self.atom_list, self.data['delta_corr_1'], '-s', color=blue,
                 label=r'$\Delta \alpha_{corr}^{(1)}$', markersize=7)
        if not self.group in [15, 16, 17]:
            self.ax4.plot(self.atom_list, self.data['delta_orb_so'], '--D', color=red,
                     label=r'$\Delta \alpha_{orb}^{(SO)}$', markersize=7,
                     markerfacecolor='w', markeredgewidth=1,
                     )
            self.ax4.plot(self.atom_list, self.data['delta_corr_so'], '--s', color=blue,
                     label=r'$\Delta \alpha_{corr}^{(SO)}$', markersize=7,
                     markerfacecolor='w', markeredgewidth=1,
                     )

        self.set_xylabel(self.ax4)
        self.set_legend(self.ax4, 4)

    def get_data_for_ml(self, ml0_callback, ml1_callback,
                        nr_key='gs_nr', sr_key='gs_for_nr'):
        # data for Ml0 and Ml1 of alpha_0
        nr_0_ml0 = []
        nr_0_ml1 = []
        for i in self.data[nr_key]:
            # nr_0_ml0.append(i.ml0.ccsd_p_T_e[0] - i.ml0.scf_e[0])
            nr_0_ml0.append(ml0_callback(i))
            nr_0_ml1.append(ml1_callback(i))

        # data for Ml0 and Ml1 of alpha_1
        sr_1_ml0 = []
        sr_1_ml1 = []
        for i in self.data[sr_key]:
            sr_1_ml0.append(ml0_callback(i))
            sr_1_ml1.append(ml1_callback(i))
        return nr_0_ml0, nr_0_ml1, sr_1_ml0, sr_1_ml1

    def set_xylim(self):
        # Set the axis limits
        for i, ax in enumerate(self.axs):
            ylim = ax.get_ylim()
            xlim = ax.get_xlim()

            ylim_l = ylim[1] - ylim[0]
            xlim_l = xlim[1] - xlim[0]
            scale = 0.05
            ylim = (ylim[0] - scale*ylim_l, ylim[1]+scale*ylim_l)
            xlim = (xlim[0] - scale*xlim_l, xlim[1]+scale*xlim_l)

            ylim_l = ylim[1] - ylim[0]
            xlim_l = xlim[1] - xlim[0]

            if self.group in [11]:
                x_center = xlim_l*0.8 + xlim[0]
                y_center = ylim_l*0.5 + ylim[0]
            else:
                x_center = xlim_l*0.1 + xlim[0]
                y_center = ylim_l*0.5 + ylim[0]

            ax.set_ylim(*ylim)
            ax.set_xlim(*xlim)
            assert i+1 in [1,2,3,4]
            tag = '(' + 'abcd'[i] + ')'
            ax.text(x_center, y_center, tag, horizontalalignment='center',
                    verticalalignment='center',)

    def set_tick(self):
        for ax in self.axs:
            ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')
            ax.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
            ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in',right='on')
            ax.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')

    def set_xylabel(self, ax, use_Z=False):
        if use_Z:
            ax.set_xlabel(r'Atomic number $Z$')
        else:
            ax.set_xlabel('Group-{0} elements'.format(self.group))
        ax.set_ylabel(r'$\alpha$ [a.u.]')

    def set_locator_size(self):
        ax_nb = len(self.axs)
        if self.group in [1]:
            major_size = [50, 50, 25, 100]
            minor_size = np.asarray(major_size)/5
        elif self.group in [2]:
            major_size = [50, 10, 25, 20]
            minor_size = np.asarray(major_size)/5
        elif self.group in [11]:
            major_size = [5, 10, 5, 20]
            minor_size = np.asarray(major_size)/5
        elif self.group in [12]:
            major_size = [10, 5, 10, 20]
            minor_size = np.asarray(major_size)/5
        elif self.group in [13]:
            major_size = [25, 5, 20, 20]
            minor_size = np.asarray(major_size)/5
        elif self.group in [14]:
            major_size = [25, 2, 10, 10]
            minor_size = np.asarray(major_size) / 5
        elif self.group in [15]:
            major_size = [15, 1, 4, 1]
            minor_size = np.asarray(major_size) / 5
        elif self.group in [16]:
            major_size = [15, 0.5, 4, 1]
            minor_size = np.asarray(major_size) / 5
        elif self.group in [17]:
            major_size = [15, 0.5, 4, 0.5]
            minor_size = np.asarray(major_size) / 5
        elif self.group in [18]:
            major_size = [15, 0.5, 4, 4]
            minor_size = np.asarray(major_size) / 5
        else:
            major_size = [50]*ax_nb
            minor_size = [10]*ax_nb

        for i, ax in enumerate([self.ax1, self.ax2, self.ax3, self.ax4]):
            ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(major_size[i]))
            ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(minor_size[i]))


    def set_legend(self, ax, nb):
        if nb == 2:
            if self.group in [11, 13, 14, 16, 17]:
                ax.legend(loc=0, ncol=2, borderpad=0.4, labelspacing=0.1, handlelength=2,
                          handletextpad=0.2, columnspacing=0.5, frameon=False)
                # ax2.legend(loc=0, ncol=3, borderpad=0.4, labelspacing=0.1,
                #            handlelength=1,
                #            handletextpad=0.2, columnspacing=0.5, frameon=False)
            else:
                ax.legend(frameon=False)
        elif nb == 1:
            if self.group in [13, 14, 16, 17]:
                ax.legend(loc=0, ncol=3, borderpad=0.4, labelspacing=0.1,
                          handlelength=1,
                          handletextpad=0.2, columnspacing=0.5, frameon=False)
            elif self.group in [11]:
                ax.legend(loc=0, ncol=2, borderpad=0.4, labelspacing=0.1, handlelength=2,
                          handletextpad=0.2, columnspacing=0.5, frameon=False)
            else:
                ax.legend(frameon=False)
        elif nb == 3:
            ax.legend(frameon=False)
        elif nb == 4:
            if self.group in [1, 2, 11, 12, 13,  14, 18]:
                ax.legend(loc=0, ncol=2, borderpad=0.4, labelspacing=0.1,
                          handlelength=2,
                          handletextpad=0.2, columnspacing=0.5, frameon=False)

            else:
                ax.legend(frameon=False)
        else:
            pass


class DataProcessor:
    default_fname = 'data_group_{}.pkl'

    def __init__(self, scratch_dir='.', fname=None, data_root='.'):
        self.scratch_dir = os.path.abspath(scratch_dir)
        self.fname = fname or self.default_fname
        self.data_root=data_root

    def dump_data(self, group=None, refresh=False, all=False, **kwargs):
        if all or group is None:
            group_lis = [1, 2, 11, 12, 13, 14, 15, 16, 17, 18]
        else:
            group_lis = [group]

        for gp in group_lis :
            _fname = self.fname.format(gp)
            fname = os.path.join(self.scratch_dir, _fname)
            if os.path.exists(fname) and not refresh:
                continue
            else:
                data = self.relativistic_analysis(gp, **kwargs)
                with open(fname, 'wb') as f:
                    pickle.dump(data, f)

    def load_data(self, group, refresh=False):
        if refresh:
            self.dump_data(group, refresh)

        filename = os.path.join(self.scratch_dir, self.fname.format(group))
        if not os.path.exists(filename):
            self.dump_data(group)

        with open(filename, 'rb') as f:
            data = pickle.load(f)
        return data

    def relativistic_analysis(self, group, is_quadru=False):
        if group == 1:
            atom_list = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']
        elif group == 2:
            atom_list = ['Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra']
        elif group  == 11:
            atom_list = ['Cu', 'Ag', 'Au']
        elif group == 12:
            atom_list = ['Zn', 'Cd', 'Hg', 'Cn']
        elif group == 13:
            atom_list = ['B', 'Al', 'Ga', 'In', 'Tl', 'Nh']
        elif group == 14:
            atom_list = ['C', 'Si', 'Ge', 'Sn', 'Pb', 'Fl']
        elif group == 15:
            atom_list = ['N', 'P', 'As', 'Sb', 'Bi', 'Mc']
        elif group == 16:
            atom_list = ['O', 'S', 'Se', 'Te', 'Po', 'Lv']
        elif group == 17:
            atom_list = ['F', 'Cl', 'Br', 'I', 'At', 'Ts']
        elif group == 18:
            atom_list = ['He', 'Ne', 'Ar','Kr', 'Xe', 'Rn', 'Og']
        else:
            raise TypeError('group {} should be integer [1-2, 11-18]!'.format(group))

        calc_path = Path(self.data_root).resolve()

        alpha_0_lis = []
        alpha_1_lis = []
        alpha_2_lis = []
        corr_0_lis = []
        corr_1_lis = []
        corr_2_lis = []
        delta_corr_1_lis = []
        delta_corr_so_lis = []
        delta_orb_1_lis = []
        delta_orb_so_lis = []
        delta_tot_1_lis = []
        delta_tot_so_lis = []

        gs_nr_lis = []
        gs_for_nr_lis = []
        gs_for_so_lis = []
        gs_cc_so_lis = []
        gs_ci_so_lis = []

        so_ci_processor_lis = []

        with cd(calc_path):
            for e in atom_list:
                ele = Element(e)
                # CC-SR
                if is_quadru:
                    ccnr = ele.symbol + '_q_nr'
                    ccso = ele.symbol + '_q_so'
                    ccsr = ele.symbol + '_q_theta'
                    ciso = ele.symbol + '_q_mrci'
                else:
                    ccnr = ele.symbol + '_nr'
                    ccso = ele.symbol + '_so'
                    ccsr = ele.symbol
                    ciso = ele.symbol + '_mrci'

                if os.path.exists(ccnr) and os.path.isdir(ccnr):
                    with cd(ccnr):
                        nrcc_processor = NRCCDataProcessor(Path('.').resolve(), has_header=True, only_best=False)
                        gs_nr = nrcc_processor.get_gs()
                else:
                    gs_nr = None
                gs_nr_lis.append(gs_nr)

                if os.path.exists(ele.symbol) and os.path.isdir(ele.symbol):
                    with cd(ccsr):
                        srcc_processor = SRCCDataProcessor(Path('.').resolve(), only_best=False)
                        gs_for_nr, gs_for_so = srcc_processor.get_gs()
                else:
                    gs_for_nr, gs_for_so = None, None
                gs_for_nr_lis.append(gs_for_nr)
                gs_for_so_lis.append(gs_for_so)

                if os.path.exists(ccso) and os.path.isdir(ccso):
                    with cd(ccso):
                        if not os.path.exists('cannot'):
                            socc_processor = SOCCDataProcessor(Path('.').resolve(),
                                                               only_best=False)
                            so_ci_processor_lis.append(soci_processor)
                            gs_so = socc_processor.get_gs()
                        else:
                            so_ci_processor_lis.append(None)
                            gs_so = None
                else:
                    gs_so = None
                gs_cc_so_lis.append(gs_so)

                if os.path.exists(ciso) and os.path.isdir(ciso):
                    with cd(ciso):
                        soci_processor = SOCIDataProcessor(Path('.').resolve(),
                                                           only_best=False)
                        gs_ci_so = soci_processor.get_gs()
                else:
                    gs_ci_so = None
                gs_ci_so_lis.append(gs_ci_so)

                assert gs_nr is not None
                assert gs_for_nr is not None
                assert gs_for_so is not None
                assert gs_ci_so is not None

                print('Res of {0}:'.format(ele.symbol))
                alpha_corr_0 = gs_nr.ave.ccsd_p_T_e[0] - gs_nr.ave.scf_e[0]
                alpha_corr_1 = gs_for_nr.ave.ccsd_p_T_e[0] - gs_for_nr.ave.scf_e[0]
                print(r'$\alpha_{{corr}}^{{(0)}} = {0:.2f} $'.format(alpha_corr_0))
                print(r'$\alpha_{{corr}}^{{(1)}} = {0:.2f} $'.format(alpha_corr_1))
                alpha_0_lis.append(gs_nr.ave.ccsd_p_T_e[0])
                alpha_1_lis.append(gs_for_nr.ave.ccsd_p_T_e[0])

                if gs_so is None:
                    alpha_corr_2 = np.nan
                    alpha_2_lis.append(gs_ci_so.gs)
                else:
                    alpha_corr_2 = gs_so.ccsd_p_T_e[0] - gs_so.scf_e[0]
                    alpha_2_lis.append(gs_so.ccsd_p_T_e[0])
                print(r'$\alpha_{{corr}}^{{(2)}} = {0:.2f} $'.format(alpha_corr_2))

                corr_0_lis.append(alpha_corr_0)
                corr_1_lis.append(alpha_corr_1)
                corr_2_lis.append(alpha_corr_2)

                delta_corr_1 = alpha_corr_1 - alpha_corr_0
                delta_corr_so = alpha_corr_2 - alpha_corr_1
                print(r'$\Delta \alpha_{{corr}}^{{(1)}} = {0:.2f}$'.format(delta_corr_1))
                print(r'$\Delta \alpha_{{corr}}^{{SO}} = {0:.2f}$'.format(delta_corr_so))
                delta_corr_1_lis.append(delta_corr_1)
                delta_corr_so_lis.append(delta_corr_so)

                # delta alpha by SR orb
                delta_orb_1 = gs_for_nr.ave.scf_e[0] - gs_nr.ave.scf_e[0]
                print(r'$\Delta \alpha_{{orb}}^{{(1)}} = {0:.2f} $'.format(delta_orb_1))

                if gs_so is None:
                    delta_orb_so = np.nan
                else:
                    delta_orb_so = gs_so.scf_e[0] - gs_for_so.ave.scf_e[0]
                print(r'$\Delta \alpha_{{orb}}^{{(SO)}} = {0:.2f} $'.format(delta_orb_so))
                delta_orb_1_lis.append(delta_orb_1)
                delta_orb_so_lis.append(delta_orb_so)


                # delta alpha by tot rel
                delta_tot_1 = gs_for_nr.ave.ccsd_p_T_e[0] - gs_nr.ave.ccsd_p_T_e[0]
                print(r'$\Delta \alpha_{{tot}}^{{(1)}} = {0:.2f} $'.format(delta_tot_1))

                if gs_so is None:
                    delta_tot_so = gs_ci_so.gs - gs_for_so.ave.ccsd_p_T_e[0]
                else:
                    delta_tot_so = gs_so.ccsd_p_T_e[0] - gs_for_so.ave.ccsd_p_T_e[0]
                print(delta_tot_so)
                print(r'$\Delta \alpha_{{tot}}^{{(SO)}} = {0:.2f} $'.format(delta_tot_so))
                print()
                delta_tot_1_lis.append(delta_tot_1)
                delta_tot_so_lis.append(delta_tot_so)

        delta_corr_2_lis = list(np.asarray(delta_corr_1_lis) +
                                np.asarray(delta_corr_so_lis))
        delta_orb_2_lis = list(np.asarray(delta_orb_1_lis) +
                               np.asarray(delta_orb_so_lis))
        delta_tot_2_lis = list(np.asarray(delta_tot_1_lis) +
                               np.asarray(delta_tot_so_lis))

        data = {
            'group': group,
            'atom_list': atom_list,
            'alpha_0': alpha_0_lis,
            'alpha_1': alpha_1_lis,
            'alpha_2': alpha_2_lis,
            'corr_0': corr_0_lis,
            'corr_1': corr_1_lis,
            'corr_2': corr_2_lis,
            'delta_corr_1': delta_corr_1_lis,
            'delta_corr_so': delta_corr_so_lis,
            'delta_corr_2': delta_corr_2_lis,
            'delta_orb_1': delta_orb_1_lis,
            'delta_orb_so': delta_orb_so_lis,
            'delta_orb_2': delta_orb_2_lis,
            'delta_tot_1': delta_tot_1_lis,
            'delta_tot_so': delta_tot_so_lis,
            'delta_tot_2': delta_tot_2_lis,
            'gs_nr': gs_nr_lis,
            'gs_for_nr': gs_for_nr_lis,
            'gs_for_so': gs_for_so_lis,
            'gs_cc_so': gs_cc_so_lis,
            'gs_ci_so': gs_ci_so_lis,
        }

        return data


class Tabulator:

    def __init__(self, data):
        self.data = data
        self.atom_list = self.data['atom_list']
        self.group = self.data['group']
        self._table = None

    @property
    def alpha_NR(self):
        return self.data['alpha_0']

    @property
    def alpha_SR(self):
        return self.data['alpha_1']

    @property
    def alpha_DC(self):
        return self.data['alpha_2']

    @property
    def corr_NR(self):
        return self.data['corr_0']

    @property
    def corr_SR(self):
        return self.data['corr_1']

    @property
    def corr_DC(self):
        if self.group in [1, 2, 11, 12, 13, 14, 18]:
            return self.data['corr_2']
        else:
            return [np.nan] * len(self.atom_list)

    def get_data_for_ml(self, ml0_callback, ml1_callback,
                        nr_key='gs_nr', sr_key='gs_for_nr'):
        # data for Ml0 and Ml1 of alpha_0
        nr_0_ml0 = []
        nr_0_ml1 = []
        for i in self.data[nr_key]:
            # nr_0_ml0.append(i.ml0.ccsd_p_T_e[0] - i.ml0.scf_e[0])
            nr_0_ml0.append(ml0_callback(i))
            nr_0_ml1.append(ml1_callback(i))

        # data for Ml0 and Ml1 of alpha_1
        sr_1_ml0 = []
        sr_1_ml1 = []
        for i in self.data[sr_key]:
            sr_1_ml0.append(ml0_callback(i))
            sr_1_ml1.append(ml1_callback(i))
        return nr_0_ml0, nr_0_ml1, sr_1_ml0, sr_1_ml1

    @staticmethod
    def get_sym_root(group):
        if group in [1, 2, 11, 12, 18]:
            syms = (3,) if group % 2 else (1,)
            nb_roots = (1,)
            terms = [Tabulator.get_gs_term(group)]
        elif group in [13]:
            syms = (3,)
            nb_roots = (3,)
            terms = [Tabulator.get_gs_term(group),
                     '^2P_{3/2}, M_J=1/2', '^2P_{3/2}, M_J=3/2']
        elif group in [14]:
            syms = (1,)
            nb_roots = (2,)
            terms = [Tabulator.get_gs_term(group)]
        elif group in [15]:
            syms = (3,)
            nb_roots = (2,)
            terms = [Tabulator.get_gs_term(group) + ', M_J=1/2',
                     Tabulator.get_gs_term(group) + ', M_J=3/2']
        elif group in [16]:
            syms = (1, 2)
            # nb_roots = (3, 2) # this is because \pm 1 state are degenerated
            nb_roots = (2, 1)
            terms = [Tabulator.get_gs_term(group) + ', M_J=0' ,
                     Tabulator.get_gs_term(group) + ', M_J=1',
                     Tabulator.get_gs_term(group) + ', M_J=2']
        elif group in [17]:
            syms = (3,)
            nb_roots = (3,)
            terms = [Tabulator.get_gs_term(group) + ', M_J=3/2',
                     Tabulator.get_gs_term(group) + ', M_J=1/2',
                     '^2P_{1/2}, M_J=1/2']
        else:
            raise RuntimeError(
                'Group-{} elements does not support!'.format(group))

        sym_root_lis = []
        for sym, nb_root in zip(syms, nb_roots):
            sym_root_lis += ['sym_{}_root_{}'.format(s, r)
                             for s, r in
                             zip([sym] * nb_root, range(1, nb_root + 1))]
        return sym_root_lis, terms

    @staticmethod
    def get_gs_term(group, scheme='jj'):
        if group in [1, 11]:
            t = '^1S_{1/2}'
            t_LS = '^1S'
        elif group in [2, 12, 18]:
            t = '^1S_0'
            t_LS = '^1S'
        elif group in [13]:
            t = '^2P_{1/2}'
            t_LS = '^2P'
        elif group in [14]:
            t = '^3P_0'
            t_LS = '^3P'
        elif group in [15]:
            t = '^4S_{3/2}'
            t_LS = '^4S'
        elif group in [16]:
            t = '^3P_2'
            t_LS = '^3P'
        elif group in [17]:
            t = '^2P_{3/2}'
            t_LS = '^2P'
        else:
            raise NotImplementedError('Not support for other group {0}'.format(group))

        if scheme == 'jj':
            return t
        else:
            return t_LS


    @property
    def table(self):
        if self._table:
            return self._table

        new_data = OrderedDict()

        group = self.data['group']
        # NR and SR
        nr_sr_new_keys = ['NR', 'SR-NR', 'SR-DC']
        ave_tag = Tabulator.get_gs_term(group, scheme='LS')
        if group in [1, 2, 11, 12, 15, 18]:
            ml_component = ['ml0', 'ave']
            new_ml_component = ['M_L=0', ave_tag]
        else:
            ml_component = ['ml0', 'ml1', 'ave']
            new_ml_component = ['M_L=0', 'M_L=1', ave_tag]
        methods = ['scf', 'mp2', 'ccsd', 'ccsd_p_T']
        new_methods = ['DHF', 'MP2', 'CCSD', 'CCSD(T)']

        # ----------------------
        # copy data from old API
        # ----------------------
        nr_sr_old_keys = ['gs_nr', 'gs_for_nr', 'gs_for_so']
        for old_key, new_key in zip(nr_sr_old_keys, nr_sr_new_keys):
            for entity in data[old_key]:
                for kk, new_kk in zip(ml_component, new_ml_component):
                    for k, new_k in zip(methods, new_methods):
                        index = r'{}@{}@{}'.format(new_key, new_kk, new_k)

                        if index not in new_data:
                            new_data[index] = []
                        new_data[index].append(entity[kk][k+'_e'][0])
                # for kk in ml_component:
                #     for k in methods:
                #         index = '{}_{}_{}'.format(new_key, kk, k)
                #         new_data[index].append(entity[kk][k + '_e'][0])

        # copy DC results
        gs_term = Tabulator.get_gs_term(group)
        if group in [1, 2, 11, 12, 13, 14, 18]:
            # get gs term
            for entity in data['gs_cc_so']:
                for k, new_k in zip(methods, new_methods):
                    index = 'DC@{}@{}'.format(gs_term, new_k)

                    if index not in new_data:
                        new_data[index] = []
                    if entity is None:
                        new_data[index].append(np.nan)
                    else:
                        new_data[index].append(entity[k + '_e'][0])
        else:
            for entity in data['gs_ci_so']:
                index = 'DC@{}@MRCISD'.format(gs_term)

                if index not in new_data:
                    new_data[index] = []
                new_data[index].append(entity['gs'])

        # copy CI results
        for entity in data['gs_ci_so']:
            sym_root_lis, state_terms = Tabulator.get_sym_root(group)
            for sym_root, st in zip(sym_root_lis, state_terms):
                index = 'DC@{}@MRCISD'.format(st)
                if index not in new_data:
                    new_data[index] = []
                new_data[index].append(entity[sym_root][0])

        # copy CI excitation state res
        group_13_32_ave_index = r'DC@^2P_{3/2}@MRCISD'
        group_13_12_gs_index = r'DC@^2P_{1/2}@MRCISD'
        group_15_32_32_index = r'DC@^4S_{3/2}, M_J=3/2@MRCISD'
        group_15_32_12_index = r'DC@^4S_{3/2}, M_J=1/2@MRCISD'
        group_17_32_gs_index = r'DC@^2P_{3/2}@MRCISD'
        group_17_12_12_index = r'DC@^2P_{1/2}, M_J=1/2@MRCISD'
        if group in [13]:
            for idx in [group_13_32_ave_index, ]:
                if idx not in new_data:
                    new_data[idx] = []

            state_32_12 = [ i['sym_3_root_2'][0] for i in data['gs_ci_so']]
            state_32_32 = [ i['sym_3_root_3'][0] for i in data['gs_ci_so']]
            new_data[group_13_32_ave_index] = np.average([state_32_12, state_32_32], axis=0)
        # elif group in [15]:
        #     for idx in [group_15_32_ave_index]:
        #         if idx not in new_data:
        #             new_data[idx] = []
        #
        #     state_32_12 = [ i['sym_3_root_1'][0] for i in data['gs_ci_so']]
        #     state_32_32 = [ i['sym_3_root_2'][0] for i in data['gs_ci_so']]
        #     new_data[group_15_32_ave_index] = np.average([state_32_12, state_32_32], axis=0)



        # derivate results
        def get_tex_key(key):
            info = key.split('_')
            if len(info) == 2:
                alpha, order = info[:]
                if alpha == 'alpha':
                    new_key = r'\alpha^{{({})}}'.format(order)
                else:
                    new_key = r'\alpha_{{{}}}^{{({})}}'.format(alpha, order)
            elif len(info) == 3:
                d, comp, order = info[:]
                new_key = r'\Delta \alpha_{{{0}}}^{{(\text{{{1}}})}}'.format(comp, order)
            else:
                raise RuntimeError('Unknown key {0}'.format(key))
            return new_key

        calc_section_key = 'Calc.'
        # obtain Ml0 and Ml1 res
        if group in [13, 14, 16, 17]:
            selected_state = ave_tag
            for state in ['M_L=0', 'M_L=1']:
                nr_dhf = 'NR@{}@DHF'.format(state)
                nr_ccsd_t = 'NR@{}@CCSD(T)'.format(state)
                sr_nr_dhf = 'SR-NR@{}@DHF'.format(state)
                sr_nr_ccsd_t = 'SR-NR@{}@CCSD(T)'.format(state)

                tmp_d = OrderedDict()
                tmp_d['alpha_0'] = new_data[nr_ccsd_t]
                tmp_d['alpha_1'] = new_data[sr_nr_ccsd_t]
                tmp_d['corr_0'] = np.asarray(new_data[nr_ccsd_t]) - np.asarray(new_data[nr_dhf])
                tmp_d['corr_1'] = np.asarray(new_data[sr_nr_ccsd_t]) - np.asarray(new_data[sr_nr_dhf])
                tmp_d['delta_corr_1'] = tmp_d['corr_1'] - tmp_d['corr_0']
                tmp_d['delta_orb_1'] = np.asarray(new_data[sr_nr_dhf]) - np.asarray(new_data[nr_dhf])
                tmp_d['delta_tot_1'] = tmp_d['delta_corr_1'] + tmp_d['delta_orb_1']

                ml_keys = ['{}@{}@'.format(calc_section_key, state) +  get_tex_key(k) for k in tmp_d.keys()]
                for k, new_k in zip(tmp_d.keys(), ml_keys):
                    new_data[new_k] = tmp_d[k]
        else:
            selected_state = 'M_L=0'

        der_keys = ['alpha_0', 'alpha_1', 'alpha_2', 'corr_0', 'corr_1', 'corr_2',
                    'delta_corr_1', 'delta_corr_so', 'delta_corr_2', 'delta_orb_1',
                    'delta_orb_so', 'delta_orb_2', 'delta_tot_1',
                    'delta_tot_so', 'delta_tot_2']
        prefix = '{}@{}@'.format(calc_section_key,selected_state)
        new_der_keys = [prefix + get_tex_key(k)  for k in der_keys]
        for k, new_k in zip(der_keys, new_der_keys):
            new_data[new_k] = data[k]

        # add excite-state properties
        if group in [13]:
            delta_32_12 = r'@\bar{\alpha}_{J=3/2}-\alpha_{J=1/2}'
            new_data['{}@{}'.format(calc_section_key, selected_state) + delta_32_12] = \
                np.asarray(new_data[group_13_32_ave_index]) - \
                np.asarray(new_data[group_13_12_gs_index])
            # alpha_32_12 = r'\alpha_{3/2,1/2}'
            # alpha_32_32 = r'\alpha_{3/2,3/2}'
            # tmp_state = r'^2P_{3/2}'
            # new_data['{}@{}@{}'.format('DC', tmp_state, alpha_32_12)] =
        elif group in [15]:
            delta_32_12 = r'@\alpha_{J=3/2}-\alpha_{J=1/2}'
            new_data['{}@{}'.format(calc_section_key, selected_state) + delta_32_12] = \
                np.asarray(new_data[group_15_32_32_index]) - \
                np.asarray(new_data[group_15_32_12_index])
        elif group in [17]:
            delta_32_12 = r'@\bar{\alpha}_{J=3/2}-\alpha_{J=1/2}'
            new_data['{}@{}'.format(calc_section_key, selected_state) + delta_32_12] = \
                np.asarray(new_data[group_17_32_gs_index]) - \
                np.asarray(new_data[group_17_12_12_index])


        self._table = pd.DataFrame.from_dict(new_data, orient='index',
                                    columns=self.data['atom_list'])
        return self._table

    def to_latex(self):
        table = self.table
        index = table.index

        group = data['group']
        calc_section_key = 'Calc.'
        first_idx = []
        second_idx = []
        third_idx = []

        ave_tag = '${}$'.format(Tabulator.get_gs_term(group, scheme='LS'))
        for i in index:
            a, b, c = i.split('@')
            first_idx.append(a)
            if a in ['DC', 'NR', 'SR-NR', 'SR-DC']:
                new_b = '${}$'.format(b)
            else:
                new_b = b

            if a in [calc_section_key]:
                if b in ['M_L=0', 'M_L=1', ave_tag.strip('$')]:
                    new_b = '${}$'.format(b)
                new_c = '${}$'.format(c)
            else:
                new_c = c
            second_idx.append(new_b)
            third_idx.append(new_c)

        new_index = np.asarray([first_idx, second_idx, third_idx])
        tuples = list(zip(*new_index))
        multi_index = pd.MultiIndex.from_tuples(tuples,
                                                names=['$\hat{H}$', 'State',
                                                       'Method'])
        table.index = multi_index

        gs_term = '$' + Tabulator.get_gs_term(group) + '$'

        group_13_32_ave_index = r'$^2P_{3/2}$'
        group_13_12_gs_index = r'$^2P_{1/2}$'
        group_13_32_12_index = r'$^2P_{3/2}, M_J=1/2$'
        group_13_32_32_index = r'$^2P_{3/2}, M_J=3/2$'
        group_17_32_32_index = r'$^2P_{3/2}, M_J=3/2$'
        group_17_32_12_index = r'$^2P_{3/2}, M_J=1/2$'
        group_17_12_12_index = r'$^2P_{1/2}, M_J=1/2$'
        group_15_32_12_index = r'$^4S_{3/2}, M_J=1/2$'
        group_15_32_32_index = r'$^4S_{3/2}, M_J=3/2$'
        if group in [1, 2, 11, 12, 15, 18]:
            df1 = table.loc[(['NR', 'SR-NR', 'SR-DC', 'DC'],
                             ['$M_L=0$', gs_term, group_15_32_12_index, group_15_32_32_index],
                             # ['DHF', 'CCSD(T)'],
                             ['CCSD(T)', 'MRCISD', ],
                             )]
        else:
            df1 = table.loc[(['NR', 'SR-NR', 'SR-DC', 'DC'],
                             # ['$M_L=0$', '$M_L=1$', 'Average', gs],
                             [ave_tag, gs_term, group_13_32_12_index,
                              group_13_32_32_index, group_13_32_ave_index,
                              group_17_32_32_index, group_17_32_12_index,
                              group_17_12_12_index],
                             # ['DHF', 'CCSD(T)'],
                             ['CCSD(T)', 'MRCISD'],
                             )]
        # df2 = table.loc[[ 'DC']]
        df3 = table.loc[[calc_section_key]]

        atom_list = data['atom_list']
        df_ref = pd.DataFrame.from_dict({
            ('Others', gs_term, r'Ref. \cite{Schwerdtfeger2019}') :
                self.reference_value()},
            orient='index', columns=df1.columns)
        # df_tot = pd.concat([df1, df2, df3])
        df_tot = pd.concat([df1, df_ref, df3])

        caption = r'Static dipole polarizabilities with non-relativistic, ' \
                  r'scalar-relativistic, full Dirac-Coulomb relativistic effects ' \
                  r'of group-{} atoms, as well as their derived properties ' \
                  r'that are defined in Sec. \ref{{sec:deriv_value}}'.format(self.group)
        label = r'tab:dipole_group_{}'.format(self.group)

        tex = df_tot.to_latex(
            longtable=True, float_format="{:0.2f}".format,
            escape=False, caption=caption, label=label).replace('nan', '$--$').replace('NaN', '$--$')
        return tex

    def reference_value(self):
        if self.group == 1:
            ref = ['164.11', '162.7\pm0.5', '289.7\pm0.3', '319.8\pm0.3',
                    '400.9\pm0.7', '317.8\pm2.4']
        elif self.group == 2:
            ref =  ['37.74\pm0.03', '71.2\pm0.4', '160.8\pm4', '197.2\pm0.2',
                    '272.9\pm10', '246.0\pm4']
        elif self.group == 11:
            ref =  ['46.5\pm0.5', '55.0\pm8', '36.0\pm3']
        elif self.group == 12:
            ref =  ['38.67\pm0.3', '46\pm2', '33.91\pm0.34', '28\pm2']
        elif self.group == 13:
            ref =  ['20.5\pm1', '57.8\pm1.0', '50.3\pm3', '65\pm4', '50.0\pm2',
                    '29.2\pm2']
        elif self.group == 14:
            ref = ['11.3\pm0.2', '37.3\pm 0.7', '40.0\pm1', '53.0\pm6',
                    '47.0\pm3', '31.0\pm4']
        elif self.group == 15:
            ref =  ['7.4\pm0.2', '25.0\pm1', '30\pm1', '43\pm2', '48\pm4', '71.0\pm20']
        elif self.group == 16:
            ref = ['5.3\pm0.2', '19.4\pm0.1', '28.9\pm1.0', '38.0\pm4', '44.0\pm4', '--']
        elif self.group == 17:
            ref =  ['3.74\pm0.08', '14.6\pm0.1', '21.0\pm1', '32.9\pm1.3', '42.0\pm4', '76.0\pm15']
        elif self.group == 18:
            ref =  ['1.38', '2.66', '11.08', '16.78\pm0.02', '27.32\pm0.2', '35.0\pm2', '58.0\pm6']
        else:
            raise NotImplementedError('Does not support group-{} elements'.format(self.group))

        return ['${}$'.format(i) for i in ref]


def main():
    import sys
    import os
    import argparse

    parser = argparse.ArgumentParser(description='Script to generate polarizability tex for atomic db')
    parser.add_argument('-q','--quadru', dest='quadru', action='store_true',
                        help='Is this a quadrupole calculation? (default: False)')
    parser.add_argument('--only_best', action='store_true',
                        help='Only show best results (default: False)')
    parser.add_argument('atom_list', nargs='+',
                        help='A list of atoms, atomic nuclei number or atomic symbol')

    args = parser.parse_args()
    data_root = '/Users/yxcheng/PhD/dirac/data_backup/backup_2020_Sep_21'
    calc_path = Path(data_root).resolve()

    with cd(calc_path):
        for e in args.atom_list:
            ele = Element(e)
            # CC-SR
            atomic_tex = []
            footer = ''
            if args.quadru:
                ccnr = ele.symbol + '_q_nr'
                ccso = ele.symbol + '_q_so'
                ccsr = ele.symbol + '_q_theta'
                ciso = ele.symbol + '_q_mrci'
            else:
                ccnr = ele.symbol + '_nr'
                ccso = ele.symbol + '_so'
                ccsr = ele.symbol
                ciso = ele.symbol + '_mrci'

            if os.path.exists(ccnr) and os.path.isdir(ccnr):
                with cd(ccnr):
                    nrcc_processor = NRCCDataProcessor(Path('.').resolve(), has_header=True, only_best=args.only_best)
                    header_tex, body_tex, footer_tex = nrcc_processor.get_tex()
                    # header_tex, body_tex, footer_tex = CCsr_main(Path('.').resolve(), has_header=True, only_best=args.only_best)
                    atomic_tex.append(header_tex)
                    atomic_tex.append(body_tex)
                    # if len(tot_tex) < 1:
                    #     tot_tex.append(header_tex)
                    # tot_tex.append(body_tex)

            if os.path.exists(ele.symbol) and os.path.isdir(ele.symbol):
                with cd(ccsr):
                    srcc_processor = SRCCDataProcessor(Path('.').resolve(), only_best=args.only_best)
                    _header_tex, sr_body_tex, _footer = srcc_processor.get_tex()
                    # _header_tex, sr_body_tex, _footer = CCsr_main(Path('.').resolve(), only_best=args.only_best)
                    atomic_tex.append(sr_body_tex)
            if os.path.exists(ccso) and os.path.isdir(ccso):
                with cd(ccso):
                    if not os.path.exists('cannot'):
                        socc_processor = SOCCDataProcessor(Path('.').resolve(),
                                                           only_best=args.only_best)
                        body_tex = socc_processor.get_tex()
                        # body_tex = CCso_main(Path('.').resolve(), only_best=args.only_best)
                        atomic_tex.append(body_tex)
            if os.path.exists(ciso) and os.path.isdir(ciso):
                with cd(ciso):
                    soci_processor = SOCIDataProcessor(Path('.').resolve(),
                                                       only_best=args.only_best)
                    body_tex = soci_processor.get_tex()
                    # body_tex = CIso_main(Path('.').resolve(), only_best=args.only_best)
                    atomic_tex.append(body_tex)
            atomic_tex.append(footer_tex)
            print(''.join(atomic_tex))


if __name__ == '__main__':
    # main()

    # ----------------------
    # Generate data
    # ----------------------
    data_root = '/Users/yxcheng/PhD/dirac/data_backup/backup_2020_Sep_21'
    dp = DataProcessor(scratch_dir='./figs/', data_root=data_root)
    # dp.dump_data(refresh=True)
    # exit()

    #data = dp.load_data(13, refresh=False)
    # print(Settings(data))

    # ----------------------
    # plot
    # ----------------------
    #for gp in [1, 2] + list(range(11, 19)):
    #    ploter = Ploter.from_group(gp)
    #    ploter.plot()

    # ploter = Ploter(data, save_path='./figs')
    # ploter.plot()

    # -----------------------
    # Output tables
    # -----------------------
    for gp in [1, 2] + list(range(11, 19)):
    #for gp in [17]:
        data = dp.load_data(gp, refresh=False)
        tb = Tabulator(data)
        tex_str = tb.to_latex()

        with open('./figs/summary_tab_group_{}.tex'.format(tb.group), 'w') as f:
            f.write(tex_str)


