#!/usr/bin/env python
"""store_cp2k_pp.py <FAMILY> <PREFIX> [-h] - Retrieve Pseudopotentials from FATMAN and format for cp2k.

Download the entire set of pseudopotentials of the specified family via the FATMAN API and format
them in a cp2k-compatible way. Output is written directly to STDOUT, should be redirected to an 
appropriate file.

Currently multiple pseudos per element are not supported (e.g. Li-q1 and Li-q3), only one PP is
retrieved per element and stored.

Parameters
    FAMILY   The name of the PP family stored in the FATMAN DB
    PREFIX   The PP prefix which is to be used in the cp2k format
             The pseudos will then be named PREFIX and PREFIX-qX
    -h       Show this help.
"""


from __future__ import print_function
from sys import argv
from ase.calculators.cp2k import CP2K
from ase import Atoms
from datetime import datetime
import requests, os, json, hashlib
from base64 import urlsafe_b64encode
from jinja2 import Template

SERVER = 'https://172.23.64.223'
SERVER = 'http://127.0.0.1:5001'
PSEUDOPOTENTIAL_URL   = SERVER + '/pseudo'
METHOD_URL            = SERVER + '/method'
TASK_URL              = SERVER + '/tasks/0'


atoms_db =[ 'H'                                                                                                                                   , 'He'  , 
            'Li'  , 'Be'  ,                                                                                 'B'   , 'C'   , 'N'   , 'O'   , 'F'   , 'Ne'  , 
            'Na'  , 'Mg'  ,                                                                                 'Al'  , 'Si'  , 'P'   , 'S'   , 'Cl'  , 'Ar'  , 
            'K'   , 'Ca'  , 'Sc'  , 'Ti'  , 'V'   , 'Cr'  , 'Mn'  , 'Fe'  , 'Co'  , 'Ni'  , 'Cu'  , 'Zn'  , 'Ga'  , 'Ge'  , 'As'  , 'Se'  , 'Br'  , 'Kr'  ,
            'Rb'  , 'Sr'  , 'Y'   , 'Zr'  , 'Nb'  , 'Mo'  , 'Tc'  , 'Ru'  , 'Rh'  , 'Pd'  , 'Ag'  , 'Cd'  , 'In'  , 'Sn'  , 'Sb'  , 'Te'  , 'I'   , 'Xe'  ,
            'Cs'  , 'Ba'  , 'La'  , 'Hf'  , 'Ta'  , 'W'   , 'Re'  , 'Os'  , 'Ir'  , 'Pt'  , 'Au'  , 'Hg'  , 'Tl'  , 'Pb'  , 'Bi'  , 'Po'  , 'At'  , 'Rn'  ]


def main(args):
    family = args[0]
    prefix = args[1]

    req = requests.get(PSEUDOPOTENTIAL_URL, data = dict(family = family))
    req.raise_for_status()

    pseudos = req.json()
    
    print("########################################")
    print("#       CP2K PSEUDOPOTENTIAL FILE      #")
    print("########################################")
    print("# generated on: {:>22s} #".format(datetime.now().strftime("%d-%b-%Y, %H:%M")))
    print("########################################")
    for el in atoms_db:
        if el in pseudos.keys():
            line1  = pseudos[el].split("\n")[0]
            nelec = sum([int(x) for x in line1.split()])
            print("#")
            print("{0:}  {1:}  {1:}-q{2:}".format(el, prefix,nelec))
            print(pseudos[el].rstrip())
    
def randomword(length):
    import random
    import string
    return ''.join(random.choice(string.lowercase + string.digits) for i in range(length))

if __name__ == "__main__":
    if len(argv) == 3 and '-h' not in argv:
        main(argv[1:])
    else:
        print(__doc__)
