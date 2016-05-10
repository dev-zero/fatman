#!/usr/bin/env python
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
    
    print "########################################"
    print "#       CP2K PSEUDOPOTENTIAL FILE      #"
    print "########################################"
    print "# generated on: {:>22s} #".format(datetime.now().strftime("%d-%b-%Y, %H:%M"))
    print "########################################"
    for el in atoms_db:
        if el in pseudos.keys():
            line1  = pseudos[el].split("\n")[0]
            nelec = sum([int(x) for x in line1.split()])
            print "#"
            print "{0:}  {1:}  {1:}-q{2:}".format(el, prefix,nelec)
            print pseudos[el].rstrip()
    
def randomword(length):
    import random, string
    return ''.join(random.choice(string.lowercase + string.digits) for i in range(length))

if __name__ == "__main__":
    main(argv[1:])
