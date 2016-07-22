#!/usr/bin/env python

from __future__ import print_function
from datetime import datetime

import requests
import click

PSEUDOPOTENTIAL_URL = '{}/pseudos'

atoms_db =[ 'H'                                                                                                                                   , 'He'  , 
            'Li'  , 'Be'  ,                                                                                 'B'   , 'C'   , 'N'   , 'O'   , 'F'   , 'Ne'  , 
            'Na'  , 'Mg'  ,                                                                                 'Al'  , 'Si'  , 'P'   , 'S'   , 'Cl'  , 'Ar'  , 
            'K'   , 'Ca'  , 'Sc'  , 'Ti'  , 'V'   , 'Cr'  , 'Mn'  , 'Fe'  , 'Co'  , 'Ni'  , 'Cu'  , 'Zn'  , 'Ga'  , 'Ge'  , 'As'  , 'Se'  , 'Br'  , 'Kr'  ,
            'Rb'  , 'Sr'  , 'Y'   , 'Zr'  , 'Nb'  , 'Mo'  , 'Tc'  , 'Ru'  , 'Rh'  , 'Pd'  , 'Ag'  , 'Cd'  , 'In'  , 'Sn'  , 'Sb'  , 'Te'  , 'I'   , 'Xe'  ,
            'Cs'  , 'Ba'  , 'La'  , 'Hf'  , 'Ta'  , 'W'   , 'Re'  , 'Os'  , 'Ir'  , 'Pt'  , 'Au'  , 'Hg'  , 'Tl'  , 'Pb'  , 'Bi'  , 'Po'  , 'At'  , 'Rn'  ]

@click.command()
@click.argument('family', type=str)
@click.argument('prefix', type=str)
@click.option('--url', type=str, default='http://127.0.0.1:5000',
              help='The URL where FATMAN is running (default: http://localhost:5000/)')
def store_cp2k_pseudos(family, prefix, url):
    """Retrieve Pseudopotentials from FATMAN and format for cp2k.

    Download the entire set of pseudopotentials of the specified FAMILY via
    the FATMAN API and format them in a cp2k-compatible way.
    The pseudos will be named PREFIX and PREFIX-qX.

    Output is written directly to STDOUT, should be redirected to an appropriate file.

    Currently multiple pseudos per element are not supported (e.g. Li-q1 and Li-q3), only one PP is
    retrieved per element and stored.
    """

    req = requests.get(PSEUDOPOTENTIAL_URL.format(url), data=dict(family=family, format='CP2K'))
    req.raise_for_status()

    # the combination of (family,element,format) is guaranteed to be unique
    pseudos = {p['element']: p for p in req.json()}

    print("########################################")
    print("#       CP2K PSEUDOPOTENTIAL FILE      #")
    print("########################################")
    print("# generated on: {:>22s} #".format(datetime.now().strftime("%d-%b-%Y, %H:%M")))
    print("########################################")
    for sym in atoms_db:
        if sym in pseudos.keys():
            line1 = pseudos[sym]['pseudo'].split("\n")[0]
            nelec = sum([int(x) for x in line1.split()])

            print("#")
            print("{0:}  {1:}-q{2:}  {1:}".format(sym, prefix, nelec))
            print(pseudos[sym]['pseudo'].rstrip())
        else:
            print("#")
            print("# no pseudo available for {}".format(sym))

if __name__ == "__main__":
    store_cp2k_pseudos()
