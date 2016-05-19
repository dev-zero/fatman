#!/usr/bin/env python
"""add_method.py -- Create a Method entry in the FATMAN DB

All the method settings, basis set, pseudo, etc. are hardcoed in here.
Modify the script and run it on the FATMAN server.

A method is only created if it is not present in the DB yet.

Parameters:
    -h          show this help.
"""


from __future__ import print_function
from sys import argv
from fatman.models import Pseudopotential, BasisSet, PseudopotentialFamily, BasissetFamily, Method
from ase.units import Ry


def main():
    basis  = BasissetFamily.get(BasissetFamily.name=='DZVP-MOLOPT-SR-GTH')
    pseudo = PseudopotentialFamily.get(PseudopotentialFamily.name=='GTH-PBE')
    code   = "cp2k"

    M, c = Method.create_or_get(basis_set        = basis, 
                         pseudopotential  = pseudo, 
                         code             = code, 
                         settings         = {   "cutoff_rho"  :   1000.*Ry, 
                                           #    "rel_settings":  {"method":"ZORA",
                                           #                      "zora_type":"MP",
                                           #                      "transformation" : "ATOM"}  , 
                                                "kpoint_shift":  [0.5,0.5,0.5],
                                                "qs_settings" :  {"method":"GAPW",
                                                                  "extrapolation": "USE_GUESS",
                                                                  "eps_default": "1.0E-10",
                                                                  "epsiso"     : "1.0E-06"},
                                           #    "kind_settings": {"lebedev_grid": "80",
                                           #                      "radial_grid" : "80" }
                                             }
                         )

    if c:
        print ("CREATED METHOD: {:}".format(M.id))
    else:
        print ("Nothing created, id of existing method: {:}".format(M.id))

if __name__ == "__main__":
    if len(argv) == 1 and '-h' not in argv:
        main(argv[1:])
    else:
        print (__doc__)
