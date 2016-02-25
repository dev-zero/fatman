#!/usr/bin/env python

from fatman.models import Pseudopotential, BasisSet, PseudopotentialFamily, BasissetFamily, Method
from ase.units import Ry

basis  = BasissetFamily.get(BasissetFamily.name=='TZVP-pob')
pseudo = PseudopotentialFamily.get(PseudopotentialFamily.name=='ALL')
code   = "cp2k"

Method.create_or_get(basis_set        = basis, 
                     pseudopotential  = pseudo, 
                     code             = code, 
                     settings         = {   "cutoff_rho"  :  1000.*Ry, 
                                            "rel_settings": {"method":"ZORA",
                                                             "zora_type":"MP"     }  , 
                                            "kpoint_shift": [0.5,0.5,0.5],
                                            "qs_settings" : {"method":"GAPW",
                                                             "extrapolation": "USE_GUESS",
                                                             "eps_default": "1.0E-10"},
                                         }
                     )
