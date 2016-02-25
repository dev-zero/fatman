#!/usr/bin/env python

from fatman.models import Pseudopotential, BasisSet, PseudopotentialFamily, BasissetFamily, Method
from ase.units import Ry

basis  = BasissetFamily.get(BasissetFamily.name=='DZVP-MOLOPT-GTH')
pseudo = PseudopotentialFamily.get(PseudopotentialFamily.name=='GTH-PBE')
code   = "cp2k"

Method.create_or_get(basis_set = basis, pseudopotential=pseudo, code=code, settings={"cutoff_rho":1000.*Ry, "GAPW": "false", "relativistic": "false", "kpoint_shift": [0.5,0.5,0.5]})
