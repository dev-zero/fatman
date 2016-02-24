#!/usr/bin/env python

from fatman.models import Pseudopotential, BasisSet, PseudopotentialFamily, BasissetFamily, Method

basis  = BasissetFamily.get(BasissetFamily.name=='DZVP-MOLOPT-GTH')
pseudo = PseudopotentialFamily.get(PseudopotentialFamily.name=='GTH-PBE')
code   = "cp2k"

Method.create_or_get(basis_set = basis, pseudopotential=pseudo, code=code, settings={"cutoff_rho":19047.976958862117, "GAPW": "false", "relativistic": "false"})
