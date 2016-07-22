#!/usr/bin/env python

import os
import json
import subprocess
from os import path
from tempfile import mkdtemp

import click
import requests

from jinja2 import Template
from fatman.tools import randomword

PSEUDOPOTENTIAL_URL = '{}/pseudos'
METHOD_URL          = '{}/methods'
TASK_URL            = '{}/tasks'

input_template_str = """
&GLOBAL
  PROGRAM_NAME ATOM
&END GLOBAL
&ATOM
  ELEMENT {{ element }}              

  RUN_TYPE PSEUDOPOTENTIAL_OPTIMIZATION

  CALCULATE_STATES {{ states|length }}
  {% for elconf in states -%}
  ELECTRON_CONFIGURATION {{ elconf }}
  {% endfor -%}
  CORE {{ core }}
  MAX_ANGULAR_MOMENTUM 2

  COULOMB_INTEGRALS ANALYTIC
  EXCHANGE_INTEGRALS ANALYTIC

  &METHOD
     METHOD_TYPE  KOHN-SHAM
     RELATIVISTIC DKH(2)
     &XC
       &XC_FUNCTIONAL PBE 
       &END XC_FUNCTIONAL
     &END XC
  &END METHOD
  &OPTIMIZATION
    EPS_SCF 1.e-10
  &END

  &AE_BASIS
     BASIS_TYPE GEOMETRICAL_GTO
  &END AE_BASIS
  &PP_BASIS
     BASIS_TYPE GEOMETRICAL_GTO
  &END PP_BASIS
  &POTENTIAL
    PSEUDO_TYPE GTH
    CONFINEMENT {{ conf|join(' ') }}
    &GTH_POTENTIAL
        {{ pseudoguess }}
    &END 
  &END POTENTIAL

  &POWELL
     ACCURACY   1.e-10 
     STEP_SIZE  0.5
  &END
&END ATOM
"""
input_template = Template(input_template_str)

upf_template_str = """
&GLOBAL
  PROGRAM_NAME ATOM
&END GLOBAL
&ATOM
  ELEMENT {{ element }} 
  ELECTRON_CONFIGURATION CORE {{ elconf }}
  CORE {{ core }}
  &METHOD
     METHOD_TYPE  KOHN-SHAM
     &XC
       &XC_FUNCTIONAL PBE
       &END XC_FUNCTIONAL
     &END XC
  &END METHOD
  &PP_BASIS
  &END PP_BASIS
  &POTENTIAL
    POTENTIAL_NAME {{ family }}
    PSEUDO_TYPE GTH
    &GTH_POTENTIAL
        {{ pseudoguess }}
    &END GTH_POTENTIAL
  &END POTENTIAL
  &PRINT
    &UPF_FILE ON
       FILENAME =output.UPF
    &END
  &END PRINT
&END ATOM
"""
upf_template = Template(upf_template_str)

ATOMS_DB = {
    'largecore': {
        'Li' : ["{}/Li/q1/GTH-PARAMETER"   ,"2s1"                ,"[He]"          ],
        'Be' : ["{}/Be/q2/GTH-PARAMETER"   ,"2s2"                ,"[He]"          ],
        'Na' : ["{}/Na/q1/GTH-PARAMETER"   ,"3s1"                ,"[Ne]"          ],
        'Mg' : ["{}/Mg/q2/GTH-PARAMETER"   ,"3s2"                ,"[Ne]"          ],
        'K'  : ["{}/K/q1/GTH-PARAMETER"    ,"4s1"                ,"[Ar]"          ],
        'Ca' : ["{}/Ca/q2/GTH-PARAMETER"   ,"4s2"                ,"[Ar]"          ],
        'Sc' : ["{}/Sc/q3/GTH-PARAMETER"   ,"4s2 3d1"            ,"[Ar]"          ],
        'Ti' : ["{}/Ti/q4/GTH-PARAMETER"   ,"4s2 3d2"            ,"[Ar]"          ],
        'V'  : ["{}/V/q5/GTH-PARAMETER"    ,"4s2 3d3"            ,"[Ar]"          ],
        'Cr' : ["{}/Cr/q6/GTH-PARAMETER"   ,"4s2 3d4"            ,"[Ar]"          ],
        'Mn' : ["{}/Mn/q7/GTH-PARAMETER"   ,"4s2 3d5"            ,"[Ar]"          ],
        'Fe' : ["{}/Fe/q8/GTH-PARAMETER"   ,"4s2 3d6"            ,"[Ar]"          ],
        'Co' : ["{}/Co/q9/GTH-PARAMETER"   ,"4s2 3d7"            ,"[Ar]"          ],
        'Ni' : ["{}/Ni/q10/GTH-PARAMETER"  ,"4s2 3d8"            ,"[Ar]"          ],
        'Cu' : ["{}/Cu/q1/GTH-PARAMETER"   ,"4s1"                ,"[Ar] 3d10"     ],
        'Zn' : ["{}/Zn/q2/GTH-PARAMETER"   ,"4s2"                ,"[Ar] 3d10"     ],
        'Ga' : ["{}/Ga/q3/GTH-PARAMETER"   ,"4s2 4p1"            ,"[Ar] 3d10"     ],
        'Rb' : ["{}/Rb/q1/GTH-PARAMETER"   ,"5s1"                ,"[Kr]"          ],
        'Sr' : ["{}/Sr/q2/GTH-PARAMETER"   ,"5s2"                ,"[Kr]"          ],
        'Y'  : ["{}/Y/q3/GTH-PARAMETER"    ,"5s2 4d1"            ,"[Kr]"          ],
        'Zr' : ["{}/Zr/q4/GTH-PARAMETER"   ,"5s2 4d2"            ,"[Kr]"          ],
        'Nb' : ["{}/Nb/q5/GTH-PARAMETER"   ,"5s1 4d4"            ,"[Kr]"          ],
        'Mo' : ["{}/Mo/q6/GTH-PARAMETER"   ,"5s1 4d5"            ,"[Kr]"          ],
        'Tc' : ["{}/Tc/q7/GTH-PARAMETER"   ,"5s2 4d5"            ,"[Kr]"          ],
        'Ru' : ["{}/Ru/q8/GTH-PARAMETER"   ,"5s1 4d7"            ,"[Kr]"          ],
        'Rh' : ["{}/Rh/q9/GTH-PARAMETER"   ,"5s1 4d8"            ,"[Kr]"          ],
        'Pd' : ["{}/Pd/q10/GTH-PARAMETER"  ,"5s1 4d9"            ,"[Kr]"          ],
        'Ag' : ["{}/Ag/q1/GTH-PARAMETER"   ,"5s1"                ,"[Kr] 4d10"     ],
        'Cd' : ["{}/Cd/q2/GTH-PARAMETER"   ,"5s2"                ,"[Kr] 4d10"     ],
        'In' : ["{}/In/q3/GTH-PARAMETER"   ,"5s2 5p1"            ,"[Kr] 4d10"     ],
        'Cs' : ["{}/Cs/q1/GTH-PARAMETER"   ,"6s1"                ,"[Xe]"          ],
        'Ba' : ["{}/Ba/q2/GTH-PARAMETER"   ,"6s2"                ,"[Xe]"          ],
        'Ta' : ["{}/Ta/q5/GTH-PARAMETER"   ,"5d3 6s2"            ,"[Kr] 4d10 4f14 5s2 5p6"],
        'W'  : ["{}/W/q6/GTH-PARAMETER"    ,"5d4 6s2"            ,"[Kr] 4d10 4f14 5s2 5p6"],
        'Re' : ["{}/Re/q7/GTH-PARAMETER"   ,"5d5 6s2"            ,"[Kr] 4d10 4f14 5s2 5p6"],
        'Os' : ["{}/Os/q8/GTH-PARAMETER"   ,"5d6 6s2"            ,"[Kr] 4d10 4f14 5s2 5p6"],
        'Ir' : ["{}/Ir/q9/GTH-PARAMETER"   ,"5d7 6s2"            ,"[Kr] 4d10 4f14 5s2 5p6"],
        'Pt' : ["{}/Pt/q10/GTH-PARAMETER"  ,"5d9 6s1"            ,"[Kr] 4d10 4f14 5s2 5p6"],
        'Au' : ["{}/Au/q1/GTH-PARAMETER"   ,"6s1"                ,"[Kr] 4d10 4f14 5s2 5p6 5d10"],
        'Hg' : ["{}/Hg/q2/GTH-PARAMETER"   ,"6s2"                ,"[Kr] 4d10 4f14 5s2 5p6 5d10"],
        'Tl' : ["{}/Tl/q3/GTH-PARAMETER"   ,"6s2 6p1"            ,"[Kr] 4d10 4f14 5s2 5p6 5d10"],
        'Pb' : ["{}/Pb/q4/GTH-PARAMETER"   ,"6s2 6p2"            ,"[Kr] 4d10 4f14 5s2 5p6 5d10"],
        'Bi' : ["{}/Bi/q5/GTH-PARAMETER"   ,"6s2 6p3"            ,"[Kr] 4d10 4f14 5s2 5p6 5d10"],
        },

    'mediumcore': {
        'Cu' : ["{}/Cu/q11/GTH-PARAMETER"  ,"4s1 3d10"           ,"[Ar]"          ],
        'Zn' : ["{}/Zn/q12/GTH-PARAMETER"  ,"4s2 3d10"           ,"[Ar]"          ],
        'Ag' : ["{}/Ag/q11/GTH-PARAMETER"  ,"5s1 4d10"           ,"[Kr]"          ],
        'Au' : ["{}/Au/q11/GTH-PARAMETER"  ,"6s1 5d10"           ,"[Kr] 4d10 4f14 5s2 5p6"],
        },

    'smallcore': {
        'H'  : ["{}/H/q1/GTH-PARAMETER"    ,"1s1"                ,"none"          ],
        'He' : ["{}/He/q2/GTH-PARAMETER"   ,"1s2"                ,"none"          ],
        'Li' : ["{}/Li/q3/GTH-PARAMETER"   ,"1s2 2s1"            ,"none"          ],
        'Be' : ["{}/Be/q4/GTH-PARAMETER"   ,"1s2 2s2"            ,"none"          ],
        'B'  : ["{}/B/q3/GTH-PARAMETER"    ,"2s2 2p1"            ,"[He]"          ],
        'C'  : ["{}/C/q4/GTH-PARAMETER"    ,"2s2 2p2"            ,"[He]"          ],
        'N'  : ["{}/N/q5/GTH-PARAMETER"    ,"2s2 2p3"            ,"[He]"          ],
        'O'  : ["{}/O/q6/GTH-PARAMETER"    ,"2s2 2p4"            ,"[He]"          ],
        'F'  : ["{}/F/q7/GTH-PARAMETER"    ,"2s2 2p5"            ,"[He]"          ],
        'Ne' : ["{}/Ne/q8/GTH-PARAMETER"   ,"2s2 2p6"            ,"[He]"          ],
        'Na' : ["{}/Na/q9/GTH-PARAMETER"   ,"2s2 2p6 3s1"        ,"[He]"          ],
        'Mg' : ["{}/Mg/q10/GTH-PARAMETER"  ,"2s2 2p6 3s2"        ,"[He]"          ],
        'Al' : ["{}/Al/q3/GTH-PARAMETER"   ,"3s2 3p1"            ,"[Ne]"          ],
        'Si' : ["{}/Si/q4/GTH-PARAMETER"   ,"3s2 3p2"            ,"[Ne]"          ],
        'P'  : ["{}/P/q5/GTH-PARAMETER"    ,"3s2 3p3"            ,"[Ne]"          ],
        'S'  : ["{}/S/q6/GTH-PARAMETER"    ,"3s2 3p4"            ,"[Ne]"          ],
        'Cl' : ["{}/Cl/q7/GTH-PARAMETER"   ,"3s2 3p5"            ,"[Ne]"          ],
        'Ar' : ["{}/Ar/q8/GTH-PARAMETER"   ,"3s2 3p6"            ,"[Ne]"          ],
        'K'  : ["{}/K/q9/GTH-PARAMETER"    ,"3s2 3p6 4s1"        ,"[Ne]"          ],
        'Ca' : ["{}/Ca/q10/GTH-PARAMETER"  ,"3s2 3p6 4s2"        ,"[Ne]"          ],
        'Sc' : ["{}/Sc/q11/GTH-PARAMETER"  ,"3s2 3p6 4s2 3d1"    ,"[Ne]"          ],
        'Ti' : ["{}/Ti/q12/GTH-PARAMETER"  ,"3s2 3p6 4s2 3d2"    ,"[Ne]"          ],
        'V'  : ["{}/V/q13/GTH-PARAMETER"   ,"3s2 3p6 4s2 3d3"    ,"[Ne]"          ],
        'Cr' : ["{}/Cr/q14/GTH-PARAMETER"  ,"3s2 3p6 4s2 3d4"    ,"[Ne]"          ],
        'Mn' : ["{}/Mn/q15/GTH-PARAMETER"  ,"3s2 3p6 4s2 3d5"    ,"[Ne]"          ],
        'Fe' : ["{}/Fe/q16/GTH-PARAMETER"  ,"3s2 3p6 4s2 3d6"    ,"[Ne]"          ],
        'Co' : ["{}/Co/q17/GTH-PARAMETER"  ,"3s2 3p6 4s2 3d7"    ,"[Ne]"          ],
        'Ni' : ["{}/Ni/q18/GTH-PARAMETER"  ,"3s2 3p6 4s2 3d8"    ,"[Ne]"          ],
        'Cu' : ["{}/Cu/q19/GTH-PARAMETER"  ,"3s2 3p6 4s1 3d10"   ,"[Ne]"          ],
        'Zn' : ["{}/Zn/q20/GTH-PARAMETER"  ,"3s2 3p6 4s2 3d10"   ,"[Ne]"          ],
        'Ga' : ["{}/Ga/q13/GTH-PARAMETER"  ,"3d10 4s2 4p1"       ,"[Ar]"          ],
        'Ge' : ["{}/Ge/q4/GTH-PARAMETER"   ,"4s2 4p2"            ,"[Ar] 3d10"     ],
        'As' : ["{}/As/q5/GTH-PARAMETER"   ,"4s2 4p3"            ,"[Ar] 3d10"     ],
        'Se' : ["{}/Se/q6/GTH-PARAMETER"   ,"4s2 4p4"            ,"[Ar] 3d10"     ],
        'Br' : ["{}/Br/q7/GTH-PARAMETER"   ,"4s2 4p5"            ,"[Ar] 3d10"     ],
        'Kr' : ["{}/Kr/q8/GTH-PARAMETER"   ,"4s2 4p6"            ,"[Ar] 3d10"     ],
        'Rb' : ["{}/Rb/q9/GTH-PARAMETER"   ,"4s2 4p6 5s1"        ,"[Ar] 3d10"     ],
        'Sr' : ["{}/Sr/q10/GTH-PARAMETER"  ,"4s2 4p6 5s2"        ,"[Ar] 3d10"     ],
        'Y'  : ["{}/Y/q11/GTH-PARAMETER"   ,"4s2 4p6 5s2 4d1"    ,"[Ar] 3d10"     ],
        'Zr' : ["{}/Zr/q12/GTH-PARAMETER"  ,"4s2 4p6 5s2 4d2"    ,"[Ar] 3d10"     ],
        'Nb' : ["{}/Nb/q13/GTH-PARAMETER"  ,"4s2 4p6 5s2 4d3"    ,"[Ar] 3d10"     ],
        'Mo' : ["{}/Mo/q14/GTH-PARAMETER"  ,"4s2 4p6 5s1 4d5"    ,"[Ar] 3d10"     ],
        'Tc' : ["{}/Tc/q15/GTH-PARAMETER"  ,"4s2 4p6 5s2 4d5"    ,"[Ar] 3d10"     ],
        'Ru' : ["{}/Ru/q16/GTH-PARAMETER"  ,"4s2 4p6 5s2 4d6"    ,"[Ar] 3d10"     ],
        'Rh' : ["{}/Rh/q17/GTH-PARAMETER"  ,"4s2 4p6 5s2 4d7"    ,"[Ar] 3d10"     ],
        'Pd' : ["{}/Pd/q18/GTH-PARAMETER"  ,"4s2 4p6 5s2 4d8"    ,"[Ar] 3d10"     ],
        'Ag' : ["{}/Ag/q19/GTH-PARAMETER"  ,"5s1 4d10"           ,"[Kr]"          ],
        'Cd' : ["{}/Cd/q12/GTH-PARAMETER"  ,"5s2 4d10"           ,"[Kr]"          ],
        'In' : ["{}/In/q13/GTH-PARAMETER"  ,"5s2 5p1 4d10"       ,"[Kr]"          ],
        'Sn' : ["{}/Sn/q4/GTH-PARAMETER"   ,"5s2 5p2"            ,"[Kr] 4d10"     ],
        'Sb' : ["{}/Sb/q5/GTH-PARAMETER"   ,"5s2 5p3"            ,"[Kr] 4d10"     ],
        'Te' : ["{}/Te/q6/GTH-PARAMETER"   ,"5s2 5p4"            ,"[Kr] 4d10"     ],
        'I'  : ["{}/I/q7/GTH-PARAMETER"    ,"5s2 5p5"            ,"[Kr] 4d10"     ],
        'Xe' : ["{}/Xe/q8/GTH-PARAMETER"   ,"5s2 5p6"            ,"[Kr] 4d10"     ],
        'Cs' : ["{}/Cs/q9/GTH-PARAMETER"   ,"5s2 5p6 6s1"        ,"[Kr] 5d10"     ],
        'Ba' : ["{}/Ba/q10/GTH-PARAMETER"  ,"5s2 5p6 6s2"        ,"[Kr] 5d10"     ],

        'La' : ["{}/La/q11/GTH-PARAMETER"  ," 4f1 5s2 5p6 5d0 6s2"  ,"[Kr] 4d10"     ],
        'Ce' : ["{}/Ce/q12/GTH-PARAMETER"  ," 4f1 5s2 5p6 5d1 6s2"  ,"[Kr] 4d10"     ],
        'Pr' : ["{}/Pr/q13/GTH-PARAMETER"  ," 4f3 5s2 5p6 5d0 6s2"  ,"[Kr] 4d10"     ],
        'Nd' : ["{}/Nd/q14/GTH-PARAMETER"  ," 4f4 5s2 5p6 5d0 6s2"  ,"[Kr] 4d10"     ],
        'Pm' : ["{}/Pm/q15/GTH-PARAMETER"  ," 4f5 5s2 5p6 5d0 6s2"  ,"[Kr] 4d10"     ],
        'Sm' : ["{}/Sm/q16/GTH-PARAMETER"  ," 4f6 5s2 5p6 5d0 6s2"  ,"[Kr] 4d10"     ],
        'Eu' : ["{}/Eu/q17/GTH-PARAMETER"  ," 4f7 5s2 5p6 5d0 6s2"  ,"[Kr] 4d10"     ],
        'Gd' : ["{}/Gd/q18/GTH-PARAMETER"  ," 4f7 5s2 5p6 5d1 6s2"  ,"[Kr] 4d10"     ],
        'Tb' : ["{}/Tb/q19/GTH-PARAMETER"  ," 4f9 5s2 5p6 5d0 6s2"  ,"[Kr] 4d10"     ],
        'Dy' : ["{}/Dy/q20/GTH-PARAMETER"  ," 4f10 5s2 5p6 5d0 6s2" ,"[Kr] 4d10"     ],
        'Ho' : ["{}/Ho/q21/GTH-PARAMETER"  ," 4f11 5s2 5p6 5d0 6s2" ,"[Kr] 4d10"     ],
        'Er' : ["{}/Er/q22/GTH-PARAMETER"  ," 4f12 5s2 5p6 5d0 6s2" ,"[Kr] 4d10"     ],
        'Tm' : ["{}/Tm/q23/GTH-PARAMETER"  ," 4f13 5s2 5p6 5d0 6s2" ,"[Kr] 4d10"     ],
        'Yb' : ["{}/Yb/q24/GTH-PARAMETER"  ," 4f14 5s2 5p6 5d0 6s2" ,"[Kr] 4d10"     ],

        'Hf' : ["{}/Hf/q12/GTH-PARAMETER"  ,"5s2 5p6 5d2 6s2"    ,"[Kr] 4d10 4f14"     ],
        'Ta' : ["{}/Ta/q13/GTH-PARAMETER"  ,"5s2 5p6 5d3 6s2"    ,"[Kr] 4d10 4f14"     ],
        'W'  : ["{}/W/q14/GTH-PARAMETER"   ,"5s2 5p6 5d4 6s2"    ,"[Kr] 4d10 4f14"     ],
        'Re' : ["{}/Re/q15/GTH-PARAMETER"  ,"5s2 5p6 5d5 6s2"    ,"[Kr] 4d10 4f14"     ],
        'Os' : ["{}/Os/q16/GTH-PARAMETER"  ,"5s2 5p6 5d6 6s2"    ,"[Kr] 4d10 4f14"     ],
        'Ir' : ["{}/Ir/q17/GTH-PARAMETER"  ,"5s2 5p6 5d7 6s2"    ,"[Kr] 4d10 4f14"     ],
        'Pt' : ["{}/Pt/q18/GTH-PARAMETER"  ,"5s2 5p6 5d9 6s1"    ,"[Kr] 4d10 4f14"     ],
        'Au' : ["{}/Au/q19/GTH-PARAMETER"  ,"5s2 5p6 5d10 6s1"   ,"[Kr] 4d10 4f14"     ],
        'Hg' : ["{}/Hg/q12/GTH-PARAMETER"  ,"5d10 6s2"           ,"[Kr] 4d10 4f14 5s2 5p6"],
        'Tl' : ["{}/Tl/q13/GTH-PARAMETER"  ,"5d10 6s2 6p1"       ,"[Kr] 4d10 4f14 5s2 5p6"],
        'Pb' : ["{}/Pb/q14/GTH-PARAMETER"  ,"5d10 6s2 6p2"       ,"[Kr] 4d10 4f14 5s2 5p6"],
        'Bi' : ["{}/Bi/q15/GTH-PARAMETER"  ,"5d10 6s2 6p3"       ,"[Kr] 4d10 4f14 5s2 5p6"],
        'Po' : ["{}/Po/q6/GTH-PARAMETER"   ,"6s2 6p4"            ,"[Kr] 4d10 4f14 5s2 5p6 5d10"],
        'At' : ["{}/At/q7/GTH-PARAMETER"   ,"6s2 6p5"            ,"[Kr] 4d10 4f14 5s2 5p6 5d10"],
        'Rn' : ["{}/Rn/q8/GTH-PARAMETER"   ,"6s2 6p6"            ,"[Kr] 4d10 4f14 5s2 5p6 5d10"],
        },
    }

DELTATEST_ELEMENTS = {
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Rn",
    }

@click.command()
@click.argument('ppname', type=str)
@click.argument('pplocation', type=click.Path(exists=True))
@click.argument('atomtype', type=click.Choice(['smallcore', 'mediumcore', 'largecore']))
@click.argument('elements', type=str, nargs=-1, metavar='[ELEMENTS]')
@click.option('--url', type=str, default='http://localhost:5000',
              help='The URL where FATMAN is running (default: http://localhost:5000/)')
@click.option('--online', is_flag=True, help='actually submit something')
@click.option('--cp2k-exe', type=str, default='cp2k.sopt', help='path the CP2K executable (default: cp2k.sopt)')
@click.option('--ignore-missing-pp', is_flag=True, help='ignore elements for which pseudo coefficients are missing')
def run(ppname, pplocation, atomtype, elements, url, online, cp2k_exe, ignore_missing_pp):
    """take pseudopotentials optimized by cp2k and upload them for deltatesting

    This script goes through a list of given elements, takes the stored cp2k-optimized pseudo from the filesystem,
    uploads it, and creates a FATMAN Method that uses this PP and new FATMAN Tasks for deltatesting the pseudo.
    Most options are hardcoded, including the electron configuration of the elements and the location of
    the optimized pseudos (a file usually called `GTH-PARAMETER`, in cp2k format).

    Some remarks:
        - In principle, existing pseudopotentials in the same family will not be overwritten, instead raising an error.
          This can be changed manually in the code.
        - Make sure that all instances of the PP family name are consistent (`ppname`) and reflect the
          core size of the desired PP (`atomtype`).
        - The script runs cp2k with a hardcoded input template to convert the cp2k PP to the UPF format. cp2k has to be
          available and runnable. The conversion has been shown to work correctly, but more rigorous testing wouldn't hurt.
    """

    sess = requests.Session()
    sess.verify = False

    atoms_db = ATOMS_DB[atomtype]

    workdir_base = mkdtemp(prefix="newpseudos_", dir=pplocation)
    print("converted pseudos will be located in {}".format(workdir_base))

    cp2k_command = "{} -i {{}} -o {{}}".format(cp2k_exe)

    if len(elements) == 0:
        elements = [e for e in atoms_db.keys() if e in DELTATEST_ELEMENTS]

    print("process elements", ', '.join(elements))
    for el in elements:
       #retrieve a default pseudopotential for that element as a starting guess
       #myguess = sess.get(PSEUDOPOTENTIAL_URL, data={'family':'GTH-PBE', 'element':el}).json()[el]

        #prepare the input file
        settings = {'element' : el,
                    'family'  : ppname,
                    'pseudoguess' : '',
                    'elconf'      : atoms_db[el][1],
                    'core'        : atoms_db[el][2],
                   }
     # #myinput = input_template.render(settings)
     #  with open(atoms_db[el][0].format(pplocation)) as infile:
     #      myinput_lines = infile.readlines()
     #      myinput = "".join(myinput_lines)
     #

        #prepare the working directory
        workdir = path.join(workdir_base, "{:s}".format(el))

        try:
            os.makedirs(workdir)
        except OSError:
            pass

        os.chdir(workdir)
     #  filename = path.dirname(atoms_db[el][0].format(pplocation))
     #  os.system("cp {:}/GTH-PARAMETER {:}".format(filename, workdir))


     #  #write the input to a file
     #  of = open("job.inp", "w")
     #  of.write(myinput)
     #  of.close()

        #and run it.
        #os.popen(cp2k_command + " -i job.inp -o job.out")
        #myinput_lines = []
        #open the resulting pseudopotential data
        #infile = open("GTH-PARAMETER")
        #pseudo_data = myguess #"".join(infile.readlines()[1:] + ['# ' + x for x in myinput_lines])

        try:
            infile = open(atoms_db[el][0].format(pplocation), 'r')
        except FileNotFoundError:
            if ignore_missing_pp:
                print("{}: ignoring element due to missing pseudopotential coefficients".format(el))
                continue
            else:
                raise

        pseudo_data = "".join(infile.readlines()[1:])

        #upload the PP data.
        if online:
            req = sess.post(PSEUDOPOTENTIAL_URL.format(url),
                            data={'family': ppname,
                                  'element': el,
                                  'pseudo': pseudo_data,
                                  'format': 'CP2K',
                                  'overwrite': False,
                                 })
            req.raise_for_status()
        else:
            print("{}: would have submitted a pseudopotential in CP2K format".format(el))

        settings['pseudoguess'] = pseudo_data

        #now for the UPF conversion
        myinput = upf_template.render(settings)
        of = open("upf.inp", "w")
        of.write(myinput)
        of.close()

        #and run it.
        print("{}: converting pseudopotential from CP2K to UPF".format(el))
        subprocess.check_call(cp2k_command.format("upf.inp", "upf.out"), shell=True)

        #open the resulting UPF file
        infile = open("output.UPF", 'r')
        upf_data = infile.read()

        if online:
            req = sess.post(PSEUDOPOTENTIAL_URL.format(url),
                            data={'family': ppname,
                                  'element': el,
                                  'pseudo': upf_data,
                                  'format': 'UPF',
                                  'overwrite': True})
            req.raise_for_status()
            upf_id = req.json()['id']

            #now, create a 'method' with this PP.
            req = sess.post(METHOD_URL.format(url),
                            data={'code':'espresso',
                                  'basis_set': 55,
                                  'pseudopotential': upf_id,
                                  'settings': json.dumps(
                                      {'smearing': 'marzari-vanderbilt',
                                       'xc': 'PBE',
                                       'cutoff_pw': 3401.4244569396633,
                                       'sigma': 0.027211395655517306,
                                       'max_kpoints': [20, 20, 20]})})
            req.raise_for_status()
            method_id = req.json()['id']

            #create the corresponding deltatest tasks for this method.
            req = sess.post(TASK_URL.format(url),
                            data={'method':method_id,
                                  'test': 'deltatest_{}'.format(el),
                                  'priority': 51})
            req.raise_for_status()
            tasklist = req.json()

            print("{:2s}: {:}".format(el, ", ".join([str(x) for x in tasklist])))

        else:
            print("{}: would have submitted the UPF pseudo, a new method entry and tasks".format(el))

if __name__ == "__main__":
    run()
