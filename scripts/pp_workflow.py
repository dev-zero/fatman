#!/usr/bin/env python
from sys import argv
from ase.calculators.cp2k import CP2K
from ase import Atoms
import requests, os, json

SERVER = 'https://172.23.64.223'
SERVER = 'http://127.0.0.1:5001'
PSEUDOPOTENTIAL_URL   = SERVER + '/pseudo'
METHOD_URL            = SERVER + '/method'
TASK_URL              = SERVER + '/tasks/0'


input_template = """
&GLOBAL
  PROGRAM_NAME ATOM
&END GLOBAL
&ATOM
  ELEMENT {element:s} 

  RUN_TYPE PSEUDOPOTENTIAL_OPTIMIZATION

  ELECTRON_CONFIGURATION {elconf:s}
  CORE {core:s}
  MAX_ANGULAR_MOMENTUM {maxang:d}

  COULOMB_INTEGRALS ANALYTIC
  EXCHANGE_INTEGRALS ANALYTIC

  &METHOD
     METHOD_TYPE  KOHN-SHAM
     RELATIVISTIC DKH(2)
     &XC
       &XC_FUNCTIONAL {xc:s} 
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
    &GTH_POTENTIAL
        {pseudoguess:s}
    &END 
  &END POTENTIAL

  &POWELL
     ! ACCURACY   1.e-10
     ACCURACY   1.e-0            !probably needs fixing!
     STEP_SIZE  0.5
  &END
&END ATOM
"""

upf_template = """
&GLOBAL
  PROGRAM_NAME ATOM
&END GLOBAL
&ATOM
  ELEMENT {element:s} 
  ELECTRON_CONFIGURATION  {elconf:s}
  CORE {core:s} 
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
    POTENTIAL_NAME {family:s}
    PSEUDO_TYPE GTH
    &GTH_POTENTIAL
        {pseudoguess:s}
    &END GTH_POTENTIAL
  &END POTENTIAL
  &PRINT
    &UPF_FILE ON
       FILENAME =output.UPF
    &END
  &END PRINT
&END ATOM
"""



rc_factor = conf_n_d  = conf_a_d = 1.

                       #0         1         2         3         4            5                             6                 7                        8                9  10  11    12
atoms_db =          {  #s         p         d         f                                                                      valence                  pseudo         n[s   p   d]   core 
            'H'  : [    1    ,    0    ,    0    ,    0    ,    conf_a_d   , 0.31  /0.5291*rc_factor   ,  6+conf_n_d*0  ,   "1s1"                ,   "GTH-PBE-q1"    , 1 , 0 , 0 , "none"           ] , 
            'He' : [    1    ,    0    ,    0    ,    0    ,    conf_a_d   , 0.28  /0.5291*rc_factor   ,    conf_n_d    ,   "1s2"                ,   "GTH-PBE-q2"    , 1 , 0 , 0 , "none"           ] , 
            'Li' : [    2    ,    0    ,    0    ,    0    ,    conf_a_d   , 1.28  /0.5291*rc_factor   ,    conf_n_d    ,   "1s2 2s1"            ,   "GTH-PBE-q3"    , 2 , 0 , 0 , "none"           ] , 
            'Be' : [    2    ,    0    ,    0    ,    0    ,    conf_a_d   , 0.96  /0.5291*rc_factor   ,    conf_n_d    ,   "1s2 2s2"            ,   "GTH-PBE-q4"    , 2 , 0 , 0 , "none"           ] , 
            'B'  : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 0.84  /0.5291*rc_factor   ,    conf_n_d    ,   "2s2 2p1"            ,   "GTH-PBE-q3"    , 2 , 2 , 0 , "[He]"           ] , 
            'C'  : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 0.76  /0.5291*rc_factor   ,    conf_n_d    ,   "2s2 2p2"            ,   "GTH-PBE-q4"    , 2 , 2 , 0 , "[He]"           ] , 
            'N'  : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 0.71  /0.5291*rc_factor   ,    conf_n_d    ,   "2s2 2p3"            ,   "GTH-PBE-q5"    , 2 , 2 , 0 , "[He]"           ] , 
            'O'  : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 0.66  /0.5291*rc_factor   ,    conf_n_d    ,   "2s2 2p4"            ,   "GTH-PBE-q6"    , 2 , 2 , 0 , "[He]"           ] , 
            'F'  : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 0.57  /0.5291*rc_factor   ,    conf_n_d    ,   "2s2 2p5"            ,   "GTH-PBE-q7"    , 2 , 2 , 0 , "[He]"           ] , 
            'Ne' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 0.58  /0.5291*rc_factor   ,    conf_n_d    ,   "2s2 2p6"            ,   "GTH-PBE-q8"    , 2 , 2 , 0 , "[He]"           ] , 
            'Na' : [    2    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.66  /0.5291*rc_factor   ,    conf_n_d    ,   "2s2 2p6 3s1"        ,   "GTH-PBE-q9"    , 3 , 2 , 0 , "[He]"           ] , 
            'Mg' : [    2    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.41  /0.5291*rc_factor   ,    conf_n_d    ,   "2s2 2p6 3s2"        ,   "GTH-PBE-q10"   , 3 , 2 , 0 , "[He]"           ] , 
            'Al' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.21  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p1"            ,   "GTH-PBE-q3"    , 3 , 3 , 0 , "[Ne]"           ] , 
            'Si' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.11  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p2"            ,   "GTH-PBE-q4"    , 3 , 3 , 0 , "[Ne]"           ] , 
            'P'  : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.07  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p3"            ,   "GTH-PBE-q5"    , 3 , 3 , 0 , "[Ne]"           ] , 
            'S'  : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.05  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p4"            ,   "GTH-PBE-q6"    , 3 , 3 , 0 , "[Ne]"           ] , 
            'Cl' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.02  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p5"            ,   "GTH-PBE-q7"    , 3 , 3 , 0 , "[Ne]"           ] , 
            'Ar' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.06  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p6"            ,   "GTH-PBE-q8"    , 3 , 3 , 0 , "[Ne]"           ] , 
            'K'  : [    2    ,    1    ,    0    ,    0    ,    conf_a_d   , 2.03  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p6 4s1"        ,   "GTH-PBE-q9"    , 4 , 3 , 0 , "[Ne]"           ] , 
            'Ca' : [    2    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.76  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p6 4s2"        ,   "GTH-PBE-q10"   , 4 , 3 , 0 , "[Ne]"           ] , 
            'Sc' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.70  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p6 4s2 3d1"    ,   "GTH-PBE-q11"   , 4 , 3 , 3 , "[Ne]"           ] , 
            'Ti' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.60  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p6 4s2 3d2"    ,   "GTH-PBE-q12"   , 4 , 3 , 3 , "[Ne]"           ] , 
            'V'  : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.53  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p6 4s2 3d3"    ,   "GTH-PBE-q13"   , 4 , 3 , 3 , "[Ne]"           ] , 
            'Cr' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.39  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p6 4s1 3d5"    ,   "GTH-PBE-q14"   , 4 , 3 , 3 , "[Ne]"           ] , 
            'Mn' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.39  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p6 4s2 3d5"    ,   "GTH-PBE-q15"   , 4 , 3 , 3 , "[Ne]"           ] , 
            'Fe' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.32  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p6 4s2 3d6"    ,   "GTH-PBE-q16"   , 4 , 3 , 3 , "[Ne]"           ] , 
            'Co' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.26  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p6 4s2 3d7"    ,   "GTH-PBE-q17"   , 4 , 3 , 3 , "[Ne]"           ] , 
            'Ni' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.24  /0.5291*rc_factor   ,    conf_n_d    ,   "3s2 3p6 4s2 3d8"    ,   "GTH-PBE-q18"   , 4 , 3 , 3 , "[Ne]"           ] , 
            'Cu' : [    1    ,    0    ,    1    ,    0    ,    conf_a_d   , 1.32  /0.5291*rc_factor   ,    conf_n_d    ,   "4s1 3d10"           ,   "GTH-PBE-q11"   , 4 , 0 , 3 , "[Ar]"           ] , 
            'Zn' : [    1    ,    0    ,    1    ,    0    ,    conf_a_d   , 1.22  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 3d10"           ,   "GTH-PBE-q12"   , 4 , 0 , 3 , "[Ar]"           ] , 
            'Ga' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.22  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p1"            ,   "GTH-PBE-q3"    , 4 , 4 , 0 , "[Ar] 3d10"      ] , 
            'Ge' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.20  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p2"            ,   "GTH-PBE-q4"    , 4 , 4 , 0 , "[Ar] 3d10"      ] , 
            'As' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.19  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p3"            ,   "GTH-PBE-q5"    , 4 , 4 , 0 , "[Ar] 3d10"      ] , 
            'Se' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.20  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p4"            ,   "GTH-PBE-q6"    , 4 , 4 , 0 , "[Ar] 3d10"      ] , 
            'Br' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.20  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p5"            ,   "GTH-PBE-q7"    , 4 , 4 , 0 , "[Ar] 3d10"      ] , 
            'Kr' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.16  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p6     "       ,   "GTH-PBE-q8"    , 4 , 4 , 0 , "[Ar] 3d10"      ] , 
            'Rb' : [    2    ,    1    ,    0    ,    0    ,    conf_a_d   , 2.20  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p6 5s1"        ,   "GTH-PBE-q9 "   , 5 , 4 , 0 , "[Ar] 3d10"      ] , 
            'Sr' : [    2    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.95  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p6 5s2"        ,   "GTH-PBE-q10"   , 5 , 4 , 0 , "[Ar] 3d10"      ] , 
            'Y'  : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.90  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p6 5s2 4d1"    ,   "GTH-PBE-q11"   , 5 , 4 , 4 , "[Ar] 3d10"      ] , 
            'Zr' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.75  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p6 5s2 4d2"    ,   "GTH-PBE-q12"   , 5 , 4 , 4 , "[Ar] 3d10"      ] , 
            'Nb' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.64  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p6 5s2 4d3"    ,   "GTH-PBE-q13"   , 5 , 4 , 4 , "[Ar] 3d10"      ] , 
            'Mo' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.54  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p6 5s1 4d5"    ,   "GTH-PBE-q14"   , 5 , 4 , 4 , "[Ar] 3d10"      ] , 
            'Tc' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.47  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p6 5s2 4d5"    ,   "GTH-PBE-q15"   , 5 , 4 , 4 , "[Ar] 3d10"      ] , 
            'Ru' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.46  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p6 5s2 4d6"    ,   "GTH-PBE-q16"   , 5 , 4 , 4 , "[Ar] 3d10"      ] , 
            'Rh' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.42  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p6 5s2 4d7"    ,   "GTH-PBE-q17"   , 5 , 4 , 4 , "[Ar] 3d10"      ] , 
            'Pd' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.39  /0.5291*rc_factor   ,    conf_n_d    ,   "4s2 4p6 5s2 4d8"    ,   "GTH-PBE-q18"   , 5 , 4 , 4 , "[Ar] 3d10"      ] , 
            'Ag' : [    1    ,    0    ,    1    ,    0    ,    conf_a_d   , 1.45  /0.5291*rc_factor   ,    conf_n_d    ,   "5s1 4d10"           ,   "GTH-PBE-q11"   , 5 , 0 , 4 , "[Kr]"           ] , 
            'Cd' : [    1    ,    0    ,    1    ,    0    ,    conf_a_d   , 1.44  /0.5291*rc_factor   ,    conf_n_d    ,   "5s2 4d10"           ,   "GTH-PBE-q12"   , 5 , 0 , 4 , "[Kr]"           ] , 
            'In' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.42  /0.5291*rc_factor   ,    conf_n_d    ,   "5s2 5p1"            ,   "GTH-PBE-q3"    , 5 , 5 , 0 , "[Kr] 4d10"      ] , 
            'Sn' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.39  /0.5291*rc_factor   ,    conf_n_d    ,   "5s2 5p2"            ,   "GTH-PBE-q4"    , 5 , 5 , 0 , "[Kr] 4d10"      ] , 
            'Sb' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.39  /0.5291*rc_factor   ,    conf_n_d    ,   "5s2 5p3"            ,   "GTH-PBE-q5"    , 5 , 5 , 0 , "[Kr] 4d10"      ] , 
            'Te' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.38  /0.5291*rc_factor   ,    conf_n_d    ,   "5s2 5p4"            ,   "GTH-PBE-q6"    , 5 , 5 , 0 , "[Kr] 4d10"      ] , 
            'I'  : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.39  /0.5291*rc_factor   ,    conf_n_d    ,   "5s2 5p5"            ,   "GTH-PBE-q7"    , 5 , 5 , 0 , "[Kr] 4d10"      ] , 
            'Xe' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.40  /0.5291*rc_factor   ,    conf_n_d    ,   "5s2 5p6"            ,   "GTH-PBE-q8"    , 5 , 5 , 0 , "[Kr] 4d10"      ] , 
            'Cs' : [    2    ,    1    ,    0    ,    0    ,    conf_a_d   , 2.44  /0.5291*rc_factor   ,    conf_n_d    ,   "5s2 5p6 6s1"        ,   "GTH-PBE-q9"    , 6 , 5 , 0 , "[Kr] 4d10"      ] , 
            'Ba' : [    2    ,    1    ,    0    ,    0    ,    conf_a_d   , 2.15  /0.5291*rc_factor   ,    conf_n_d    ,   "5s2 5p6 6s2"        ,   "GTH-PBE-q10"   , 6 , 5 , 0 , "[Kr] 4d10"      ] , 
            'La' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 2.07  /0.5291*rc_factor   ,    conf_n_d    ,   "5s2 5p6 6s2 5d1"    ,   "GTH-PBE-q11"   , 6 , 5 , 5 , "[Kr] 4d10"      ] , 
            'Hf' : [    2    ,    1    ,    1    ,    0    ,    conf_a_d   , 1.75  /0.5291*rc_factor   ,    conf_n_d    ,   "5s2 5p6 6s2 5d2"    ,   "GTH-PBE-q12"   , 6 , 5 , 5 , "[Kr] 4d10 4f14" ] , 
            'Ta' : [    1    ,    0    ,    1    ,    0    ,    conf_a_d   , 1.70  /0.5291*rc_factor   ,    conf_n_d    ,   "6s2 5d3"            ,   "GTH-PBE-q5"    , 6 , 0 , 5 , "[Xe] 4f14"      ] , 
            'W'  : [    1    ,    0    ,    1    ,    0    ,    conf_a_d   , 1.62  /0.5291*rc_factor   ,    conf_n_d    ,   "6s1 5d5"            ,   "GTH-PBE-q6"    , 6 , 0 , 5 , "[Xe] 4f14"      ] , 
            'Re' : [    1    ,    0    ,    1    ,    0    ,    conf_a_d   , 1.51  /0.5291*rc_factor   ,    conf_n_d    ,   "6s2 5d5"            ,   "GTH-PBE-q7"    , 6 , 0 , 5 , "[Xe] 4f14"      ] , 
            'Os' : [    1    ,    0    ,    1    ,    0    ,    conf_a_d   , 1.44  /0.5291*rc_factor   ,    conf_n_d    ,   "6s2 5d6"            ,   "GTH-PBE-q8"    , 6 , 0 , 5 , "[Xe] 4f14"      ] , 
            'Ir' : [    1    ,    0    ,    1    ,    0    ,    conf_a_d   , 1.41  /0.5291*rc_factor   ,    conf_n_d    ,   "6s2 5d7"            ,   "GTH-PBE-q9"    , 6 , 0 , 5 , "[Xe] 4f14"      ] , 
            'Pt' : [    1    ,    0    ,    1    ,    0    ,    conf_a_d   , 1.36  /0.5291*rc_factor   ,    conf_n_d    ,   "6s2 5d8"            ,   "GTH-PBE-q10"   , 6 , 0 , 5 , "[Xe] 4f14"      ] , 
            'Au' : [    1    ,    0    ,    1    ,    0    ,    conf_a_d   , 1.36  /0.5291*rc_factor   ,    conf_n_d    ,   "6s1 5d10"           ,   "GTH-PBE-q11"   , 6 , 0 , 5 , "[Xe] 4f14"      ] , 
            'Hg' : [    1    ,    0    ,    1    ,    0    ,    conf_a_d   , 1.32  /0.5291*rc_factor   ,    conf_n_d    ,   "6s2 5d10"           ,   "GTH-PBE-q12"   , 6 , 0 , 5 , "[Xe] 4f14"      ] , 
            'Tl' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.45  /0.5291*rc_factor   ,    conf_n_d    ,   "6s2 6p1"            ,   "GTH-PBE-q3"    , 6 , 6 , 0 , "[Xe] 4f14 5d10" ] , 
            'Pb' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.46  /0.5291*rc_factor   ,    conf_n_d    ,   "6s2 6p2"            ,   "GTH-PBE-q4"    , 6 , 6 , 0 , "[Xe] 4f14 5d10" ] , 
            'Bi' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.48  /0.5291*rc_factor   ,    conf_n_d    ,   "6s2 6p3"            ,   "GTH-PBE-q5"    , 6 , 6 , 0 , "[Xe] 4f14 5d10" ] , 
            'Po' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.40  /0.5291*rc_factor   ,    conf_n_d    ,   "6s2 6p4"            ,   "GTH-PBE-q6"    , 6 , 6 , 0 , "[Xe] 4f14 5d10" ] , 
            'At' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.50  /0.5291*rc_factor   ,    conf_n_d    ,   "6s2 6p5"            ,   "GTH-PBE-q7"    , 6 , 6 , 0 , "[Xe] 4f14 5d10" ] , 
            'Rn' : [    1    ,    1    ,    0    ,    0    ,    conf_a_d   , 1.50  /0.5291*rc_factor   ,    conf_n_d    ,   "6s2 6p6"            ,   "GTH-PBE-q8"    , 6 , 6 , 0 , "[Xe] 4f14 5d10" ] , 
            }



def main(args):
    elements = args

    workdir_prefix = "/data/ralph/deltatests/newpseudos/"
    cp2k_command   = "module load gcc-suite/5.3.0; /users/ralph/cp2k/cp2k-code/exe/Linux-x86-64-gfortran-local/cp2k.popt"

    workdir_prefix += randomword(6)

    for el in elements:
        #retrieve a default pseudopotential for that element as a starting guess
        myguess = requests.get(PSEUDOPOTENTIAL_URL, data={'family':'GTH-PBE', 'element':el}, verify=False).json()[el]

        #prepare the input file
        settings = { 'element' : el, 
                     'elconf'  : (atoms_db[el][12] if atoms_db[el][12]!='none' else '') + ' ' + atoms_db[el][7],
                     'core'    : atoms_db[el][12], 
                     'maxang'  : 2,
                     'xc'      : 'PBE',
                     'pseudoguess':  myguess,
                     'family'  : 'self-optimized-pbe',
                   }
        myinput = input_template.format(**settings)

        
        #prepare the working directory
        workdir = os.path.join(workdir_prefix, "{:s}".format(el))

        try:
            os.makedirs(workdir)
        except OSError:
            pass

        os.chdir(workdir)

        #write the input to a file
        of = open("job.inp", "w")
        of.write(myinput)
        of.close()
               
        #and run it.
        os.popen(cp2k_command + " -i job.inp -o job.out")

        #open the resulting pseudopotential data
        infile = open("GTH-PARAMETER")
        pseudo_data = "".join(infile.readlines()[1:])


        #upload the PP data.
        req = requests.post(PSEUDOPOTENTIAL_URL, data={'family':settings['family'], 'element':el, 'pseudo': pseudo_data, 'overwrite': True}, verify=False)
        req.raise_for_status()

        settings['pseudoguess'] = pseudo_data
        settings['elconf']      = 'CORE' + atoms_db[el][7]

        myinput = upf_template.format(**settings)
        of = open("upf.inp", "w")
        of.write(myinput)
        of.close()
               
        #and run it.
        os.popen(cp2k_command + " -i upf.inp -o upf.out")

        #open the resulting UPF file
        infile = open("output.UPF")
        upf_data = infile.read()

        req = requests.post(PSEUDOPOTENTIAL_URL, data={'family':settings['family']+'-UPF', 'element':el, 'pseudo': upf_data, 'overwrite': True}, verify=False)
        req.raise_for_status()
        upf_id = req.json()['id']


        #now, create a 'method' with this PP.
        req = requests.post(METHOD_URL, data={'code':'espresso', 
                                              'basis_set': 55, 
                                              'pseudopotential': upf_id, 
                                              'settings': json.dumps({'smearing': 'marzari-vanderbilt', 'xc': 'PBE', 'cutoff_pw': 3401.4244569396633, 'sigma': 0.027211395655517306})}, 
                            verify=False)
        req.raise_for_status()
        method_id = req.json()['id']

        #create the corresponding delta tasks for this method.
        req = requests.post(TASK_URL, data={  'method':method_id, 
                                              'test': 'deltatest_'+el }, 
                            verify=False)
        req.raise_for_status()
        tasklist = req.json()


        print "{:2s}: {:}".format(el, ", ".join([str(x) for x in tasklist]))

def randomword(length):
    import random, string
    return ''.join(random.choice(string.lowercase + string.digits) for i in range(length))

if __name__ == "__main__":
    main(argv[1:])
