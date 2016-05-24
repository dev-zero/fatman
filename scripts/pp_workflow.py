#!/usr/bin/env python
from sys import argv
from ase.calculators.cp2k import CP2K
from ase import Atoms
import requests, os, json, hashlib
from base64 import urlsafe_b64encode
from jinja2 import Template
from fatman.tools import randomword

SERVER = 'https://172.23.64.223'
SERVER = 'http://127.0.0.1:5001'
PSEUDOPOTENTIAL_URL   = SERVER + '/pseudo'
METHOD_URL            = SERVER + '/methods'
TASK_URL              = SERVER + '/tasks/0'


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

conf_a_d  = 1.0
rc_factor = 3.00
conf_n_d  = 10


atoms_db_largecore ={  
            'Li' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Li/q1/GTH-PARAMETER"   , "2s1"                ,  "[He]"          ] , 
            'Be' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Be/q2/GTH-PARAMETER"   , "2s2"                ,  "[He]"          ] , 
            'Na' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Na/q1/GTH-PARAMETER"   , "3s1"                ,  "[Ne]"          ] , 
            'Mg' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Mg/q2/GTH-PARAMETER"   , "3s2"                ,  "[Ne]"          ] , 
            'K'  : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/K/q1/GTH-PARAMETER"    , "4s1"                ,  "[Ar]"          ] ,
            'Ca' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Ca/q2/GTH-PARAMETER"   , "4s2"                ,  "[Ar]"          ] ,
            'Sc' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Sc/q3/GTH-PARAMETER"   , "4s2 3d1"            ,  "[Ar]"          ] ,
            'Ti' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Ti/q4/GTH-PARAMETER"   , "4s2 3d2"            ,  "[Ar]"          ] ,
            'V'  : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/V/q5/GTH-PARAMETER"    , "4s2 3d3"            ,  "[Ar]"          ] ,
            'Cr' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Cr/q6/GTH-PARAMETER"   , "4s2 3d4"            ,  "[Ar]"          ] ,
            'Mn' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Mn/q7/GTH-PARAMETER"   , "4s2 3d5"            ,  "[Ar]"          ] ,
            'Fe' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Fe/q8/GTH-PARAMETER"   , "4s2 3d6"            ,  "[Ar]"          ] ,
            'Co' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Co/q9/GTH-PARAMETER"   , "4s2 3d7"            ,  "[Ar]"          ] ,
            'Ni' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Ni/q10/GTH-PARAMETER"  , "4s2 3d8"            ,  "[Ar]"          ] ,
            'Cu' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Cu/q1/GTH-PARAMETER"   , "4s1"                ,  "[Ar] 3d10"     ] ,
            'Zn' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Zn/q2/GTH-PARAMETER"   , "4s2"                ,  "[Ar] 3d10"     ] ,
            'Ga' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Ga/q3/GTH-PARAMETER"   , "4s2 4p1"            ,  "[Ar] 3d10"     ] }

atoms_db_mediumcore ={  
            'Cu' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Cu/q11/GTH-PARAMETER"  , "4s1 3d10"           ,  "[Ar]"          ] ,
            'Zn' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Zn/q12/GTH-PARAMETER"  , "4s2 3d10"           ,  "[Ar]"          ] }

            
atoms_db_smallcore ={  
            'H'  : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/H/q1/GTH-PARAMETER"    , "1s1"                ,  "none"          ] , 
            'He' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/He/q2/GTH-PARAMETER"   , "1s2"                ,  "none"          ] , 
            'Li' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Li/q3/GTH-PARAMETER"   , "1s2 2s1"            ,  "none"          ] , 
            'Be' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Be/q4/GTH-PARAMETER"   , "1s2 2s2"            ,  "none"          ] , 
            'B'  : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/B/q3/GTH-PARAMETER"    , "2s2 2p1"            ,  "[He]"          ] , 
            'C'  : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/C/q4/GTH-PARAMETER"    , "2s2 2p2"            ,  "[He]"          ] , 
            'N'  : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/N/q5/GTH-PARAMETER"    , "2s2 2p3"            ,  "[He]"          ] , 
            'O'  : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/O/q6/GTH-PARAMETER"    , "2s2 2p4"            ,  "[He]"          ] , 
            'F'  : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/F/q7/GTH-PARAMETER"    , "2s2 2p5"            ,  "[He]"          ] , 
            'Ne' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Ne/q8/GTH-PARAMETER"   , "2s2 2p6"            ,  "[He]"          ] , 
            'Na' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Na/q9/GTH-PARAMETER"   , "2s2 2p6 3s1"        ,  "[He]"          ] , 
            'Mg' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Mg/q10/GTH-PARAMETER"  , "2s2 2p6 3s2"        ,  "[He]"          ] , 
            'Al' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Al/q3/GTH-PARAMETER"   , "3s2 3p1"            ,  "[Ne]"          ] , 
            'Si' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Si/q4/GTH-PARAMETER"   , "3s2 3p2"            ,  "[Ne]"          ] , 
            'P'  : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/P/q5/GTH-PARAMETER"    , "3s2 3p3"            ,  "[Ne]"          ] , 
            'S'  : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/S/q6/GTH-PARAMETER"    , "3s2 3p4"            ,  "[Ne]"          ] , 
            'Cl' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Cl/q7/GTH-PARAMETER"   , "3s2 3p5"            ,  "[Ne]"          ] , 
            'Ar' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Ar/q8/GTH-PARAMETER"   , "3s2 3p6"            ,  "[Ne]"          ] , 
            'K'  : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/K/q9/GTH-PARAMETER"    , "3s2 3p6 4s1"        ,  "[Ne]"          ] ,
            'Ca' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Ca/q10/GTH-PARAMETER"  , "3s2 3p6 4s2"        ,  "[Ne]"          ] ,
            'Sc' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Sc/q11/GTH-PARAMETER"  , "3s2 3p6 4s2 3d1"    ,  "[Ne]"          ] ,
            'Ti' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Ti/q12/GTH-PARAMETER"  , "3s2 3p6 4s2 3d2"    ,  "[Ne]"          ] ,
            'V'  : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/V/q13/GTH-PARAMETER"   , "3s2 3p6 4s2 3d3"    ,  "[Ne]"          ] ,
            'Cr' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Cr/q14/GTH-PARAMETER"  , "3s2 3p6 4s2 3d4"    ,  "[Ne]"          ] ,
            'Mn' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Mn/q15/GTH-PARAMETER"  , "3s2 3p6 4s2 3d5"    ,  "[Ne]"          ] ,
            'Fe' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Fe/q16/GTH-PARAMETER"  , "3s2 3p6 4s2 3d6"    ,  "[Ne]"          ] ,
            'Co' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Co/q17/GTH-PARAMETER"  , "3s2 3p6 4s2 3d7"    ,  "[Ne]"          ] ,
            'Ni' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Ni/q18/GTH-PARAMETER"  , "3s2 3p6 4s2 3d8"    ,  "[Ne]"          ] ,
            'Cu' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Cu/q19/GTH-PARAMETER"  , "3s2 3p6 4s1 3d10"   ,  "[Ne]"          ] ,
            'Zn' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Zn/q20/GTH-PARAMETER"  , "3s2 3p6 4s2 3d10"   ,  "[Ne]"          ] ,
            'Ga' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Ga/q13/GTH-PARAMETER"  , "3d10 4s2 4p1"       ,  "[Ar]"          ] ,
            'Ge' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Ge/q4/GTH-PARAMETER"   , "4s2 4p2"            ,  "[Ar] 3d10"     ] ,
            'As' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/As/q5/GTH-PARAMETER"   , "4s2 4p3"            ,  "[Ar] 3d10"     ] ,
            'Se' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Se/q6/GTH-PARAMETER"   , "4s2 4p4"            ,  "[Ar] 3d10"     ] ,
            'Br' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Br/q7/GTH-PARAMETER"   , "4s2 4p5"            ,  "[Ar] 3d10"     ] ,
            'Kr' : [   "/users/ralph/work/fatman/PBE_23May2016/PBE/Kr/q8/GTH-PARAMETER"   , "4s2 4p6"            ,  "[Ar] 3d10"     ] ,
            'Rb' : [   ""                                                                , "4s2 4p6 5s1"        ,  "[Ar] 3d10"     ] ,
            'Sr' : [   ""                                                                , "4s2 4p6 5s2"        ,  "[Ar] 3d10"     ] ,
            'Y'  : [   ""                                                                , "4s2 4p6 5s2 4d1"    ,  "[Ar] 3d10"     ] ,
            'Zr' : [   ""                                                                , "4s2 4p6 5s2 4d2"    ,  "[Ar] 3d10"     ] ,
            'Nb' : [   ""                                                                , "4s2 4p6 5s2 4d3"    ,  "[Ar] 3d10"     ] ,
            'Mo' : [   ""                                                                , "4s2 4p6 5s1 4d5"    ,  "[Ar] 3d10"     ] ,
            'Tc' : [   ""                                                                , "4s2 4p6 5s2 4d5"    ,  "[Ar] 3d10"     ] ,
            'Ru' : [   ""                                                                , "4s2 4p6 5s2 4d6"    ,  "[Ar] 3d10"     ] ,
            'Rh' : [   ""                                                                , "4s2 4p6 5s2 4d7"    ,  "[Ar] 3d10"     ] ,
            'Pd' : [   ""                                                                , "4s2 4p6 5s2 4d8"    ,  "[Ar] 3d10"     ] ,
            'Ag' : [   ""                                                                , "5s1 4d10"           ,  "[Kr]"          ] ,
            'Cd' : [   ""                                                                , "5s2 4d10"           ,  "[Kr]"          ] ,
            'In' : [   ""                                                                , "5s2 5p1"            ,  "[Kr] 4d10"     ] ,
            'Sn' : [   ""                                                                , "5s2 5p2"            ,  "[Kr] 4d10"     ] ,
            'Sb' : [   ""                                                                , "5s2 5p3"            ,  "[Kr] 4d10"     ] ,
            'Te' : [   ""                                                                , "5s2 5p4"            ,  "[Kr] 4d10"     ] ,
            'I'  : [   ""                                                                , "5s2 5p5"            ,  "[Kr] 4d10"     ] ,
            'Xe' : [   ""                                                                , "5s2 5p6"            ,  "[Kr] 4d10"     ] ,
            'Cs' : [   ""                                                                , "5s2 5p6 6s1"        ,  "[Kr] 4d10"     ] ,
            'Ba' : [   ""                                                                , "5s2 5p6 6s2"        ,  "[Kr] 4d10"     ] ,
            'La' : [   ""                                                                , "5s2 5p6 6s2 5d1"    ,  "[Kr] 4d10"     ] ,
            'Hf' : [   ""                                                                , "5s2 5p6 6s2 5d2"    ,  "[Kr] 4d10 4f14"] ,
            'Ta' : [   ""                                                                , "6s2 5d3"            ,  "[Xe] 4f14"     ] ,
            'W'  : [   ""                                                                , "6s1 5d5"            ,  "[Xe] 4f14"     ] ,
            'Re' : [   ""                                                                , "6s2 5d5"            ,  "[Xe] 4f14"     ] ,
            'Os' : [   ""                                                                , "6s2 5d6"            ,  "[Xe] 4f14"     ] ,
            'Ir' : [   ""                                                                , "6s2 5d7"            ,  "[Xe] 4f14"     ] ,
            'Pt' : [   ""                                                                , "6s2 5d8"            ,  "[Xe] 4f14"     ] ,
            'Au' : [   ""                                                                , "6s1 5d10"           ,  "[Xe] 4f14"     ] ,
            'Hg' : [   ""                                                                , "6s2 5d10"           ,  "[Xe] 4f14"     ] ,
            'Tl' : [   ""                                                                , "6s2 6p1"            ,  "[Xe] 4f14 5d10"] ,
            'Pb' : [   ""                                                                , "6s2 6p2"            ,  "[Xe] 4f14 5d10"] ,
            'Bi' : [   ""                                                                , "6s2 6p3"            ,  "[Xe] 4f14 5d10"] ,
            'Po' : [   ""                                                                , "6s2 6p4"            ,  "[Xe] 4f14 5d10"] ,
            'At' : [   ""                                                                , "6s2 6p5"            ,  "[Xe] 4f14 5d10"] ,
            'Rn' : [   ""                                                                , "6s2 6p6"            ,  "[Xe] 4f14 5d10"] ,
            }



def main(args):
    elements = args
    atoms_db = atoms_db_largecore

    workdir_prefix = "/data/ralph/deltatests/newpseudos/"
    cp2k_command   = "module load gcc-suite/5.3.0; /users/ralph/cp2k/cp2k-code/exe/Linux-x86-64-gfortran-local/cp2k.popt"

    workdir_prefix += randomword(6)
    q = hashlib.sha1( ''.join(input_template_str.split()) + ''.join(str(atoms_db).split()) )
    pp_unique_name = 'family_23_05_16_largecore' #urlsafe_b64encode(q.digest())[:7]

    for el in elements:
        #retrieve a default pseudopotential for that element as a starting guess
       #myguess = requests.get(PSEUDOPOTENTIAL_URL, data={'family':'GTH-PBE', 'element':el}, verify=False).json()[el]

       ##prepare the input file
        settings = { 'element' : el, 
                     'family'  : pp_unique_name,
                     'pseudoguess' : '',
                     'elconf'      : atoms_db[el][1],
                     'core'        : atoms_db[el][2],
                   }
     # #myinput = input_template.render(settings)
     #  with open(atoms_db[el][0]) as infile:
     #      myinput_lines = infile.readlines()
     #      myinput = "".join(myinput_lines)
     #  
     #  #prepare the working directory
        workdir = os.path.join(workdir_prefix, "{:s}".format(el))

        try:
            os.makedirs(workdir)
        except OSError:
            pass

        filename = os.path.dirname(atoms_db[el][0])
        os.chdir(workdir)
     #  os.system("cp {:}/GTH-PARAMETER {:}".format(filename, workdir))


     #  #write the input to a file
     #  of = open("job.inp", "w")
     #  of.write(myinput)
     #  of.close()
               
        #and run it.
        #os.popen(cp2k_command + " -i job.inp -o job.out")
        myinput_lines = []
        #open the resulting pseudopotential data
        #infile = open("GTH-PARAMETER")
        infile = open(atoms_db[el][0])
        #pseudo_data = myguess #"".join(infile.readlines()[1:] + ['# ' + x for x in myinput_lines])
        pseudo_data = "".join(infile.readlines()[1:] )

        #upload the PP data.
        req = requests.post(PSEUDOPOTENTIAL_URL, data={'family':pp_unique_name, 'element':el, 'pseudo': pseudo_data, 'overwrite': True}, verify=False)
        req.raise_for_status()

        settings['pseudoguess'] = pseudo_data

        myinput = upf_template.render(settings)
        of = open("upf.inp", "w")
        of.write(myinput)
        of.close()
               
        #and run it.
        os.popen(cp2k_command + " -i upf.inp -o upf.out")

        #open the resulting UPF file
        infile = open("output.UPF")
        upf_data = infile.read()

        req = requests.post(PSEUDOPOTENTIAL_URL, data={'family':pp_unique_name+'-UPF', 'element':el, 'pseudo': upf_data, 'overwrite': True}, verify=False)
        req.raise_for_status()
        upf_id = req.json()['id']


        #now, create a 'method' with this PP.
        req = requests.post(METHOD_URL, data={'code':'espresso', 
                                              'basis_set': 55, 
                                              'pseudopotential': upf_id, 
                                              'settings': json.dumps({'smearing': 'marzari-vanderbilt', 'xc': 'PBE', 'cutoff_pw': 3401.4244569396633, 'sigma': 0.027211395655517306, 'max_kpoints':[20,20,20]})}, 
                            verify=False)
        req.raise_for_status()
        method_id = req.json()['id']

        #create the corresponding delta tasks for this method.
        req = requests.post(TASK_URL, data={  'method':method_id, 
                                              'test': 'deltatest_'+el,
                                              'priority': 1}, 
                            verify=False)
        req.raise_for_status()
        tasklist = req.json()

        print("{:2s}: {:}".format(el, ", ".join([str(x) for x in tasklist])))

if __name__ == "__main__":
    main(argv[1:])
