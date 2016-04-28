#!/usr/bin/env python
from sys import argv
from ase.calculators.cp2k import CP2K
from ase import Atoms
import requests, os, json, hashlib
from base64 import urlsafe_b64encode
from jinja2 import Template

SERVER = 'https://172.23.64.223'
SERVER = 'http://127.0.0.1:5001'
PSEUDOPOTENTIAL_URL   = SERVER + '/pseudo'
METHOD_URL            = SERVER + '/method'
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



            
atoms_db ={  
            'H'  : [   "/users/ralph/work/fatman/PBE/H/q1/H.inp"  , "1s1"                ,  "none"          ] , 
            'He' : [   "/users/ralph/work/fatman/PBE/He/q2/He.inp", "1s2"                ,  "none"          ] , 
            'Li' : [   "/users/ralph/work/fatman/PBE/Li/q3/Li.inp", "1s2 2s1"            ,  "none"          ] , 
            'Be' : [   "/users/ralph/work/fatman/PBE/Be/q4/Be.inp", "1s2 2s2"            ,  "none"          ] , 
            'B'  : [   "/users/ralph/work/fatman/PBE/B/q3/B.inp"  , "2s2 2p1"            ,  "[He]"          ] , 
            'C'  : [   "/users/ralph/work/fatman/PBE/C/q4/C.inp"  , "2s2 2p2"            ,  "[He]"          ] , 
            'N'  : [   "/users/ralph/work/fatman/PBE/N/q5/N.inp"  , "2s2 2p3"            ,  "[He]"          ] , 
            'O'  : [   "/users/ralph/work/fatman/PBE/O/q6/O.inp"  , "2s2 2p4"            ,  "[He]"          ] , 
            'F'  : [   "/users/ralph/work/fatman/PBE/F/q7/F.inp"  , "2s2 2p5"            ,  "[He]"          ] , 
            'Ne' : [   "/users/ralph/work/fatman/PBE/Ne/q8/Ne.inp", "2s2 2p6"            ,  "[He]"          ] , 
            'Na' : [   "/users/ralph/work/fatman/PBE/H/q1/H.inp"  ,        ] , 
            'Mg' : [   "/users/ralph/work/fatman/PBE/H/q1/H.inp"  ,        ] , 
            'Al' : [   "/users/ralph/work/fatman/PBE/H/q1/H.inp"  ,        ] , 
            'Si' : [   "/users/ralph/work/fatman/PBE/H/q1/H.inp"  ,        ] , 
            'P'  : [    ] , 
            'S'  : [    ] , 
            'Cl' : [    ] , 
            'Ar' : [    ] , 
            'K'  : [    ] , 
            'Ca' : [    ] , 
            'Sc' : [    ] , 
            'Ti' : [    ] , 
            'V'  : [    ] , 
            'Cr' : [    ] , 
            'Mn' : [    ] , 
            'Fe' : [    ] , 
            'Co' : [    ] , 
            'Ni' : [    ] , 
            'Cu' : [    ] , 
            'Zn' : [    ] , 
            'Ga' : [    ] , 
            'Ge' : [    ] , 
            'As' : [    ] , 
            'Se' : [    ] , 
            'Br' : [    ] , 
            'Kr' : [    ] , 
            'Rb' : [    ] , 
            'Sr' : [    ] , 
            'Y'  : [    ] , 
            'Zr' : [    ] , 
            'Nb' : [    ] , 
            'Mo' : [    ] , 
            'Tc' : [    ] , 
            'Ru' : [    ] , 
            'Rh' : [    ] , 
            'Pd' : [    ] , 
            'Ag' : [    ] , 
            'Cd' : [    ] , 
            'In' : [    ] , 
            'Sn' : [    ] , 
            'Sb' : [    ] , 
            'Te' : [    ] , 
            'I'  : [    ] , 
            'Xe' : [    ] , 
            'Cs' : [    ] , 
            'Ba' : [    ] , 
            'La' : [    ] , 
            'Hf' : [    ] , 
            'Ta' : [    ] , 
            'W'  : [    ] , 
            'Re' : [    ] , 
            'Os' : [    ] , 
            'Ir' : [    ] , 
            'Pt' : [    ] , 
            'Au' : [    ] , 
            'Hg' : [    ] , 
            'Tl' : [    ] , 
            'Pb' : [    ] , 
            'Bi' : [    ] , 
            'Po' : [    ] , 
            'At' : [    ] , 
            'Rn' : [    ] , 
            }



def main(args):
    elements = args

    workdir_prefix = "/data/ralph/deltatests/newpseudos/"
    cp2k_command   = "module load gcc-suite/5.3.0; /users/ralph/cp2k/cp2k-code/exe/Linux-x86-64-gfortran-local/cp2k.popt"

    workdir_prefix += randomword(6)
    q = hashlib.sha1( ''.join(input_template_str.split()) + ''.join(str(atoms_db).split()) )
    pp_unique_name = 'family_28_04_16' #urlsafe_b64encode(q.digest())[:7]

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
       #myinput = input_template.render(settings)
        with open(atoms_db[el][0]) as infile:
            myinput_lines = infile.readlines()
            myinput = "".join(myinput_lines)
        
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
        pseudo_data = "".join(infile.readlines()[1:] + ['# ' + x for x in myinput_lines])

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
                                              'settings': json.dumps({'smearing': 'marzari-vanderbilt', 'xc': 'PBE', 'cutoff_pw': 3401.4244569396633, 'sigma': 0.027211395655517306, 'cp2k_template': input_template_str})}, 
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
