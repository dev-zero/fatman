#!/usr/bin/env python

import requests, os
from fatman.tools import Json2Atoms
from codehandling import HandlerFactory
from time import sleep
from datetime import datetime
import os, sys

SERVER = 'https://172.23.64.223'
TASKS_URL = SERVER + '/fatman/tasks'
RESULTS_URL = SERVER + '/fatman/results'
BASISSET_URL = SERVER + '/fatman/basis'
PSEUDOPOTENTIAL_URL   = SERVER + '/fatman/pseudo'

def main(maxtask=1):
    """The fatman client queries the DB for new tasks and runs them until none are left"""
    ntask=0

    while ntask<maxtask:
        ntask+=1
        os.system ("rm /tmp/espresso-remote-task/*")

        #request a task
        req = requests.get(TASKS_URL, params={'limit': 1, 'status': 'new'}, verify=False)
        req.raise_for_status()

        #perhaps this could be replaced by daemon-like behavior: sleep, and check for new tasks after a few minutes
        tasks = req.json()
        if len(tasks) == 0:
            print '{:}: currently no more tasks. Goodbye.'.format(datetime.now())
            quit()

        #set the task's status to pending
        req = requests.patch(SERVER + tasks[0]['_links']['self'], data={'status': 'new', 'machine':'hpc'}, verify=False)
        req.raise_for_status()
        task = req.json()

        try: 
            #get structure and convert to ASE object
            struct_json = task['structure']['ase_structure']
            struct = Json2Atoms(struct_json)

            #which code to use with which settings?
            mymethod = task['method']
            if "kpoints" in struct.info["key_value_pairs"].keys():
                mymethod["kpoints"] = struct.info["key_value_pairs"]["kpoints"]
                #kindof a hack: kpoints are specified with each structure, but are in fact code parameters

            if "charge" in struct.info["key_value_pairs"].keys():
                mymethod["charge"] = struct.info["key_value_pairs"]["charge"]

            if "multiplicity" in struct.info["key_value_pairs"].keys():
                mymethod["multiplicity"] = struct.info["key_value_pairs"]["multiplicity"]

            #REMOTE: get dicts containing the Pseudos and Basissets for all required chemical elements
            pseudo = requests.get(PSEUDOPOTENTIAL_URL, data={'family':mymethod['pseudopotential'], 'element':set(struct.get_chemical_symbols())}, verify=False).json()

            if mymethod["code"] == "cp2k":
                basis = requests.get(BASISSET_URL, data={'family':mymethod['basis_set'], 'element':set(struct.get_chemical_symbols())}, verify=False).json()
                if "kind_settings" in mymethod["settings"].keys():
                    ks = mymethod["settings"]["kind_settings"].copy()
                else:
                    ks = {}

                mymethod["kind_settings"] = {}


                for x in set(basis.keys()) & set(pseudo.keys()):
                    mymethod["kind_settings"][x] =  {"basis_set":basis[x], "potential":pseudo[x]}
                    for k,v in ks.items():
                        mymethod["kind_settings"][x][k] = v 

            elif mymethod["code"] == "espresso":
                odir = os.path.join("/users/ralph/work/espresso/", mymethod['pseudopotential'])
                try:
                    os.makedirs(odir)
                except OSError:
                    pass

                pp_files=  []
                for e, pp in pseudo.items():
                    of = open(os.path.join(odir, e+".UPF"), 'w')
                    of.write (pp)
                    of.close()

                    os.system("cp {:}/{:}.UPF /tmp/espresso-remote-task".format(odir,e))


            #create our code-running object with the relevant settings.
            codehandler = HandlerFactory(struct, mymethod)

            #REMOTE: update the task's status to "running"
            req = requests.patch(SERVER + task['_links']['self'], data={'status': 'running-remote'}, verify=False)
            req.raise_for_status()
            task = req.json()

            #run the code
            output_file_path = codehandler.createOne()

            
            runscript_template = (open("./dora_script_espresso.txt")).read()
            params = {"workdir"    : "{:04d}/{:s}".format(mymethod["id"],struct.info["key_value_pairs"]["identifier"]), 
                      "taskupdate" : task['_links']['self'], 
                      "natom"      : len(struct.get_masses()),
                      "id"         : task['id']}

            of = open(os.path.join("/tmp/espresso-remote-task","runjob_dora.sh"),"w")
            of.write(runscript_template.format(**params))
            of.close()

            os.system("ssh rkoitz@ela.cscs.ch 'mkdir -p /users/rkoitz/deltatests/espresso/{workdir}/'".format(**params))
            os.system("scp -rp /tmp/espresso-remote-task/* rkoitz@ela.cscs.ch:/users/rkoitz/deltatests/espresso/{workdir}/".format(**params))

        except :
            req = requests.patch(SERVER + task['_links']['self'], data={'status': 'error'}, verify=False)
            req.raise_for_status()


if __name__=="__main__":
    main(int(sys.argv[1]))
