#!/usr/bin/env python

import requests, os
from fatman.tools import Json2Atoms, randomword
from codehandling import HandlerFactory
from time import sleep
from datetime import datetime
import os, sys

SERVER = 'https://172.23.64.223'
TASKS_URL = SERVER + '/fatman/tasks'
RESULTS_URL = SERVER + '/fatman/results'
BASISSET_URL = SERVER + '/fatman/basis'
PSEUDOPOTENTIAL_URL   = SERVER + '/fatman/pseudo'

def main():
    """The fatman client queries the DB for new tasks and runs them until none are left"""

    exitword = "./{:}".format(randomword(8))
    print "##########################################################################"
    print "##                       FATMAN CLIENT                                  ##"
    print "##########################################################################"
    print "Run the following command to cleanly shutdown the client after finishing a task: touch {:}".format(exitword)
    print "##########################################################################"
    sys.stdout.flush()

    while 1:
        #check if the 'exit file' exists, then exit
        if os.path.exists(exitword):
            os.remove(exitword)
            return 0

        #request a task
        req = requests.get(TASKS_URL, params={'limit': 1, 'status': 'new'}, verify=False)
        req.raise_for_status()

        #perhaps this could be replaced by daemon-like behavior: sleep, and check for new tasks after a few minutes
        tasks = req.json()
        if len(tasks) == 0:
            print '{:}: currently no tasks. Waiting 5 minutes.'.format(datetime.now())
            sleep (5*60)
            continue

        #set the task's status to pending
        req = requests.patch(SERVER + tasks[0]['_links']['self'], data={'status': 'pending', 'machine':os.uname()[1]}, verify=False)
        req.raise_for_status()
        task = req.json()

        try: 
            #get structure and convert to ASE object
            struct_json = task['structure']['ase_structure']
            struct = Json2Atoms(struct_json)

            #which code to use with which settings?
            mymethod = task['method']
            if "kpoints" in struct.info["key_value_pairs"].keys() and not "kpoints" in mymethod['settings'].keys():
                mymethod["kpoints"] = struct.info["key_value_pairs"]["kpoints"]
                #kindof a hack: kpoints are specified with each structure, but are in fact code parameters
            elif "kpoints" in mymethod['settings'].keys():
                mymethod["kpoints"] = mymethod['settings']['kpoints']
                #kpoints specified in the 'method' override those in the 'structure'

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

                for e, pp in pseudo.items():
                    of = open(os.path.join(odir, e+".UPF"), 'w')
                    of.write (pp)
                    of.close()


            #create our code-running object with the relevant settings.
            codehandler = HandlerFactory(struct, mymethod)

            #REMOTE: update the task's status to "running"
            req = requests.patch(SERVER + task['_links']['self'], data={'status': 'running'}, verify=False)
            req.raise_for_status()
            task = req.json()

            #run the code
            e,output_file_path = codehandler.runOne()

            #REMOTE: throw the energy into the DB, get the id of the stored result
            req = requests.post(RESULTS_URL, data={'energy': e, 'task_id': task['id']}, verify=False)
            req.raise_for_status()
            done_task = req.json()
            
            #REMOTE: upload the output file
            os.system ("bzip2 {}".format(output_file_path))
            with open(output_file_path + ".bz2") as f:
                req = requests.post(SERVER + done_task["_links"]["self"]+'/file', files={'file':f}, verify=False)
                req.raise_for_status()


            #REMOTE: update the task's status to "done"
            req = requests.patch(SERVER + task['_links']['self'], data={'status': 'done'}, verify=False)
            req.raise_for_status()
            task = req.json()

        except: 
            req = requests.patch(SERVER + task['_links']['self'], data={'status': 'error'}, verify=False)
            req.raise_for_status()
            sleep (60)  #just take a short break to avoid cycling through too many tasks if errors happen


if __name__=="__main__":
    main()
