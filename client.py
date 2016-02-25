#!/usr/bin/env python

import requests, os
from tools import Json2Atoms
from codehandling import HandlerFactory

SERVER = 'http://172.23.64.223'
TASKS_URL = SERVER + '/fatman/tasks'
RESULTS_URL = SERVER + '/fatman/results'
BASISSET_URL = SERVER + '/fatman/basis'
PSEUDOPOTENTIAL_URL   = SERVER + '/fatman/pseudo'

def main():
    """The fatman client queries the DB for new tasks and runs them until none are left"""

    while 1:
        #request a task -- and set its status to "pending"
        req = requests.get(TASKS_URL, params={'limit': 1, 'status': 'new'})
        req.raise_for_status()

        tasks = req.json()
        if len(tasks) == 0:
            print('No task to run, exiting')
            return 0

        req = requests.patch(SERVER + tasks[0]['_links']['self'], data={'status': 'pending', 'machine':os.uname()[1]})
        req.raise_for_status()
        task = req.json()

        #get structure and convert to ASE object
        struct_json = task['structure']['ase_structure']
        struct = Json2Atoms(struct_json)

        #which code to use with which settings?
        mymethod = task['method']
        if "kpoints" in struct.info["key_value_pairs"].keys():
            mymethod["kpoints"] = struct.info["key_value_pairs"]["kpoints"]
            #kindof a hack: kpoints are specified with each structure, but are in fact code parameters
        basis = requests.get(BASISSET_URL, data={'family':mymethod['basis_set'], 'element':set(struct.get_chemical_symbols())}).json()
        pseudo = requests.get(PSEUDOPOTENTIAL_URL, data={'family':mymethod['pseudopotential'], 'element':set(struct.get_chemical_symbols())}).json()

        print basis 
        if not "kind_settings" in mymethod.keys():
            mymethod["kind_settings"] = {}

        for x in set(basis.keys()) & set(pseudo.keys()):
            mymethod["kind_settings"][x] =  {"basis_set":basis[x], "potential":pseudo[x]}

        #create our code-running object with the relevant settings.
        codehandler = HandlerFactory(struct, mymethod)

        #update the task's status to "running"
        req = requests.patch(SERVER + task['_links']['self'], data={'status': 'running', 'machine':os.uname()[1]})
        req.raise_for_status()
        task = req.json()

        #run the code
        try:
            e,output_file_path = codehandler.runOne()
        except RuntimeError or AssertionError:
            req = requests.patch(SERVER + task['_links']['self'], data={'status': 'error'})
            req.raise_for_status()


        #throw the energy into the DB, upload the output file
        req = requests.post(RESULTS_URL, data={'energy': e, 'task_id': task['id']})
        req.raise_for_status()
        
        #upload the output file
       #with open(output_file_path) as f:
       #    req = requests.post(RESULTS_URL, data=f)


        #update the task's status to "done"
        req = requests.patch(SERVER + task['_links']['self'], data={'status': 'done'})
        req.raise_for_status()
        task = req.json()

    #upload the data


if __name__=="__main__":
    main()
