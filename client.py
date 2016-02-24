#!/usr/bin/env python

import requests
from tools import Json2Atoms
from codehandling import HandlerFactory

SERVER = 'http://172.23.64.223'
TASKS_URL = SERVER + '/fatman/tasks'

def main():
    """The fatman client queries the DB for new tasks and runs them until none are left"""

    #request a task -- and set its status to "pending"

    req = requests.get(TASKS_URL, params={'limit': 1, 'status': 'new'})
    req.raise_for_status()

    tasks = req.json()
    if len(tasks) == 0:
        print('No task to run, exiting')
        return 0

    req = requests.patch(SERVER + tasks[0]['_links']['self'], data={'status': 'pending'})
    req.raise_for_status()
    task = req.json()
    print task

    #which structure?
    struct_json = task['structure']['ase_structure']
    struct = Json2Atoms(struct_json)

    #which code to use with which settings?
    mymethod = task['method']
    if "kpoints" in struct.info["key_value_pairs"].keys():
        mymethod["kpoints"] = struct.info["key_value_pairs"]["kpoints"]


    print mymethod

    #which settings for the code?
    codehandler = HandlerFactory(struct, mymethod)

    codehandler.runOne()

    #run the code

    #update the task's status to "running"
    req = requests.patch(SERVER + task['_links']['self'], data={'status': 'running'})
    req.raise_for_status()
    task = req.json()

    #throw the energy into the DB, upload the output file


    #update the task's status to "done"
    req = requests.patch(SERVER + task['_links']['self'], data={'status': 'done'})
    req.raise_for_status()
    task = req.json()

    #upload the data


if __name__=="__main__":
    main()
