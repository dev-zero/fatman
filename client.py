#!/usr/bin/env python

import requests
from tools import Json2Atoms

SERVER = 'http://localhost:5000'
TASKS_URL = SERVER + '/tasks'

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

    #which structure?

    #which code to use?

    #which settings for the code?



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
