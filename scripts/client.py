#!/usr/bin/env python
"""The fatman client queries the DB for new tasks and runs them until none are left"""

from __future__ import print_function

import argparse
import os
import sys
import json
import requests

from fatman.tools import Json2Atoms, randomword, get_data_from_outputfile
from codehandling import HandlerFactory
from time import sleep
from datetime import datetime

SERVER                = 'https://172.23.64.223'
TASKS_URL             = SERVER + '/fatman/tasks'
RESULTS_URL           = SERVER + '/fatman/results'
BASISSET_URL          = SERVER + '/fatman/basis'
PSEUDOPOTENTIAL_URL   = SERVER + '/fatman/pseudo'

DAEMON_SLEEPTIME      = 5*60


def main(args):

    mywd = os.getcwd()
    exitword = os.path.join(mywd, "exit.{:}".format(randomword(6)))


    hpc = args.hpc_mode
    maxtask = args.maxtask
    ntask = 0
    if hpc:
        machinename = 'hpc'
    else:
        machinename = os.uname()[1]

    print("############################################################################")
    print("##                         FATMAN CLIENT                                  ##")
    print("############################################################################")
    print("# {:<30s} {:}".format("Started on: ", datetime.now()))
    print("# {:<30s} {:}".format("This machine: ", machinename))
    print("# {:<30s} {:}".format("Current working directory: ", mywd))
    print("# {:<30s} {:}".format("FATMAN Server: ", SERVER))
    if args.no_calc:
        print("# Calculations will not be run, as per user request")
    if args.no_update:
        print("# Task status will not be updated on server, as per user request")
    if maxtask < 99999:
        print("# Stopping after {:} tasks.".format(maxtask))
    if hpc:
        print("# HPC MODE: Creating and uploading tasks.")
    print("############################################################################")
    print("# To shutdown the client after finishing a task: \n#        touch {:}".format(exitword))
    print("############################################################################")
    sys.stdout.flush()

    sess = requests.Session()

    # the certificate is signed by the inofficial TC-Chem CA
    sess.verify = False

    while ntask < maxtask:
        ntask += 1

        if hpc:
            os.system("rm /tmp/espresso-remote-task/*")

        # check if the 'exit file' exists, then exit
        print("# Checking for {:}...".format(exitword), end='')
        if os.path.exists(exitword):
            os.remove(exitword)
            print("decided to quit.")
            return 0
        print("still running.")

        # request a task
        req = sess.get(TASKS_URL, params={'limit': 1, 'status': 'new'})
        req.raise_for_status()

        tasks = req.json()

        # daemon-like behavior: sleep, and check for new tasks after a few minutes
        if len(tasks) == 0:
            print('{:}: Currently no tasks. Waiting {:} seconds.'.format(datetime.now(), DAEMON_SLEEPTIME))
            sleep(DAEMON_SLEEPTIME)
            continue

        # set the task's status to pending
        if not args.no_update:
            req = sess.patch(SERVER + tasks[0]['_links']['self'], data={'status': 'pending', 'machine': machinename})
        else:
            req = sess.get(SERVER + tasks[0]['_links']['self'])

        req.raise_for_status()
        task = req.json()

        print('{:}: Received task with id = {:}.'.format(datetime.now(), task['id']))

        # Start the actual processing of the task, trying to capture/ignore all possible errors.
        try:
            # get structure and convert to ASE object
            struct_json = task['structure']['ase_structure']
            struct = Json2Atoms(struct_json)

            # which code to use with which settings?
            mymethod = task['method']
            if "kpoints" in struct.info["key_value_pairs"].keys() and "kpoints" not in mymethod['settings'].keys() and "max_kpoints" not in mymethod['settings'].keys():
                # kindof a hack: kpoints are specified with each structure, but are in fact code parameters
                mymethod["kpoints"] = struct.info["key_value_pairs"]["kpoints"]
            elif "kpoints" in mymethod['settings'].keys():
                # kpoints specified in the 'method' override those in the 'structure'
                mymethod["kpoints"] = mymethod['settings']['kpoints']
            elif "max_kpoints" in mymethod['settings'].keys():
                # if max_kpoints is specified, we use whichever kpoint setting is smaller
                if "kpoints" in struct.info["key_value_pairs"].keys():
                    nkp1 = struct.info["key_value_pairs"]["kpoints"][0]*struct.info["key_value_pairs"]["kpoints"][1] * struct.info["key_value_pairs"]["kpoints"][2]
                else:
                    nkp1 = 1e10

                nkp2 = mymethod['settings']['max_kpoints'][0] * mymethod['settings']['max_kpoints'][1] * mymethod['settings']['max_kpoints'][2]
                if nkp2 > nkp1:
                    mymethod["kpoints"] = struct.info["key_value_pairs"]["kpoints"]
                else:
                    mymethod["kpoints"] = mymethod['settings']['max_kpoints']

            if "charge" in struct.info["key_value_pairs"].keys():
                mymethod["charge"] = struct.info["key_value_pairs"]["charge"]

            if "multiplicity" in struct.info["key_value_pairs"].keys():
                mymethod["multiplicity"] = struct.info["key_value_pairs"]["multiplicity"]

            # REMOTE: get dicts containing the Pseudos and Basissets for all required chemical elements
            pseudo = sess.get(PSEUDOPOTENTIAL_URL, data={'family': mymethod['pseudopotential'], 'element': set(struct.get_chemical_symbols())}).json()

            if mymethod["code"] == "cp2k":
                # when running cp2k jobs, the basis and pseudo settings have to be rearranged into a 'kind_settings' structure, while preserving potentially existing kind_settings.
                basis = sess.get(BASISSET_URL, data={'family': mymethod['basis_set'], 'element': set(struct.get_chemical_symbols())}).json()
                if "kind_settings" in mymethod["settings"].keys():
                    ks = mymethod["settings"]["kind_settings"].copy()
                else:
                    ks = {}

                mymethod["kind_settings"] = {}

                for x in set(basis.keys()) & set(pseudo.keys()):
                    mymethod["kind_settings"][x] = {"basis_set": basis[x], "potential": pseudo[x]}
                    for k, v in ks.items():
                        mymethod["kind_settings"][x][k] = v

            elif mymethod["code"] == "espresso":
                # when running espresso jobs the pseudo file has to be written somewhere on disk
                odir = os.path.join("/users/ralph/work/espresso/", mymethod['pseudopotential'])
                try:
                    os.makedirs(odir)
                except OSError:
                    pass

                for e, pp in pseudo.items():
                    of = open(os.path.join(odir, e+".UPF"), 'w')
                    of.write(pp)
                    of.close()

            # REMOTE: update the task's status to "running"
            if not args.no_update and not args.no_calc and not hpc:
                req = sess.patch(SERVER + task['_links']['self'], data={'status': 'running'})
                req.raise_for_status()
                task = req.json()
            elif not args.no_update and not args.no_calc and hpc:
                req = sess.patch(SERVER + task['_links']['self'], data={'status': 'running-remote'})
                req.raise_for_status()
                task = req.json()

            # run the code
            if not args.no_calc and not hpc:
                print('{:}: Calculating task with id = {:}.'.format(datetime.now(), task['id']))
                sys.stdout.flush()

                # create our code-running object with the relevant settings.
                codehandler = HandlerFactory(struct, mymethod)
                e, output_file_path = codehandler.runOne()

                print('{:}: Task {:} done. Energy: {:18.12f}'.format(datetime.now(), task['id'], e))

                extradata = get_data_from_outputfile(output_file_path, mymethod['code'])

                # REMOTE: throw the energy into the DB, get the id of the stored result
                if extradata is not None:
                    req = sess.post(RESULTS_URL, data={'energy': e, 'task_id': task['id'], 'data': json.dumps(extradata, sort_keys=True)})
                else:
                    req = sess.post(RESULTS_URL, data={'energy': e, 'task_id': task['id']})
                req.raise_for_status()
                done_task = req.json()

                # REMOTE: upload the output file
                os.system("bzip2 {}".format(output_file_path))
                with open(output_file_path + ".bz2") as f:
                    req = sess.post(SERVER + done_task["_links"]["self"]+'/file', files={'file': f})
                    req.raise_for_status()

            elif hpc:
                os.system("cp {:}/{:}.UPF /tmp/espresso-remote-task".format(odir, e))

                codehandler = HandlerFactory(struct, mymethod)
                output_file_path = codehandler.createOne()

                runscript_template = (open("./dora_script_espresso.txt")).read()
                params = {"workdir"    : "{:04d}/{:s}".format(mymethod["id"], struct.info["key_value_pairs"]["identifier"]),
                          "taskupdate" : task['_links']['self'],
                          "natom"      : len(struct.get_masses()),
                          "id"         : task['id']}

                of = open(os.path.join("/tmp/espresso-remote-task", "runjob_dora.sh"), "w")
                of.write(runscript_template.format(**params))
                of.close()

                if not args.no_upload:
                    os.system("ssh rkoitz@ela.cscs.ch 'mkdir -p /users/rkoitz/deltatests/espresso/{workdir}/'".format(**params))
                    os.system("scp -rp /tmp/espresso-remote-task/* rkoitz@ela.cscs.ch:/users/rkoitz/deltatests/espresso/{workdir}/".format(**params))
                    if args.do_submit:
                        os.system("ssh rkoitz@ela.cscs.ch \"ssh dora 'cd /users/rkoitz/deltatests/espresso/{workdir}/; /apps/dora/slurm/default/bin/sbatch runjob_dora.sh'\"".format(**params))

            # REMOTE: update the task's status to "done"
            if not args.no_update and not args.no_calc and not hpc:
                req = sess.patch(SERVER + task['_links']['self'], data={'status': 'done'})
                req.raise_for_status()
                task = req.json()

        except Exception as e:
            if not args.no_update:
                req = sess.patch(SERVER + task['_links']['self'], data={'status': 'error'})
                req.raise_for_status()

            print('{:}: Encountered an error! {:}.'.format(datetime.now(), str(e)))

        if not hpc:
            sleep(20)  # just take a short break to avoid cycling through too many tasks if errors happen
        else:
            sleep(1)
        sys.stdout.flush()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-calc',
                        help='Do not run the actual calculation',
                        action='store_true',
                        dest='no_calc')

    parser.add_argument('--no-update',
                        help='Do not update the task\'s status on the server.',
                        action='store_true',
                        dest='no_update')

    parser.add_argument('--no-upload',
                        help='Do not upload the task to the hpc node.',
                        action='store_true',
                        dest='no_upload')

    parser.add_argument('--do-submit',
                        help='Submit the job to SLURM after upload. Without this option the user has to submit the job him/herself.',
                        action='store_true',
                        dest='do_submit')

    parser.add_argument('--hpc-mode',
                        help='Number of HPC calculations to upload (hpc-mode only!)',
                        action='store_true',
                        dest='hpc_mode')

    parser.add_argument('--maxtask',
                        help='Number of tasks after which to stop.',
                        type=int,
                        default=99999,
                        dest='maxtask')

    args = parser.parse_args(sys.argv[1:])

    main(args)
