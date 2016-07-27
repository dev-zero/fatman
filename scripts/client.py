#!/usr/bin/env python

from __future__ import print_function

import os
from os import path
import sys
import json
import tempfile
from time import sleep
from datetime import datetime as dt
from datetime import timedelta
try:
    from urllib.parse import urlparse
except ImportError:
    from urlparse import urlparse
import traceback
import subprocess

import requests
import click

from fatman.tools import Json2Atoms, randomword, get_data_from_outputfile
from codehandling import HandlerFactory

TASKS_URL           = '{}/tasks'
RESULTS_URL         = '{}/results'
BASISSET_URL        = '{}/basis'
PSEUDOPOTENTIAL_URL = '{}/pseudos'

DAEMON_SLEEPTIME = 5*60

# some code in run() changes the working directory, preserve it here
SCRIPTDIR = path.dirname(path.abspath(__file__))

def get_runtime_estimate(sess, url, code, machine, structure, fallback=86340):
    '''Try to calculate a runtime estimate based on previous calculations in seconds

    Search criteria for the estimate are: code-name, machine-name and the structure calculated.

    If no estimate can be obtained (because there is no history), 23:59:00 is returned.
    This fallback can be overwritten.

    Returned is 1.2x the maximum runtime found in seconds.'''

    req = sess.get(RESULTS_URL.format(url),
                   params={'code': code, 'calculated_on': machine, 'structure': structure})
    req.raise_for_status()

    print("estimating runtime for calculating {} with {} on {}".format(structure, code, machine))

    runtime_max = 0
    for result in req.json():
        try:
            runtime = result['data']['runtime']
            if runtime > runtime_max:
                runtime_max = int(runtime)
        except (KeyError, TypeError):
            # missing either data or data.runtime, ignore it
            pass

    if runtime_max == 0:
        return fallback

    return int(1.2*runtime_max) # truncate to int (seconds)

@click.command()
@click.option('--url', type=str, default='http://localhost:5000',
              help='The URL where FATMAN is running (default: http://localhost:5000/)')
@click.option('--workdir', type=click.Path(exists=True, resolve_path=True),
              default='./fatman-client', help="Work directory (default: ./fatman-client)")
@click.option('--hpc-mode/--no-hpc-mode', default=False,
              help='whether to run this client for a HPC system')
@click.option('--calc/--no-calc', default=True, help="do not run the actual calculation")
@click.option('--update/--no-update', default=True,
              help="Do not update the task's status on the server")
@click.option('--upload/--no-upload', default=True,
              help="Do not upload the task to the hpc node (only create the directory structure)")
@click.option('--submit/--no-submit', default=True,
              help="Submit the job to SLURM after upload, otherwise you have to submit it manually")
@click.option('--maxtask', '-m', type=int, default=0,
              help="Number of tasks after which to stop (default: never stop)")
@click.option('--exit-on-error/--no-exit-on-error', default=False,
              help="Exit on error (default: no)")
def run(url, workdir, hpc_mode, calc, update, upload, submit, maxtask, exit_on_error):
    """The fatman client queries the DB for new tasks and runs them until none are left"""

    os.chdir(workdir)

    exitword = path.join(workdir, "exit.{:}".format(randomword(6)))

    if hpc_mode:
        machinename = 'hpc'
    else:
        machinename = os.uname()[1]

    print("############################################################################")
    print("##                         FATMAN CLIENT                                  ##")
    print("############################################################################")
    print("# {:<30s} {:}".format("Started on: ", dt.now()))
    print("# {:<30s} {:}".format("This machine: ", machinename))
    print("# {:<30s} {:}".format("Working directory: ", workdir))
    print("# {:<30s} {:}".format("FATMAN URL: ", url))
    if not calc:
        print("# Calculations will not be run, as per user request")
    if not update:
        print("# Task status will not be updated on server, as per user request")
    if maxtask > 0:
        print("# Stopping after {:} tasks.".format(maxtask))
    if hpc_mode:
        print("# HPC MODE: Creating and uploading tasks.")
    print("############################################################################")
    print("# To shutdown the client after finishing a task: \n#        touch {:}".format(exitword))
    print("############################################################################")
    sys.stdout.flush()

    sess = requests.Session()
    # the certificate is signed by the inofficial TC-Chem CA
    sess.verify = False

    parsed_uri = urlparse(url)
    server = '{uri.scheme}://{uri.netloc}'.format(uri=parsed_uri)

    ntask = 0
    while (maxtask == 0) or (ntask < maxtask):
        ntask += 1

        # check if the 'exit file' exists, then exit
        print("# Checking for {:}...".format(exitword), end=' ')
        if path.exists(exitword):
            os.remove(exitword)
            print("decided to quit.")
            return 0
        print("still running.")

        # request a task
        req = sess.get(TASKS_URL.format(url), params={'limit': 1, 'status': 'new'})
        req.raise_for_status()
        tasks = req.json()

        # daemon-like behavior: sleep, and check for new tasks after a few minutes
        if len(tasks) == 0:
            print('{:}: Currently no tasks. Waiting {:} seconds.'.format(dt.now(), DAEMON_SLEEPTIME))
            sleep(DAEMON_SLEEPTIME)
            continue

        # set the task's status to pending
        if update:
            req = sess.patch(server + tasks[0]['_links']['self'], data={'status': 'pending',
                                                                        'machine': machinename})
        else:
            req = sess.get(server + tasks[0]['_links']['self'])

        req.raise_for_status()
        task = req.json()

        print('{:}: Received task with id = {:}.'.format(dt.now(), task['id']))

        # Start the actual processing of the task, trying to capture/ignore all possible errors.
        try:
            # get structure and convert to ASE object
            struct_json = task['structure']['ase_structure']
            struct = Json2Atoms(struct_json)

            # which code to use with which settings?
            mymethod = task['method']
            if "kpoints" in struct.info["key_value_pairs"].keys() \
                    and "kpoints" not in mymethod['settings'].keys() \
                    and "max_kpoints" not in mymethod['settings'].keys():
                # kindof a hack: kpoints are specified with each structure, but are in fact code parameters
                mymethod["kpoints"] = struct.info["key_value_pairs"]["kpoints"]
            elif "kpoints" in mymethod['settings'].keys():
                # kpoints specified in the 'method' override those in the 'structure'
                mymethod["kpoints"] = mymethod['settings']['kpoints']
            elif "max_kpoints" in mymethod['settings'].keys():
                # if max_kpoints is specified, we use whichever kpoint setting is smaller
                if "kpoints" in struct.info["key_value_pairs"].keys():
                    nkp1 = struct.info["key_value_pairs"]["kpoints"][0] * struct.info["key_value_pairs"]["kpoints"][1] * struct.info["key_value_pairs"]["kpoints"][2]
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

            code_workdir_suffix = path.join(mymethod['code'],
                                            "method_{id:04d}".format(**mymethod),
                                            task['structure']['name'],
                                            "task_{id:04d}".format(**task))
            code_workdir = path.join(workdir, code_workdir_suffix)
            os.makedirs(code_workdir)

            # many wrappers/codes create/expect their stuff relative to their working directory
            os.chdir(code_workdir)

            # Code-specific part
            if mymethod["code"] == "cp2k":
                # REMOTE: get dicts containing the Pseudos and Basissets for all required chemical elements
                req = sess.get(PSEUDOPOTENTIAL_URL.format(url),
                               data={'family': mymethod['pseudopotential'],
                                     'element': set(struct.get_chemical_symbols()),
                                     'format': 'CP2K',
                                    })
                req.raise_for_stats()
                pseudo = req.json()

                # when running cp2k jobs, the basis and pseudo settings have to be rearranged into
                # a 'kind_settings' structure, while preserving potentially existing kind_settings.
                req = sess.get(BASISSET_URL.format(url),
                               data={'family': mymethod['basis_set'],
                                     'element': set(struct.get_chemical_symbols())})
                req.raise_for_status()
                basis = req.json()
                if "kind_settings" in mymethod["settings"].keys():
                    ks = mymethod["settings"]["kind_settings"].copy()
                else:
                    ks = {}

                mymethod["kind_settings"] = {}

                # the combination of (family,element,format) is guaranteed to be unique
                pseudo = {p['element']: p['pseudo'] for p in pseudo}

                for x in set(basis.keys()) & set(pseudo.keys()):
                    mymethod["kind_settings"][x] = {"basis_set": basis[x], "potential": pseudo[x]}
                    for k, v in ks.items():
                        mymethod["kind_settings"][x][k] = v

            elif mymethod["code"] == "espresso":
                # REMOTE: get dicts containing the Pseudos and Basissets for all required chemical elements
                req = sess.get(PSEUDOPOTENTIAL_URL.format(url),
                               data={'family': mymethod['pseudopotential'],
                                     'element': set(struct.get_chemical_symbols()),
                                     'format': 'UPF'})
                req.raise_for_status()
                pseudo = req.json()

                # when running espresso jobs the pseudo file has to be written somewhere on disk
                pseudodir = path.join(code_workdir, mymethod['pseudopotential'])
                os.makedirs(pseudodir)

                for pp in pseudo:
                    with open(path.join(pseudodir, "{element}.UPF".format(**pp)), 'w') as of:
                        of.write(pp['pseudo'])

            # REMOTE: update the task's status to "running"
            if update and calc:
                if hpc_mode:
                    req = sess.patch(server + task['_links']['self'], data={'status': 'running-remote'})
                    req.raise_for_status()
                    task = req.json()
                else:
                    req = sess.patch(server + task['_links']['self'], data={'status': 'running'})
                    req.raise_for_status()
                    task = req.json()

            # run the code
            if calc and not hpc_mode:
                print('{:}: Calculating task with id = {:}.'.format(dt.now(), task['id']))
                sys.stdout.flush()

                # create our code-running object with the relevant settings.
                codehandler = HandlerFactory(struct, mymethod, code_workdir)
                e, output_file_path = codehandler.runOne()

                print('{:}: Task {:} done. Energy: {:18.12f}'.format(dt.now(), task['id'], e))

                extradata = get_data_from_outputfile(output_file_path, mymethod['code'])

                # REMOTE: throw the energy into the DB, get the id of the stored result
                if extradata is not None:
                    req = sess.post(RESULTS_URL.format(url), data={'energy': e, 'task_id': task['id'], 'data': json.dumps(extradata, sort_keys=True)})
                else:
                    req = sess.post(RESULTS_URL.format(url), data={'energy': e, 'task_id': task['id']})
                req.raise_for_status()
                done_task = req.json()

                # REMOTE: upload the output file
                os.system("bzip2 {}".format(output_file_path))
                with open(output_file_path + ".bz2", 'rb') as f:
                    req = sess.post(server + done_task["_links"]["self"]+'/file', files={'file': f})
                    req.raise_for_status()

                # REMOTE: update the task's status to "done"
                if update:
                    req = sess.patch(server + task['_links']['self'], data={'status': 'done'})
                    req.raise_for_status()
                    task = req.json()

            elif hpc_mode:
                codehandler = HandlerFactory(struct, mymethod, workdir=code_workdir)
                input_file_path = codehandler.createOne()
                input_file_name = path.relpath(input_file_path, code_workdir)

                runscript_template_path = path.join(SCRIPTDIR, "dora_script_espresso.txt")

                remote_workdir = "/users/timuel/deltatests"
                remote_scratchdir = "/scratch/daint/timuel/deltatests"
                remote_code_workdir = path.join(remote_workdir, code_workdir_suffix)

                runtime_estimate = get_runtime_estimate(sess, url, "espresso", "hpc", task['structure']['name'])

                with open(runscript_template_path, 'r') as runscript_template, \
                        open(path.join(code_workdir, "runjob_dora.sh"), "w") as runscript:

                    params = {"code_workdir": remote_code_workdir,
                              "workdir"     : remote_workdir,
                              "taskupdate"  : task['_links']['self'],
                              "natom"       : len(struct.get_masses()),
                              "id"          : task['id'],
                              "scratchdir"  : path.join(remote_scratchdir, code_workdir_suffix),
                              "projectname" : code_workdir_suffix,
                              "output_file" : path.join(code_workdir_suffix,
                                                        '{}.out.bz2'.format(input_file_name)),
                              "server"      : server,
                              "results_url" : RESULTS_URL.format(url),
                              "runtime"     : str(timedelta(seconds=runtime_estimate))
                             }

                    runscript.write(runscript_template.read().format(**params))

                if upload:
                    subprocess.check_call(["ssh", "timuel@ela.cscs.ch",
                                           "mkdir -p {}".format(remote_code_workdir)])
                    # using shell expansion in the file argument:
                    subprocess.check_call("scp -rp '{}'/* timuel@ela.cscs.ch:'{}'".format(
                        code_workdir, remote_code_workdir), shell=True)

                    if submit:
                        subprocess.check_call(["ssh", "timuel@ela.cscs.ch",
                                               "ssh dora 'cd {}; sbatch runjob_dora.sh'".format(remote_code_workdir)])

        except Exception as e:
            if update:
                req = sess.patch(server + task['_links']['self'], data={'status': 'error'})
                req.raise_for_status()

            print('{}: Encountered an error!'.format(dt.now()))

            if exit_on_error:
                raise
            else:
                traceback.print_exc(file=sys.stdout)

        if not hpc_mode:
            sleep(20) # take a short break to avoid cycling through too many tasks if errors happen

        sys.stdout.flush()

if __name__ == "__main__":
    run()
