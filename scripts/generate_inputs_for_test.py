#!/usr/bin/env python
"""This script takes as parameters a Test and a Method and creates input files on disc.
Filenames of these input files are then printed to STDOUT."""

from __future__ import print_function

import argparse
import os
import requests
import sys
import json

from fatman.tools import Json2Atoms, randomword, get_data_from_outputfile
from codehandling import HandlerFactory
from time import sleep
from datetime import datetime

SERVER                = 'http://127.0.0.1:5000'
TASKS_URL             = SERVER + '/tasks'
RESULTS_URL           = SERVER + '/results'
BASISSET_URL          = SERVER + '/basis'
PSEUDOPOTENTIAL_URL   = SERVER + '/pseudo'
TESTS_URL             = SERVER + '/tests'
METHODS_URL            = SERVER + '/methods'

DAEMON_SLEEPTIME      = 5*60


def main(args):

    mywd = os.getcwd()

    method_id = args.method
    testname = args.test
    
    # request a task
    req = requests.get(METHODS_URL + '/{}'.format(method_id), verify=False)
    req.raise_for_status()
    mymethod = req.json()

    # request a task
    req = requests.get(TESTS_URL + '/{}'.format(testname), verify=False)
    req.raise_for_status()
    tests = req.json()

    for t in tests:
        # Start the actual processing of the task, trying to capture/ignore all possible errors.
        try:
            # get structure and convert to ASE object
            struct_json = t['structure']['ase_structure']
            struct = Json2Atoms(struct_json)

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

            if mymethod["code"] == "cp2k":
                pseudo = requests.get(PSEUDOPOTENTIAL_URL,
                                      data={
                                          'family': mymethod['pseudopotential'],
                                          'element': set(struct.get_chemical_symbols()),
                                          'format': 'CP2K',
                                          },
                                      verify=False).json()

                # when running cp2k jobs, the basis and pseudo settings have to be rearranged into a 'kind_settings' structure, while preserving potentially existing kind_settings.
                basis = requests.get(BASISSET_URL, data={'family': mymethod['basis_set'], 'element': set(struct.get_chemical_symbols())}, verify=False).json()
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
                pseudo = requests.get(PSEUDOPOTENTIAL_URL,
                                      data={
                                          'family': mymethod['pseudopotential'],
                                          'element': set(struct.get_chemical_symbols()),
                                          'format': 'UPF',
                                          },
                                      verify=False).json()

                # when running espresso jobs the pseudo file has to be written somewhere on disk
                odir = os.path.join("/users/ralph/work/espresso/", mymethod['pseudopotential'])
                try:
                    os.makedirs(odir)
                except OSError:
                    pass

                for pp in pseudo:
                    of = open(os.path.join(odir, "{element}.UPF".format(**pp)), 'w')
                    of.write(pp['pseudo'])
                    of.close()

            # create our code-running object with the relevant settings.
            codehandler = HandlerFactory(struct, mymethod)
            output_file_path = codehandler.createOne()
            print(output_file_path)

        except Exception as e:
            print('{:}: Encountered an error! {:}.'.format(datetime.now(), str(e)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-m', '--method',
                        help='ID number of the desired method.',
                        type=int,
                        required=True,
                        dest='method')

    parser.add_argument('-t', '--test',
                        help='Name of the desired test',
                        type=str,
                        required=True,
                        dest='test')

    args = parser.parse_args(sys.argv[1:])

    main(args)
