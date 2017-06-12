#!/usr/bin/env python3
"""
Check the output of (CP2K) calculation tasks with the error status
for cases where the fdaemon reported an error due to missing accounting data.
"""

import bz2
from urllib.parse import urlsplit

from fatman import app, resultfiles, db
from fatman.models import Calculation, Code, TaskStatus
from fatman.tools import parsers
from fatman.tasks import generate_calculation_results, generate_test_result

with app.app_context():
    # can currently only check CP2K
    for calc in Calculation.query.join(Code).filter(Code.name == "CP2K"):
        task = calc.current_task
        # ignore non-errored tasks
        if task.status.name != "error":
            continue

        resultfile_name = "calc.out"

        try:
            artifact = [a for a in task.outfiles if a.name == resultfile_name][0]
        except IndexError:
            # ok, there is no output file, leave this task as it is
            print("skipping task {}: no {} file found".format(task.id, resultfile_name))
            continue

        try:
            scheme, nwloc, path, _, _ = urlsplit(artifact.path)

            #print("... checking {}".format(path))

            if scheme == 'fkup' and nwloc == 'results':
                filepath = resultfiles.path(path[1:])

                # Artifact.name contains a full path, including possible subdirs
                compressed = artifact.mdata.get('compressed', None)

                if compressed is None:
                    with open(filepath, 'r') as fhandle:
                        results = parsers.get_data_from_output(fhandle, calc.code.name)
                elif compressed == 'bz2':
                    with bz2.open(filepath, mode='rt') as fhandle:
                        results = parsers.get_data_from_output(fhandle, calc.code.name)
                else:
                    # ignore unknown compression formats
                    print("skipping task {}: invalid compression format for {} found".format(task.id, resultfile_name))
                    continue

        except parsers.OutputParseError:
            # ignore tasks where the result is actually not parseable
            print("skipping task {}: output {} not parseable".format(task.id, resultfile_name))
            continue

        print("task {}: result for calculation {}, structure {}, collection {} is actually available".format(
            task.id, calc.id, calc.structure.name, calc.collection.name))

        print(".. setting task to done")
        task.status = TaskStatus.query.filter_by(name="done").one()
        db.session.commit()

        print(".. kick of result/testresult generation")
        generate_calculation_results.apply_async((calc.id, ), link=generate_test_result.si(calc.id))
