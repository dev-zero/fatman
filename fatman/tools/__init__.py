#!/usr/bin/env python

""" Various helper functions used in FATMAN """


import numpy as np
import re
import datetime as dt


def atoms2json(structure, additional_information=None):
    """Serialize an ASE Structure definition to JSON and return it as a string"""

    import json, os
    from ase.db.row import AtomsRow
    from ase.db.core import now
    from ase.io.jsonio import MyEncoder as AseJsonEncoder

    row = AtomsRow(structure)  # this is what ASE would store in its DB
    row.ctime = mtime = now()  # the Row object has an attribute ctime, but not mtime,
                               # we have to wiggle it into the dict later
    row.user = os.getenv("USER")

    dct = row.__dict__.copy()
    del dct["_keys"], dct["_data"], dct["_constraints"]  # containing useless default entries that shouldn't be stored
    dct["mtime"] = mtime
    dct["key_value_pairs"] = additional_information if additional_information else {}

    return json.dumps(dct, sort_keys=True, cls=AseJsonEncoder)


def json2atoms(jsonstring):
    """Read a JSON string and return an Atoms object"""

    from ase.io.jsonio import decode
    from ase.db.row import AtomsRow

    dct = decode(jsonstring)
    row = AtomsRow(dct)

    return row.toatoms(attach_calculator=False, add_additional_information=True)


def nodehours_from_job_data(jobdata):
    """Get the number of node hours as timedelta based on JSON-ified data from sacct.

    alloctres='cpu=576,mem=488000M,node=8'
    """

    try:
        nodes = int(re.search(r"node=(\d+)", jobdata['alloctres']).group(1))
        times = [int(i) for i in jobdata['elapsed'].split(':')]
        return nodes * dt.timedelta(hours=times[0], minutes=times[1], seconds=times[2])
    except (KeyError, AttributeError):
        return None


def mergedicts(dict1, dict2, array_merge_strategy=None):
    """
    Original version from http://stackoverflow.com/a/7205672/1400465
    """
    for k in set(dict1.keys()).union(dict2.keys()):
        if k in dict1 and k in dict2:
            if isinstance(dict1[k], dict) and isinstance(dict2[k], dict):
                yield (k, dict(mergedicts(dict1[k], dict2[k], array_merge_strategy)))
            elif dict2[k] is None:
                # blank-out values in dict1 if their value in dict2 is None
                pass
            elif isinstance(dict1[k], list) and isinstance(dict2[k], list):
                if array_merge_strategy is None:
                    raise RuntimeError("merging of arrays requires an array merge strategy")
                yield (k, array_merge_strategy(k, dict1[k], dict2[k]))
            else:
                # If one of the values is not a dict, you can't continue merging it.
                # Value from second dict overrides one in first and we move on.
                yield (k, dict2[k])
                # Alternatively, replace this with exception raiser
                # to alert you of value conflicts
        elif k in dict1:
            yield (k, dict1[k])
        else:
            yield (k, dict2[k])


def test():
    """Simple test function to test atoms2json and json2atoms function"""

    from ase import Atoms
    teststructure = Atoms('N2', [(0, 0, 0), (0, 0, 1.1)])

    teststring = atoms2json(teststructure, additional_information={"abc":"def"}) #convert to the string

    newstructure = json2atoms(teststring)  #and convert back

    assert (newstructure.get_atomic_numbers() == teststructure.get_atomic_numbers()).all()
    assert (newstructure.get_positions() == teststructure.get_positions()).all()
    assert (newstructure.get_cell() == teststructure.get_cell()).all()
    assert (newstructure.get_pbc() == teststructure.get_pbc()).all()
    assert newstructure.info["key_value_pairs"]["abc"] == "def"

    print("JSON Conversion in both directions: apparently no issue")


if __name__ == "__main__":
    #if called execute "unit" tests
    test()

#  vim: set ts=4 sw=4 tw=0 :
