#!/usr/bin/env python

import numpy as np
import re
import datetime as dt



def Atoms2Json(structure, additional_information={}):
    """Serialize an ASE Structure definition to JSON and return it as a string"""

    import json,os
    from ase.db.row import AtomsRow
    from ase.db.core import now
    from ase.io.jsonio import MyEncoder as AseJsonEncoder

    row = AtomsRow(structure)       #this is what ASE would store in its DB.
    row.ctime = mtime = now()       #the Row object has an attribute ctime, but not mtime -- we have to wiggle it into the dict later.
    row.user = os.getenv("USER")

    dct = row.__dict__.copy()
    del dct["_keys"], dct["_data"],dct["_constraints"]   #containing useless default entries that shouldn't be stored
    dct["mtime"] = mtime
    dct["key_value_pairs"] = additional_information 

    return json.dumps(dct, sort_keys=True, cls=AseJsonEncoder)


def Json2Atoms(jsonstring):
    """Read a JSON string and return an Atoms object"""

    from ase.io.jsonio import decode
    from ase.db.row import AtomsRow

    dct = decode(jsonstring)
    row = AtomsRow(dct)

    return row.toatoms(attach_calculator=False, add_additional_information=True)


def randomword(length=8):
    import random, string
    return ''.join(random.choice(string.ascii_lowercase + string.digits) for i in range(length))


def energies_for_test(method, test):
    V = []
    E = []
    r = result(method=method,test=test)
    for x in r:
        natom = len(Json2Atoms(x["task"]['structure']['ase_structure']).get_masses())
        V.append(Json2Atoms(x["task"]['structure']['ase_structure']).get_volume()/natom)
        E.append(x['energy']/natom)
        
    return np.array(V),np.array(E)


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


def test():
    from ase import Atoms
    teststructure = Atoms('N2',[(0,0,0),(0,0,1.1)])

    teststring = Atoms2Json(teststructure, additional_information={"abc":"def"}) #convert to the string

    newstructure = Json2Atoms(teststring)  #and convert back

    assert (newstructure.get_atomic_numbers() == teststructure.get_atomic_numbers()).all()
    assert (newstructure.get_positions()      == teststructure.get_positions()).all()
    assert (newstructure.get_cell()           == teststructure.get_cell()).all()
    assert (newstructure.get_pbc()            == teststructure.get_pbc()).all()
    assert  newstructure.info["key_value_pairs"]["abc"] == "def"  

    print("JSON Conversion in both directions: apparently no issue")


if __name__=="__main__":
    #if called execute "unit" tests
    test()

#  vim: set ts=4 sw=4 tw=0 :
