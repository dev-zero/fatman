#!/usr/bin/env python

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


def test1():
    from ase import Atoms
    teststructure = Atoms('N2',[(0,0,0),(0,0,1.1)])

    teststring = Atoms2Json(teststructure, additional_information={"abc":"def"}) #convert to the string

    newstructure = Json2Atoms(teststring)  #and convert back

    assert (newstructure.get_atomic_numbers() == teststructure.get_atomic_numbers()).all()
    assert (newstructure.get_positions()      == teststructure.get_positions()).all()
    assert (newstructure.get_cell()           == teststructure.get_cell()).all()
    assert (newstructure.get_pbc()            == teststructure.get_pbc()).all()
    assert  newstructure.info["key_value_pairs"]["abc"] == "def"  

    print("apparently no issue")

if __name__=="__main__":
    #if called execute "unit" tests
    test1()
