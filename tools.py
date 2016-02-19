#!/usr/bin/python

def Atoms2Json(structure, additional_information={}):
    """Serialize an ASE Structure definition to JSON and return it as a string"""

    import json,os
    from ase.db.row import AtomsRow
    from ase.db.core import now
    from ase.io.jsonio import MyEncoder as AseJsonEncoder

    row=AtomsRow(structure)
    row.ctime = mtime = now()
    row.user = os.getenv("USER")

    dct = row.__dict__.copy()
    del dct["_keys"], dct["_data"],dct["_constraints"]
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
