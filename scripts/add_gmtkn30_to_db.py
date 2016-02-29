#!/usr/bin/env python

from ase import io
from string import find
import os

from tools import Atoms2Json, Json2Atoms

from fatman import db
from fatman.models import Structure, Test, TestStructure

def main():
    """Add all the structures and associated tests for the 29 subsets of the GMTKN30 database"""
    os.chdir("gmtkn30")

    #some default values for cell sizes for the various subdatabases
    volumes ={ "ACONF": 20., "ADIM6": 22., "ALK6": 15., "AL2X": 15.,
               "BH76": 15., "BHPERI": 25., "BSR36": 20., "CYCONF": 15.,
               "DARC": 15., "DC9": 20., "G21EA": 15., "G21IP": 15.,
               "G2RC": 15., "HEAVY28": 20., "IDISP": 40., "ISO34": 22.,
               "ISOL22": 40., "MB08-165":22., "NBPRC": 20., "O3ADD6": 15.,
               "PA": 22., "PCONF": 25., "RG6": 15., "RSE43": 15.,
               "S22": 25., "SCONF": 20., "SIE11": 25., "W4-08":20., "WATER27": 30.,}

    for adir in os.listdir("./"):
        if os.path.isdir(adir):
            db_id = adir[0:-10]

            #The entire subdatabase is considered 1 test
            test, created = Test.create_or_get(name="GMTKN_{}".format(db_id), description="GMTKN30 Subset {}".format(db_id))

            for afile in os.listdir(adir):
                if afile.endswith(".xyz"):
                    fn = os.path.join(adir, afile)
                    struct_id = afile[0:-4]

                    charge = multiplicity = None
                    with  open(os.path.join(adir, struct_id) +".config") as infile:
                        for line in infile:
                            foo = line.split()
                            if foo[0] == "CHARGE": charge=int(foo[1])
                            if foo[0] == "MULTIPLICITY": multiplicity =int(foo[1])
                    assert charge != None and multiplicity != None

                    #load the xyz and assign a cubic cell based on the volume above
                    mystruc = io.read(fn)
                    mystruc.set_cell([volumes[db_id]]*3,scale_atoms=False)
                    
                    #serialize to json using a derivative of the native ASE methods.
                    #more elegantly one should probably extend the Atoms object with its own routine to output json.
                    ase_structure = Atoms2Json(mystruc,
                                               additional_information = {"dataset":"gmtkn30", "identifier":"gmtkn30_{}_{}".format(db_id,struct_id),"multiplicity":multiplicity, "charge":charge })

                    #assumes the ase_structure field to be of type STRING/VARCHAR/...
                    struct, created = Structure.create_or_get(name="gmtkn30_{}_{}".format(db_id,struct_id),ase_structure=ase_structure)
                    teststructure, created = TestStructure.create_or_get(structure = struct.id, test = test.id)


if __name__=="__main__":
    main()
