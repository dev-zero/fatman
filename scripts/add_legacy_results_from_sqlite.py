#!/usr/bin/env/python
import sys
import sqlite3
from fatman.models import *
from ase.units import Ry
from datetime import datetime
from ase.test.tasks import dcdft

def main(fn):
    #open file, etc
    conn = sqlite3.connect(fn)
    cur = conn.cursor()

    abinit_results(conn,cur)



def abinit_results(conn,cur):
    edb = dcdft.DeltaCodesDFTCollection()

    basis  = BasissetFamily.get(BasissetFamily.name=='planewave')
    pseudo = PseudopotentialFamily.get(PseudopotentialFamily.name=='GTH-PBE')
        

    sql = "SELECT * FROM abinit WHERE status=\"done\" AND pseudo=\"hgh.k\""

    cur.execute(sql)

    for line in cur:
        M, created = Method.get_or_create(basis_set        = basis, 
                             pseudopotential  = pseudo, 
                             code             = "abinit", 
                             settings         = {   "cutoff_pw"  : line[6], 
                                                    "fband"      : 1.5,
                                                    "ecutsm"     : 0.0,
                                                    "tolsym"     : 1e-12,
                                                    "tsmear"     : 0.01,
                                                    "xc"         : "PBE",
                                                    "occopt"     : 3,
                                                    "pawovlp"    : -1,
                                                    "chksymbreak": 0,
                                                    "prtwf"      : 0,
                                                    "prtden"     : 0,
                                                }
                             )

        element = line[4]
        natom = len(edb[element].get_atomic_numbers())
        struct = Structure.get(Structure.name == "deltatest_{}_{:4.2f}".format(line[4],line[5]))
        test = Test.get(Test.name == "deltatest_{}".format(line[4]))
        status = TaskStatus.get(TaskStatus.name =="done")
        ctime = line[12]#datetime.fromtimestamp(line[12])
        mtime = datetime.now()
        machine = "legacy"
        method = M

        T, Tc = Task.get_or_create(
                           structure = struct,
                           method = method,
                           status = status,
                           test = test,
                           machine = machine,
                           defaults = {'ctime':ctime,
                                       'mtime':mtime})
        
        R, Rc = Result.get_or_create(
                           energy = line[11]*natom,
                           task = T,
                           filename = "not available")
        print "Created Task: ", Tc, "     Created Result: ", Rc

main(sys.argv[1])
