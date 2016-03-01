#!/usr/bin/env python

from fatman.models import *
from datetime import datetime
from fatman.tools import Json2Atoms
from ase.test.tasks import dcdft
from ase.units import kJ, Ry


def main():
    """retrieve all the datapoints for the various tests that have been done and compute the deltatest scores"""
    created_count = 0
    attempt_count = 0

    q = Test.select().join(Task).join(Result).group_by(Test)
    available_tests = list(q)

    q = Method.select().join(Task).join(Result).group_by(Method)
    available_methods = list(q)
    
    for t in available_tests:
        #print "TEST: ", t.name
        for m in available_methods:
            #print "  METHOD: ", m, 

            q = Result.select().join(Task).where((Task.test == t)&(Task.method==m))

            if len(q)>0:   #it can happen that a particular no data is available for a particular method 
                created = store_test_result(list(q))
            else:
                created = False

            if created:
                created_count+=1
            attempt_count +=1
            #print "created:", created 

    print "CREATED {} NEW ENTRIES FROM {} RESULTS.".format(created_count, attempt_count)            


def store_test_result(testset):
    """Take a list of 'result' rows and process them according to the type of test that was performed.
       Store a Row object in the TestResults table"""

    testname = testset[0].task.test.name

    if "deltatest" in testname:

        if len(testset) < 5: return None
        #for a deltatest we have to have (at least) 5 structure/volume pairs, to which we fit a curve

        energies = []
        volumes = []
        natom = 0

        method = testset[0].task.method
        ctime = datetime.now()
        test = testset[0].task.test

        for x in range(len(testset)):
            struct = Json2Atoms(testset[x].task.structure.ase_structure)
            natom = len(struct.get_atomic_numbers())

            energies.append(testset[x].energy/natom)
            volumes.append(struct.get_volume()/natom)
        
        v,e,B0,B1,R = ev_curve(volumes, energies)
        if isinstance(v,complex):
            result_data = {"_status" : "unfittable"}
        else:
            result_data = {"_status" : "fitted", 
                           "V"       : v,
                           "E0"      : e,
                           "B0"      : B0,
                           "B1"      : B1,
                           "R"       : R }
        
        res, created = TestResult.get_or_create(test=test, 
                                 method=method,
                                 defaults={'ctime': ctime,
                                           'result_data': result_data })


    else:
        raise RuntimeError("Can't generate a Result entry for test %s" % testname)


    return created 

def ev_curve(volumes,energies):                                                                                                                                                                   
    structuredata = dcdft.DeltaCodesDFTCollection()                                                                                                                                                     
    eos = dcdft.FullEquationOfState(volumes, energies)                                                                                                                                                  
                                                                                                                                                                                                        
    try:                                                                                                                                                                                                     
        v,e,B0, B1, R = eos.fit()                                                                                                                                                                            
    except ValueError:
        print "failure"
    else:
        return v,e,B0/kJ * 1.0e24, B1, R[0] 




if __name__=="__main__":

    main()
