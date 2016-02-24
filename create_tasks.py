#!/usr/bin/env python
from datetime import datetime
from fatman import db
from fatman.models import Task, Test, TestStructure, Method, TaskStatus

def main():
    """Given a list of test names and a method, add new tasks to the task list"""

    available_elements = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"]
    desired_tests = ["deltatest_"+x for x in available_elements]

    desired_method = Method.get(Method.id == 6)
    status_new = TaskStatus.get(TaskStatus.name == "new")

    
    all_teststructures = TestStructure.select()

    for x in all_teststructures:
        if x.test.name in desired_tests[:2]:
            Task.create(structure = x.structure, 
                        method    = desired_method, 
                        status    = status_new,
                        ctime     = datetime.now(),
                        mtime     = datetime.now(),
                        machine   = "-")


if __name__=="__main__":
    main()
