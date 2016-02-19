#!/usr/bin/env python

from ase.test.tasks import dcdft
from ase.db import connect

from fatman import db
from fatman.models import Structure, Test, TestStructure
from tools import Atoms2Json

def main():
    """Add all the structures and associated tests for the DELTATEST framework to the database"""

    deltastructures = dcdft.DeltaCodesDFTCollection()
    element_list = deltastructures.keys()
    
    con = connect("postgresql:///fatman")

    for el in element_list:
        test, created = Test.create_or_get(name="deltatest_{}".format(el), description="Deltatest for element {}".format(el))

        for strain in [0.98, 0.99, 1.00, 1.01, 1.02]:
            data = deltastructures[el]

            eq_cell = data.get_cell()
            
            struct = data.copy()
            struct.set_cell(eq_cell*strain, scale_atoms=True)

            ase_structure = Atoms2Json(struct, 
                                       additional_information={"dataset":"deltatest","identifier":"deltatest_{}_{:4.2f}".format(el,strain)})

            struct, created = Structure.create_or_get(name="deltatest_{}_{:4.2f}".format(el,strain),ase_structure=ase_structure)
            teststructure, created = TestStructure.create_or_get(structure = struct.id, test = test.id)


if __name__=="__main__":
    main()
