#!/usr/bin/env python
"""add_deltatest_to_db.py -YES -- Populate the FATMAN database with all deltatest-related structures.

Use this script on the FATMAN server to populate the database with all chemical structures
that are required for the suite of deltatests.

The structures are taken from the 'dcdft' module in ASE. For each solid-state structure,
5 systems are created with a uniform compression/strain of -2% to +2%.
These 5 systems make up the deltatest for that particular element. The Test is created in
the Test table, and the systems are added to the System table.
Moreover, for each element k-point settings are stored from a hardcoded list based on the
settings used for HGHk/abinit in the Lejaeghere, Science (2016) paper (SI).

Parameters:
    -h     Show this help.
    -YES   Call the program with this switch to add the entries. 
           Otherwise nothing happens to prevent accidental DB clusterfuck.
"""

from __future__ import print_function
from sys import argv

from ase.test.tasks import dcdft

from fatman import db
from fatman.models import Structure, Test, TestStructure
from fatman.tools import atoms2json


def main():
    """Add all the structures and associated tests for the DELTATEST framework to the database"""
    kpts = {
     "H":  (28, 28, 20),  "He": (40, 40, 22),  "Li": (38, 38, 38),  "Be": (52, 52, 28),  "B":  (26, 26, 24),  "C":  (48, 48, 12),  "N":  (16, 16, 16),  "O":  (26, 24, 24),
     "F":  (16, 28, 14),  "Ne": (22, 22, 22),  "Na": (32, 32, 32),  "Mg": (36, 36, 20),  "Al": (24, 24, 24),  "Si": (32, 32, 32),  "P":  (30, 8,  22),  "S":  (38, 38, 38),
     "Cl": (12, 24, 12),  "Ar": (16, 16, 16),  "K":  (20, 20, 20),  "Ca": (18, 18, 18),  "Sc": (34, 34, 20),  "Ti": (40, 40, 22),  "V":  (34, 34, 34),  "Cr": (36, 36, 36),
     "Mn": (28, 28, 28),  "Fe": (36, 36, 36),  "Co": (46, 46, 24),  "Ni": (28, 28, 28),  "Cu": (28, 28, 28),  "Zn": (44, 44, 20),  "Ga": (22, 12, 22),  "Ge": (30, 30, 30),
     "As": (30, 30, 10),  "Se": (26, 26, 20),  "Br": (12, 24, 12),  "Kr": (16, 16, 16),  "Rb": (18, 18, 18),  "Sr": (16, 16, 16),  "Y":  (32, 32, 18),  "Zr": (36, 36, 20),
     "Nb": (30, 30, 30),  "Mo": (32, 32, 32),  "Tc": (42, 42, 22),  "Ru": (42, 42, 24),  "Rh": (26, 26, 26),  "Pd": (26, 26, 26),  "Ag": (24, 24, 24),  "Cd": (38, 38, 18),
     "In": (30, 30, 20),  "Sn": (26, 26, 26),  "Sb": (26, 26,  8),  "Te": (26, 26, 16),  "I":  (12, 22, 10),  "Xe": (14, 14, 14),  "Cs": (16, 16, 16),  "Ba": (20, 20, 20),
     "Hf": (36, 36, 20),  "Ta": (30, 30, 30),  "W":  (32, 32, 32),  "Re": (42, 42, 22),  "Os": (42, 42, 24),  "Ir": (26, 26, 26),  "Pt": (26, 26, 26),  "Au": (24, 24, 24),
     "Hg": (24, 24, 28),  "Tl": (32, 32, 18),  "Pb": (20, 20, 20),  "Bi": (26, 26, 8),   "Po": (30, 30, 30),  "Rn": (14, 14, 14)}

    deltastructures = dcdft.DeltaCodesDFTCollection()
    element_list = deltastructures.keys()


    for el in element_list:
        if el not in kpts.keys():
            continue
        test, created = Test.create_or_get(name="deltatest_{}".format(el), description="Deltatest for element {}".format(el))

        for strain in [0.98, 0.99, 1.00, 1.01, 1.02]:
            data = deltastructures[el]

            eq_cell = data.get_cell()

            struct = data.copy()
            struct.set_cell(eq_cell*strain, scale_atoms=True)

            ase_structure = atoms2json(struct,
                                       additional_information={"dataset": "deltatest",
                                                               "identifier": "deltatest_{}_{:4.2f}".format(el, strain),
                                                               "kpoints": kpts[el]
                                                               }
                                       )

            struct, created = Structure.create_or_get(name="deltatest_{}_{:4.2f}".format(el, strain), ase_structure=ase_structure)
            teststructure, created = TestStructure.create_or_get(structure=struct.id, test=test.id)

if __name__ == "__main__":
    if len(argv) == 2 and argv[1] == "-YES" and '-h' not in argv:
        main()
    else:
        print (__doc__)
