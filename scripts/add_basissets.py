#!/usr/bin/env python
"""add_basissets.py <filename> - parse a cp2k BASIS file and add the basissets to FATMAN

Use this script on the FATMAN server to populate the database with orbital basis sets
from a file with pseudos in cp2k-compatible format.

Parameters:
    - filename: path to the basisset file to be parsed

The basisset family is not automatically created in the db, thus it must exist upfront.
The list of these 'allowed' basisset families is hardcoded ('available_bases').
"""

from sys import argv
from fatman.models import BasissetFamily, BasisSet


def main(args):
    fn = args[0]

    infile = open(fn)
    available_bases = ['SZV-GTH', 'DZV-GTH', 'DZVP-GTH',
                       'TZVP-GTH', 'TZV2P-GTH', 'QZV2P-GTH',
                       'QZV3P-GTH', 'aug-DZVP-GTH', 'aug-TZVP-GTH',
                       'aug-TZV2P-GTH', 'aug-QZV2P-GTH', 'aug-QZV3P-GTH',
                       '6-31G*', '6-311ppG3f2d', '6-31ppG3f2d',
                       'TZVP-pob', 'DZVP-MOLOPT-SR-GTH', 'DZVP-MOLOPT-GTH',
                       'DZVP-ALL', 'DZ-ANO', 'TZVP-MOLOPT-GTH',
                       'TZV2PX-MOLOPT-GTH', 'pc-1', 'pc-2',
                       'pc-3', 'pc-4']
    # only add basissets belonging to these families; must exist in
    # the BasissetFamily table in the db.

    stored_basis = None
    element = None
    bases = None

    for line in infile:
        if line.startswith("#") or line.startswith("!"):
            continue

        # find a line with an element name at the beginning
        if line[0] not in [str(x) for x in range(10)] + [" "]:
            if element is not None and stored_basis is not None:

                # here we save the previously found basis to the DB
                fam = [x for x in bases if x in available_bases]
                if len(fam) == 1:
                    basis = BasissetFamily.get(name=fam[0])

                    b, created = BasisSet.get_or_create(family=basis, element=element, defaults=dict(basis=stored_basis))
                    if created:
                        print "created: ", element, fam[0]

            s = line.split()
            element = s[0]
            bases = s[1:]

            ok_basis = [x in available_bases for x in bases]

            stored_basis = ""

            continue

        stored_basis += line

if __name__ == "__main__":
    if len(argv) == 2:
        main(argv[1:])
    else:
        print __doc__
