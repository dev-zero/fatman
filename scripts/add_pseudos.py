#!/usr/bin/env python
"""add_pseudos.py <filename> - parse a cp2k PP file and add the pseudos to FATMAN

Use this script on the FATMAN server to populate the database with pseudopotentials
from a file with pseudos in cp2k-compatible format.

Parameters:
    - filename: path to the pseudopotential file to be parsed

The PP family is not automatically created in the db, thus it must exist upfront.
The list of these 'allowed' pseudopotential families is hardcoded ('available_pseudos').
"""

from sys import argv
from fatman.models import PseudopotentialFamily, Pseudopotential


def main(args):
    fn = args[0]

    available_pseudos = ['GTH-PBE', 'GTH-NLCC-PBE', 'GTH-NLCC2015-PBE', 'ALL']
    # only add pseudos belonging to these families; must exist in
    # the PseudopotentialFamily table in the db.

    stored_pseudo = None
    element = None
    pseudos = None

    with open(fn) as infile:
        for line in infile:
            if line.startswith("#") or line.startswith("!"):
                continue

            # find a line with an element name at the beginning
            if line[0] not in [str(x) for x in range(10)] + [" "]:
                if element is not None and stored_pseudo is not None:

                    # here we save the previously found pseudo to the DB
                    fam = [x for x in pseudos if x in available_pseudos]
                    if len(fam) == 1:
                        pseudo = PseudopotentialFamily.get(name=fam[0])

                        Pseudopotential.create_or_get(family=pseudo, element=element, pseudo=stored_pseudo)

                s = line.split()
                element = s[0]
                pseudos = s[1:]

                stored_pseudo = ""

                continue

            stored_pseudo += line

if __name__ == "__main__":
    if len(argv) == 2:
        main(argv[1:])
    else:
        print __doc__
