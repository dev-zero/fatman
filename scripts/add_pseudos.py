#!/usr/bin/env python

import click

from fatman.models import PseudopotentialFamily, Pseudopotential

AVAILABLE_PSEUDOS = ['GTH-PBE', 'GTH-NLCC-PBE', 'GTH-NLCC2015-PBE', 'ALL']

@click.command()
@click.argument('infile', type=click.File('r'), metavar='<potentialfile>')
@click.argument('elements', type=str, nargs=-1, metavar='[elements..]')
def main(infile, elements):
    """Parse a cp2k PP file and add the pseudos to FATMAN

    Use this script on the FATMAN server to populate the database with pseudopotentials
    from a file with pseudos in cp2k-compatible format.

    The PP family is not automatically created in the db, thus it must exist upfront.
    The list of these 'allowed' pseudopotential families is hardcoded ('available_pseudos').
    """

    # only add pseudos belonging to these families; must exist in
    # the PseudopotentialFamily table in the db.

    stored_pseudo = ""
    element = None
    pseudos = []

    for line in infile:
        if line.startswith("#") or line.startswith("!"):
            continue

        # find a line with an element name at the beginning
        if line[0] not in [str(x) for x in range(10)] + [" "]:

            # are we done parsing the entry?
            if element is not None and len(stored_pseudo) > 0:

                fam = [x for x in pseudos if x in AVAILABLE_PSEUDOS]

                # here we save the previously found pseudo to the DB
                # if one of its aliases is in the list of pseudos
                # and the list of wanted elements is empty (aka 'all')
                # or the element is in the list of elements
                if len(fam) == 1 and (not elements or element in elements):
                    pseudofamily = PseudopotentialFamily.get(name=fam[0])

                    Pseudopotential.create_or_get(family=pseudofamily,
                                                  element=element,
                                                  pseudo=stored_pseudo,
                                                  format='ABINIT')

            s = line.split()
            element = s[0]
            pseudos = s[1:]

            stored_pseudo = ""

            continue

        stored_pseudo += line

if __name__ == '__main__':
    main()

#  vim: set ts=4 sw=4 tw=0 :
