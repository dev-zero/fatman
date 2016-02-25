#!/usr/bin/env python
from fatman.models import PseudopotentialFamily,Pseudopotential

def main(args):
    fn = args[0]

    infile = open(fn)
    available_pseudos= ['GTH-PBE', 'GTH-NLCC-PBE', 'GTH-NLCC2015-PBE', 'ALL']

    stored_pseudo = None
    element = None
    pseudos = None

    for line in infile:
        if line.startswith("#") or line.startswith("!"): continue
        
        #find a line with an element name at the beginning
        if line[0] not in [str(x) for x in range(10)] + [" "]:
            if element is not None and stored_pseudo is not None:

                #here we save the previously found pseudo to the DB
                fam = [x for x in pseudos if x in available_pseudos]
                if len(fam)==1:
                    pseudo = PseudopotentialFamily.get(name=fam[0])

                    Pseudopotential.create_or_get(family=pseudo, element=element, pseudo=stored_pseudo)



            s = line.split()
            element = s[0]
            pseudos = s[1:]

            stored_pseudo=""
                
            continue 

        stored_pseudo += line

if __name__=="__main__":
    from sys import argv
    main(argv[1:])










