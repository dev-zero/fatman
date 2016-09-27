
export interface PseudoRef {
  id: string;
  format: string;
  links: { [desc: string]: string };
}

export interface Pseudo {
  id: string;
  element: string;
  family: string;
  format: string;
  converted_from: Pseudo;
  converted_pseudos: Pseudo[];
  links: { [desc: string]: string };
}

export function toPseudoRef(r: any) : PseudoRef {
  return <PseudoRef>({
    id: r.id,
    format: r.format,
    links: r._links,
  });
}

export function toPseudo(r: any) : Pseudo {
  return <Pseudo>({
    id: r.id,
    element: r.element,
    family: r.family,
    format: r.format,
    converted_from: Object.getOwnPropertyNames(r.converted_from).length > 0 ?  toPseudoRef(r.converted_from) : null,
    converted_pseudos: [],
    links: r._links,
  });
}

export interface PseudoFamily {
  name: string;  
}

export function toPseudoFamily(r: any) : PseudoFamily {
  return <PseudoFamily>({
    name: r.name,
  });
}


//  vim: set ts=2 sw=2 tw=0 :
