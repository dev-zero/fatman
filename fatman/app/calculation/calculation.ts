
export interface Calculation {
  id: string;
  code: string;
  collection: string;
  test: string;
  structure: string;
}

export function toCalculation(r: any) : Calculation {
  return <Calculation>({
    id: r.id,
    code: r.code,
    collection: r.collection,
    test: r.test,
    structure: r.structure
  });
}
