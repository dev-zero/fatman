
export interface Calculation {
  id: string;
  code: string;
  collection: string;
  test: string;
  structure: string;
  _links: { [name: string]: string };
  settings?: any;
  tasks?: any;
}
