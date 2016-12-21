
export interface Artifact {
  id: string;
  name: string;
  _links: { [name: string]: string };
}

export interface Task2 {
  id: string;
  ctime: Date;
  mtime: Date;
  machine: string;
  priority: number;
  _links: { [name: string]: string };

  outfiles: Artifact[];
  infiles: Artifact[];
}
