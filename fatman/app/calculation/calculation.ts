
import { Task2 } from '../task2/task';

export interface TestResult {
  test: string;
  data: any;
}

export interface Structure {
  id: string;
  name: string;
}

export interface Calculation {
  id: string;
  code: string;
  collection: string;
  test: string;
  structure: Structure;
  _links: { [name: string]: string };
  settings?: any;
  results_available: boolean;

  current_task?: Task2;
  tasks?: any;
  testresults?: TestResult[];
}
