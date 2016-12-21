
import { Task2 } from '../task2/task';

export interface Calculation {
  id: string;
  code: string;
  collection: string;
  test: string;
  structure: string;
  _links: { [name: string]: string };
  settings?: any;
  latest_task: Task2;
  results_available: boolean;
  tasks?: any;
}
