import {Component, Input, OnChanges, SimpleChange} from '@angular/core';

import {TaskService} from './tasks.service';

@Component({
  selector: 'tasks-list',
  template: `
  <div *ngIf="statusname">
    <table class="tasks table table-bordered table-striped table-condensed">
      <thead>
        <tr>
            <th> Id </th>
            <th> Test </th>
            <th> Machine </th>
            <th> Last update </th>
        </tr>
      </thead>
      <tbody>
          <tr *ngFor="let task of tasks">
              <td> {{task.id}} </td>
              <td> {{task.structure.name}} ({{task.method.pseudopotential}})</td>
              <td> {{task.machine}} </td>
              <td> {{task.mtime}} </td>
          </tr>
      </tbody>
    </table>
  </div>
  `
})

export class TasksList implements OnChanges {
  @Input()
  statusname: string;

  errorMessage: string;
  tasks: Object[];

  constructor(private _taskService: TaskService) { }

  updateTasks() {
    this._taskService.getTasks(this.statusname).subscribe(
      tasks => this.tasks = tasks,
      error => this.errorMessage = <any>error);
  }

  ngOnChanges(changes: {[propertyName: string]: SimpleChange}) {
    for (let propName in changes) {

      if (propName === "statusname")
          this.updateTasks();

    }
  }
}
