import {Component, OnInit} from '@angular/core';

import {TaskService} from './tasks.service';

@Component({
  template: `
  <table class="tasks table table-bordered table-striped table-condensed">
    <thead>
      <tr>
          <th> Id </th>
          <th> Status </th>
          <th> Test </th>
          <th> Machine </th>
          <th> Last update </th>
      </tr>
    </thead>
    <tbody>
        <tr *ngFor="let task of tasks">
            <td> {{task.id}} </td>
            <td> {{task.status}} </td>
            <td> {{task.structure.name}} </td>
            <td> {{task.machine}} </td>
            <td> {{task.mtime}} </td>
        </tr>
    </tbody>
  </table>
  `,
  providers: [TaskService]
})

export class Tasks implements OnInit {
  title = 'Tasks';

  errorMessage: string;
  tasks: Object[];

  constructor(private _taskService: TaskService) { }

  getTasks() {
    this._taskService.getTasks().subscribe(
      tasks => this.tasks = tasks,
      error => this.errorMessage = <any>error);
  }

  ngOnInit() {
    this.getTasks();
  }
}
