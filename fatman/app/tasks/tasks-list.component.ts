import {Component, OnInit} from '@angular/core';

import {TaskService} from './tasks.service';

@Component({
  template: `
  <div style="position: relative; left: 30px; /*border: 1px solid #73AD21; padding:5px; display: inline-block;*/">
  <table  class="tasks">
    <thead>
    <tr style="border-bottom: 2px solid; border-top: 2px solid;">
        <th style="padding-right:20px"> <b>Id</b> </th>
        <th style="padding-right:20px"> <b>Status</b> </th>
        <th style="padding-right:20px"> <b>Test</b> </th>
        <th style="padding-right:20px"> <b>Machine</b> </th>
        <th style="padding-right:20px"> <b>Last update</b> </th>
    </tr>
    </thead>
    <tbody>
        <tr *ngFor="let task of tasks" style="border-bottom: 1px solid;">
            <td style="padding-right:20px"> {{task['id']}} </td>
            <td style="padding-right:20px"> {{task['status']}} </td>
            <td style="padding-right:20px"> {{task['structure']['name']}} </td>
            <td style="padding-right:20px"> {{task['machine']}} </td>
            <td style="padding-right:20px"> {{task['mtime']}} </td>
        </tr>
    </tbody>
  </table></div>`,
  providers: [TaskService]
})

export class Tasks implements OnInit {
  title = 'Tasks';
  tasks: Object[];

  constructor(private _taskService: TaskService) { }

  getTasks() {
    this._taskService.getTasks().then(tasks => this.tasks = tasks);
  }

  ngOnInit() {
    this.getTasks();
  }
}
