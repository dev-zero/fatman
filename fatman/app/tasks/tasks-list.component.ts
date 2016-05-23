import {Component, OnInit} from '@angular/core';

import {TaskService} from './tasks.service';

@Component({
  template: `
    <ul class="tasks">
        <li *ngFor="let task of tasks">
            {{task['id']}}: {{task['status']}}
        </li>
    </ul>`,
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
