import {Component, OnInit} from '@angular/core';

import {TaskService} from '../../services/task';

@Component({
  selector: 'tasks',
  templateUrl: 'components/tasks/tasks.html',
  styleUrls: ['components/tasks/tasks.css'],
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
