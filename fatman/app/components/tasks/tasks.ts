import {Component, OnInit} from 'angular2/core';
import {TaskService} from '../../services/task';

@Component({
  selector: 'tasks',
  templateUrl: 'components/tasks/tasks.html',
  styleUrls: ['components/tasks/tasks.css'],
})

export class Tasks implements OnInit {
  private tasks: Array<Object> = new Array<Object>();
  private _interval: any;

  constructor(private _taskService: TaskService) { }

  ngOnInit() {
    this._taskService.tasks$.subscribe(latestTasks => {
      this.tasks = latestTasks;
    });

    this._taskService.load();
    this._interval = setInterval(() => { this._taskService.load(); }, 10*1000); // every 10 seconds
  }

  ngOnDestroy() {
    clearInterval(this._interval);
  }

}
