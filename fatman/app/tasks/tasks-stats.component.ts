import {Component, OnInit} from '@angular/core';

import {TaskService} from './tasks.service';
import {TasksList} from './tasks-list.component';

@Component({
  template: `
  <div class="list-group">
    <button
      *ngFor="let tstat of taskstats"
      type="button" class="list-group-item"
      [class.active]="tstat.status === selectedStatus"
      [disabled]="tstat.status === 'total'"
      (click)="onSelect(tstat.status)">
      {{tstat.status}} <span class="badge">{{tstat.count}}</span>
    </button>
    <button type="button" class="list-group-item" disabled>
      Total <span class="badge">{{total}}</span>
    </button>
  </div>
  <tasks-list [statusname]="selectedStatus"></tasks-list>
  `,
  directives: [TasksList]
})

export class TasksStats implements OnInit {
  title = 'Task Stats';

  errorMessage: string;
  taskstats: Object[];
  total: number = 0;
  selectedStatus: string;

  constructor(private _taskService: TaskService) { }

  getStats() {
    this._taskService.getStats().subscribe(
      taskstats => {
        this.total = taskstats.reduce(function(sum: number, entry) {
          return sum + entry['count'];
        }, 0);
        this.taskstats = taskstats;
      },
      error => this.errorMessage = <any>error);

  }

  onSelect(stat: string) { this.selectedStatus = stat; }

  ngOnInit() {
    this.getStats();
  }
}
