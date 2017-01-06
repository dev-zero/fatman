
import { Component, Input } from '@angular/core';

import { Task2 }         from './task';
import { Task2Service }  from './task.service';

@Component({
  selector: 'task2-detail',
  template: `
    <div *ngIf="task" class="panel panel-default"
      [ngClass]="{
        'panel-danger': task.status == 'error',
        'panel-warning': task.status == 'running',
        'panel-success': task.status == 'done'
        }">

      <div class="panel-heading">
        <h4 class="panel-title">
          Task {{ task.id }} <span class="pull-right">{{ task.status | uppercase }}</span>
        </h4>
      </div>

      <div class="panel-body">
        <dl class="dl-horizontal">
          <dt>Last Update</dt><dd>{{ task.mtime }}</dd>
          <dt>Created</dt><dd>{{ task.ctime }}</dd>
          <dt>Machine</dt><dd>{{ task.machine }}</dd>
        </dl>
      </div>

      <div class="list-group">
        <div class="list-group-item">
          <div class="list-group-item-text row">
            <div class="col-xs-6">
              <h4>Input Files</h4>
              <ul>
                <li *ngFor="let artifact of task.infiles">
                  <a target="_blank" [href]="artifact._links.download">{{ artifact.name }}</a>
                </li>
              </ul>
            </div>
            <div class="col-xs-6">
              <h4>Output files</h4>
              <ul>
                <li *ngFor="let artifact of task.outfiles">
                  <a target="_blank" [href]="artifact._links.download">{{ artifact.name }}</a>
                </li>
              </ul>
            </div>
        </div>
      </div>

    </div>
  `
})
export class Task2DetailComponent {
  @Input()
  task: Task2;
}
