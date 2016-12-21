
import { Component, OnInit } from '@angular/core';
import { Router, ActivatedRoute, Params } from '@angular/router';

import { Task2 }         from './task';
import { Task2Service }  from './task.service';

@Component({
  template: `
  <div *ngIf="task">
    hello, world!
  </div>
  `,
  providers: [
    Task2Service
  ]
})
export class Task2DetailComponent implements OnInit {
  task: Task2;

  constructor(
    private _route: ActivatedRoute,
    private _router: Router,
    private _service: Task2Service
  ) {}

  ngOnInit() {
    this._route.params
      .switchMap((params: Params) => this._service.getTask(params['id']))
      .subscribe((task: Task2) => this.task = task);
  }
}
