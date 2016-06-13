import {Component, OnInit, ViewChild, AfterViewChecked, Pipe, PipeTransform} from '@angular/core';
import {OnActivate, Router, RouteSegment, ROUTER_DIRECTIVES} from '@angular/router';

import {DetailsService} from './details.service';
import {MethoddetailsComponent} from '../details/methoddetails.component';

@Pipe({
  name: 'truncate'
})
export class TruncatePipe {
  transform(value: string, args: string[]) : string {
    let limit = args.length > 0 ? parseInt(args[0], 10) : 16;

    return value.length > limit ? value.substring(0, limit) : value;
  }
}

@Component({
  selector:'testdetails', 
  inputs: ['method_id', 'test'],
  directives: [ROUTER_DIRECTIVES, MethoddetailsComponent],
  pipes: [TruncatePipe],
  providers: [DetailsService],
  template: `
  <div class='panel panel-info'>
    <div class='panel-heading'>
      <methoddetails method_id="{{ method_id }}"></methoddetails>
    </div>
    <ul class='list-group'>
      <li class='list-group-item' *ngFor='let res of results'>
        <span class='text-muted'><a href='../results/{{ res.id }}/file'> {{ res.task.structure.name }}</a>:</span> <span class='small'><b>{{ res.energy }} eV</b>, update: <b>{{ res.task.mtime |truncate:16 }}</b> on <b>{{ res.task.machine }}</b></span>
      </li>
    </ul>
  </div>
  `,
})
export class TestdetailsComponent implements OnInit {
  constructor(private _service: DetailsService) { }

  errorMessage: string;
  results: Object[];

  public method_id: number;
  public test: string;

  getResults(method_id, test) {
    this._service.getResults(method_id,test).subscribe(
      results => this.results = results,
      error => this.errorMessage = <any>error);
  };

  ngOnInit(){
    this.getResults(this.method_id, this.test);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
