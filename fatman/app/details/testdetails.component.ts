import {Component, OnInit, Pipe, PipeTransform} from '@angular/core';
import {Router} from '@angular/router';

import {DetailsService} from './details.service';
import {MethoddetailsComponent} from '../details/methoddetails.component';

@Pipe({name: 'truncate'})
export class TruncatePipe {
  transform(value: string, args: string[]) : string {
    let limit = args.length > 0 ? parseInt(args[0], 10) : 16;

    return value.length > limit ? value.substring(0, limit) : value;
  }
}

@Component({
  selector:'testdetails', 
  inputs: ['method_id', 'test'],
  directives: [MethoddetailsComponent],
  pipes: [TruncatePipe],
  providers: [DetailsService],
  template: `
  <div class='panel panel-info'>
    <div class='panel-heading'>
      <div class='row'>
        <div class='col-md-10'>
          <methoddetails method_id="{{ method_id }}"></methoddetails>
        </div>
        <div *ngIf='testresults' class='col-md-2'>
          <table>
            <tr><td style='padding-right: 10px'>V<sub>0</sub></td><td> <b>{{ testresults[0][test].V  | number: ".3"}}</b></td></tr>
            <tr><td style='padding-right: 10px'>B<sub>0</sub></td><td> <b>{{ testresults[0][test].B0 | number: ".3"}}</b></td></tr>
            <tr><td style='padding-right: 10px'>B<sub>1</sub></td><td> <b>{{ testresults[0][test].B1 | number: ".3"}}</b></td></tr>
          </table>
        </div>
      </div>
    </div>
    <ul class='list-group'>
      <li class='list-group-item' *ngFor='let res of results'>
        <span class='text-muted'>
          <a [href]="'../results/' + res.id + '/file'" target="_blank"> {{ res.task.structure.name }}</a>:
          </span> 
          <span class='small'>
            <b>{{ res.energy }} eV</b>, update: 
            <b>{{ res.task.mtime |truncate:16 }}</b> on 
            <b>{{ res.task.machine }}</b> 
            <span *ngIf='res.data?.mpiranks>0'>with <b>{{ res.data?.mpiranks }}</b> ranks</span>
            <span *ngIf='res.data?.runtime>0'>, runtime <b>{{ res.data?.runtime | number:'.0-0'}} s</b></span>
            .
            <span *ngIf='res.data?.user'> User: <b>{{ res.data?.username }}</b></span>
            <span *ngIf='res.data?.nkpoints>0'> k-points: <b>{{ res.data?.nkpoints }}</b></span>
        </span>
      </li>
    </ul>
  </div>
  `,
})
export class TestdetailsComponent implements OnInit {
  constructor(private _service: DetailsService) { }

  errorMessage: string;
  results: Object[];
  testresults: Object[];

  public method_id: number;
  public test: string;

  getResults(method_id, test) {
    this._service.getResults(method_id,test).subscribe(
      results => this.results = results,
      error => this.errorMessage = <any>error);
  };

  getTestresults(method_id, test) {
    this._service.getTestresults(method_id,test).subscribe(
      testresults => this.testresults = testresults,
      error => this.errorMessage = <any>error);
  };

  ngOnInit(){
    this.getResults(this.method_id, this.test);
    this.getTestresults(this.method_id, this.test);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
