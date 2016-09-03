import { Component, OnInit } from '@angular/core';

import { DetailsService } from './details.service';

@Component({
  selector:'comparematrix', 
  inputs: ['methods','test'],
  providers: [DetailsService],
  template: `
    <span>
      <div class='panel panel-info'>
        <div class='panel-heading'>
          {{ test }} <a (click)="ngOnInit()">&#8634;</a>
        </div>
        <table class="table table-bordered table-striped table-condensed">
          <thead>
            <th></th>
            <th *ngFor='let i of methods'><b>{{ i }}</b></th>
          </thead>
          <tr *ngFor='let m of methodmatrix; let i = index'><td><b>{{ methods[i] }}</b></td> <td *ngFor='let j of m'> {{ j | number:".3"}} &nbsp;</td></tr>
        </table>
      </div>
    </span>
  `,
})
export class ComparematrixComponent implements OnInit {
  constructor(private _service: DetailsService) { }

  errorMessage: string;
  method: Object[];

  public methods = [];
  public methodmatrix = [];
  public craplist = [];
  public test = "";

  getComparematrix(methods, test) {
    this.methodmatrix = [];
    for (var i=0; i<this.methods.length; i++) {
      this.methodmatrix.push(Array.apply(null, Array(this.methods.length)).map(Number.prototype.valueOf,0.));
    }
    for (var i=0; i<this.methods.length; i++) {
      for (var j=0; j<this.methods.length; j++){
        this._service.getComparisonOneTest(this.methods[i],this.methods[j],this.test).subscribe(
          ret => this.craplist.push([ret['methods'][0], ret['methods'][1], ret['summary']['avg']]),
          error => this.errorMessage = <any>error,
            () => {this.updateMethodmatrix()});
      }
    }
  //this._service.getMethodDetails(methods).subscribe(
  //  method => this.method = method,
  //  error => this.errorMessage = <any>error);
  };

  updateMethodmatrix() {
    for (var x=0; x<this.craplist.length; x++) {
        var i = this.methods.indexOf(this.craplist[x][0]);
        var j = this.methods.indexOf(this.craplist[x][1]);
        var d = this.craplist[x][2];

        this.methodmatrix[i][j] = d;
      }
  }

  ngOnInit(){
    this.getComparematrix(this.methods, this.test);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
