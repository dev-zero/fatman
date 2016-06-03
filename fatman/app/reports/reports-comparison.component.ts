import {Component} from '@angular/core';
import {OnActivate, Router, RouteSegment, ROUTER_DIRECTIVES} from '@angular/router';

import {ReportService} from './reports.service';

@Component({
  directives: [ROUTER_DIRECTIVES],
  template: `
    <h3 *ngIf="method1">Comparing methods {{ method1 }} and  {{ method2 }}.</h3>
    <a [routerLink]="['/periodictable', method1, method2]">Go to periodic table view</a>
    <table class="table table-bordered table-striped table-condensed">
      <thead>
        <th (click)="re_sort(0)">z</th>
        <th (click)="re_sort(1)">Element</th>
        <th (click)="re_sort(2)">V0</th>
        <th (click)="re_sort(3)">B0</th>
        <th (click)="re_sort(4)">B1</th>
        <th (click)="re_sort(5)">V0,r</th>
        <th (click)="re_sort(6)">B0,r</th>
        <th (click)="re_sort(7)">B1,r</th>
        <th (click)="re_sort(8)">Delta</th>
      </thead>
      <tbody>
        <tr *ngFor="let line of comparelist">
          <td> {{ line[0] }} </td>
          <td> {{ line[1] }} </td>
          <td> {{ line[2] | number:'.4-4'  }} </td>
          <td> {{ line[3] | number:'.4-4'  }} </td>
          <td> {{ line[4] | number:'.4-4'  }} </td>
          <td> {{ line[5] | number:'.4-4'  }} </td>
          <td> {{ line[6] | number:'.4-4'  }} </td>
          <td> {{ line[7] | number:'.4-4'  }} </td>
          <td> {{ line[8] | number:'.4-4' }} </td>
        </tr>
      </tbody>
    </table>
  `,
})

export class ReportsComparison implements OnActivate {

  constructor(
    private _service: ReportService,
    private _router: Router) { }

  elements: Object[];
  errorMessage:  string;
  method1:       number;
  method2:       number;
  comparetable:      Object[];
  comparelist = [];

  routerOnActivate(curr: RouteSegment): void {
    this.method1 = +curr.getParam('id1');
    this.method2 = +curr.getParam('id2');

    // TODO
    this.getComparison(this.method1, this.method2);
  }

  getComparison(method1, method2) {
    this._service.getComparison(method1, method2).subscribe(
      comparetable => this.comparetable = comparetable,
      error => this.errorMessage = <any>error,
         () => this.myComplete());
  }

  re_sort(col) {
     this.comparelist = this.comparelist.sort(function(a,b) { return a[col] - b[col]});
  };

  myComplete() { 
    var elements =  { "H":1, "He":2, "Li":3, "Be":4, "B":5, "C":6, "N":7, "O":8, "F":9, "Ne":10, "Na":11, "Mg":12, "Al":13, "Si":14, "P":15, "S":16, "Cl":17, "Ar":18, "K":19, "Ca":20, "Sc":21, "Ti":22, "V":23, "Cr":24, "Mn":25, "Fe":26, "Co":27, "Ni":28, "Cu":29, "Zn":30, "Ga":31, "Ge":32, "As":33, "Se":34, "Br":35, "Kr":36, "Rb":37, "Sr":38, "Y":39, "Zr":40, "Nb":41, "Mo":42, "Tc":43, "Ru":44, "Rh":45, "Pd":46, "Ag":47, "Cd":48, "In":49, "Sn":50, "Sb":51,  "Te":52, "I":53, "Xe":54, "Cs":55, "Ba":56, "Hf":72, "Ta":73, "W":74, "Re":75, "Os":76, "Ir":77, "Pt":78, "Au":79, "Hg":80, "Tl":81, "Pb":82, "Bi":83,  "Po":84, "Rn":86 };
     for (var t in this.comparetable['test']) {
        var el = t.replace('deltatest_','');
        var a = [elements[el], el].concat(this.comparetable['test'][t]);
        this.comparelist.push(a);
     };
     this.re_sort(0);

  };

  gotoReports() {
    this._router.navigate(['/reports']);
  }

  gotoPT() {
    this._router.navigate(['/periodictable', this.method1, this.method2]);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
