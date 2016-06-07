import {Component, ViewChild, AfterViewChecked} from '@angular/core';
import {OnActivate, Router, RouteSegment, ROUTER_DIRECTIVES} from '@angular/router';

import {ReportService} from './reports.service';

@Component({
  directives: [ROUTER_DIRECTIVES],
  template: `
    <div style="font-size: 16pt">
      Comparing Methods:
        <select #meth1 (change)="onSelect(meth1.value,meth2.value);">
          <option *ngFor="let method of methods" value="{{ method.id }}">{{ method.id }} ({{ method.pseudopotential }})</option>
        </select> and 
        <select #meth2 (change)="onSelect(meth1.value,meth2.value);">
            <option *ngFor="let method of methods" value="{{ method.id }}">{{ method.id }} ({{ method.pseudopotential }})</option>
        </select>.
    </div>
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
          <td><a href="reports/elementcomparison/{{ method1 }}/deltatest_{{ line[1] }}"> {{ line[1] }} </a></td>
          <td> {{ line[2] | number:'.4-4'  }} </td>
          <td> {{ line[3] | number:'.4-4'  }} </td>
          <td> {{ line[4] | number:'.4-4'  }} </td>
          <td> {{ line[5] | number:'.4-4'  }} </td>
          <td> {{ line[6] | number:'.4-4'  }} </td>
          <td> {{ line[7] | number:'.4-4'  }} </td>
          <td><a href="details/{{ method1 }}/{{ method2 }}/deltatest_{{ line[1] }}"> {{ line[8] | number:'.4-4' }}</a></td>
        </tr>
      </tbody>
    </table>
  <div style="font-size: 16pt">
      Average Delta: {{ summary.avg | number: ".3" }} &plusmn; {{ summary.stdev | number: ".3" }} (N = {{ summary.N }})
  </div>
  `,
})

export class ReportsComparison implements OnActivate {
  @ViewChild('meth1') mymeth1;
  @ViewChild('meth2') mymeth2;

  constructor(
    private _service: ReportService,
    private _router: Router) { }

  elements: Object[];
  errorMessage:  string;
  method1:       number;
  method2:       number;
  comparetable:      Object[];
  comparelist = [];
  methods:       Object[];
  summary = {'N':0, 'avg': 0., 'stdev':0. };

  routerOnActivate(curr: RouteSegment): void {
    this.method1 = +curr.getParam('id1');
    this.method2 = +curr.getParam('id2');

    this.getMethods();
    this.getComparison(this.method1, this.method2);
  }

  getMethods() {
    this._service.getMethods().subscribe(
      methods => this.methods = methods,
      error => this.errorMessage = <any>error,
        () => this.CompleteMethod());
  };

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

     this.summary['N'] = this.comparetable['summary']['N'];
     this.summary['avg'] = this.comparetable['summary']['avg'];
     this.summary['stdev'] = this.comparetable['summary']['stdev'];
  };

  CompleteMethod(){
     var opts = this.mymeth1.nativeElement.options;
     for(var j = 0; j<opts.length; j++) {
        var opt = opts[j];
        if(opt.value == this.method1) {
            this.mymeth1.nativeElement.selectedIndex = j;
        }
        if(opt.value == this.method2) {
            this.mymeth2.nativeElement.selectedIndex = j;
        }
    }
  };

  ngAfterViewChecked() {
    this.CompleteMethod();
  }

  gotoReports() {
    this._router.navigate(['/reports']);
  }

  gotoPT() {
    this._router.navigate(['/periodictable', this.method1, this.method2]);
  }

  onSelect(method1, method2){
    this._router.navigate(['/periodictable', method1, method2]);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
