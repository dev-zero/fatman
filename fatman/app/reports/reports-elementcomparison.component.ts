import {Component, OnInit, OnDestroy, ViewChild, AfterViewChecked} from '@angular/core';
import {Router, ActivatedRoute} from '@angular/router';

import {ReportService} from './reports.service';
import {MethoddetailsComponent} from '../details/methoddetails.component';

@Component({
  directives: [MethoddetailsComponent],
  template: `
      <div style="font-size: 16pt">
        <p>
        Comparing all 
          <select #testselect (change)="onSelect(meth1.value,testselect.value);">
              <option *ngFor="let test of tests" value="{{ test[1] }}">{{ test[1] }}</option>
          </select> to reference method 
          <select #meth1 (change)="onSelect(meth1.value,testselect.value);">
              <option *ngFor="let method of methods" value="{{ method.id }}">{{ method.id }} ({{ method.pseudopotential }})</option>
          </select>.
        </p>
      </div>
    <p>
     <methoddetails method_id="{{ method1 }}"></methoddetails>
    </p>
    <table class="table table-bordered table-striped table-condensed">
      <thead>
        <th (click)="re_sort(0)">Method ID</th>
        <th (click)="re_sort(1)">Method</th>
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
          <td><a href="reports/comparison/{{ method1 }}/{{ line[0] }}"> {{ line[0] }} </a></td>
          <td> {{ line[1] }} </td>
          <td> {{ line[2] | number:'.4-4'  }} </td>
          <td> {{ line[3] | number:'.4-4'  }} </td>
          <td> {{ line[4] | number:'.4-4'  }} </td>
          <td> {{ line[5] | number:'.4-4'  }} </td>
          <td> {{ line[6] | number:'.4-4'  }} </td>
          <td> {{ line[7] | number:'.4-4'  }} </td>
          <td><a href="details/{{ test1 }};methods={{ method1 }},{{ line[0] }}"> {{ line[8] | number:'.4-4' }}</a></td>
        </tr>
      </tbody>
    </table>
  `,
})

export class ReportsElementComparison implements OnInit, OnDestroy {
  @ViewChild('meth1') mymeth1;
  @ViewChild('testselect') mytest1;

  constructor(
    private _route: ActivatedRoute,
    private _router: Router,
    private _service: ReportService) {}

  private _sub: any;

  elements: Object[];
  errorMessage: string;
  method1: number;
  test1: string;
  comparetable: Object[];
  comparelist = [];
  methods: Object[];
  tests: Object[];

  ngOnInit() {
    this._sub = this._route.params.subscribe(params => {
      this.test1 = params['test1'];
      this.method1 = +params['id1'];

      this.getMethods();
      this.getTests();

      this.getElementComparison(this.method1, this.test1);
    });
  }

  ngOnDestroy() {
    this._sub.unsubscribe();
  }

  ngAfterViewChecked() {
    this.CompleteMethod();
    this.CompleteTest();
  }

  getMethods() {
    this._service.getMethods().subscribe(
      methods => this.methods = methods,
      error => this.errorMessage = <any>error,
        () => this.CompleteMethod());
  };

  getTests() {
    this._service.getTests().subscribe(
      tests => this.tests = tests,
      error => this.errorMessage = <any>error,
        () => this.CompleteTest());
  };

  getElementComparison(method1, test1) {
    this._service.getElementComparison(method1, test1).subscribe(
      comparetable => this.comparetable = comparetable,
      error => this.errorMessage = <any>error,
         () => this.myComplete());
  }

  re_sort(col) {
     this.comparelist = this.comparelist.sort(function(a,b) { return a[col] - b[col]});
  };

  myComplete() { 
    var elements =  { "H":1, "He":2, "Li":3, "Be":4, "B":5, "C":6, "N":7, "O":8, "F":9, "Ne":10, "Na":11, "Mg":12, "Al":13, "Si":14, "P":15, "S":16, "Cl":17, "Ar":18, "K":19, "Ca":20, "Sc":21, "Ti":22, "V":23, "Cr":24, "Mn":25, "Fe":26, "Co":27, "Ni":28, "Cu":29, "Zn":30, "Ga":31, "Ge":32, "As":33, "Se":34, "Br":35, "Kr":36, "Rb":37, "Sr":38, "Y":39, "Zr":40, "Nb":41, "Mo":42, "Tc":43, "Ru":44, "Rh":45, "Pd":46, "Ag":47, "Cd":48, "In":49, "Sn":50, "Sb":51,  "Te":52, "I":53, "Xe":54, "Cs":55, "Ba":56, "Hf":72, "Ta":73, "W":74, "Re":75, "Os":76, "Ir":77, "Pt":78, "Au":79, "Hg":80, "Tl":81, "Pb":82, "Bi":83,  "Po":84, "Rn":86 };
     for (var t in this.comparetable['method']) {
        var a = [t].concat(this.comparetable['method'][t]);
        this.comparelist.push(a);
     };
     this.re_sort(0);
   };

  CompleteMethod() { 
     var opts = this.mymeth1.nativeElement.options;
     for(var j = 0; j<opts.length; j++) {
        var opt = opts[j];
        if(opt.value == this.method1) {
            this.mymeth1.nativeElement.selectedIndex = j;
        }
    }
  }
  
  CompleteTest() { 
     var opts = this.mytest1.nativeElement.options;
     for(var j = 0; j<opts.length; j++) {
        var opt = opts[j];
        if(opt.value == this.test1) {
            this.mytest1.nativeElement.selectedIndex = j;
        }
    }
   };

  gotoReports() {
    this._router.navigate(['/reports']);
  }

  onSelect(method1, test1){
    this._router.navigate(['/reports/elementcomparison', method1, test1]);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
