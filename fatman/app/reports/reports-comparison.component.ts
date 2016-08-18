import {Component, OnInit, OnDestroy, ViewChild, AfterViewChecked} from '@angular/core';
import {Router, ActivatedRoute, ROUTER_DIRECTIVES} from '@angular/router';

import {ReportService} from './reports.service';
import {MethoddetailsComponent} from '../details/methoddetails.component';

@Component({
  directives: [ROUTER_DIRECTIVES, MethoddetailsComponent],
  templateUrl: 'reports/reports-comparison.component.html',
  styleUrls: ['reports/reports.component.css'],
})

export class ReportsComparison implements OnInit, OnDestroy {
  @ViewChild('meth1') mymeth1;
  @ViewChild('meth2') mymeth2;

  constructor(
    private _route: ActivatedRoute,
    private _router: Router,
    private _service: ReportService) {}

  private _sub: any;

  elements: Object[];
  errorMessage:  string;
  method1:       number;
  method2:       number;
  comparetable:      Object[];
  comparelist = [];
  methods:       Object[];
  summary = {'N':0, 'avg': 0., 'stdev':0. };

  ngOnInit() {
    this._sub = this._route.params.subscribe(params => {
      this.method1 = +params['id1'];
      this.method2 = +params['id2'];

      this.getMethods();
      this.getComparison(this.method1, this.method2);
    });
  }

  ngOnDestroy() {
    this._sub.unsubscribe();
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
