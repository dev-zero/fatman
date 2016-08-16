import {Component, OnInit, OnDestroy, Pipe, PipeTransform, ViewChild, AfterViewChecked} from '@angular/core';
import {Router, ActivatedRoute, ROUTER_DIRECTIVES} from '@angular/router';

import {PeriodictableService} from './periodictable.service';
import {MethoddetailsComponent} from '../details/methoddetails.component';

@Pipe({name: 'deltavalue'})
export class DeltavaluePipe implements PipeTransform {
    transform(value) {
        if (value == -1.0) {
            return "-";
        }
        else {
            return value.toFixed(3);
        }
    }
}

@Component({
  selector: 'periodictable',
  directives: [ROUTER_DIRECTIVES, MethoddetailsComponent],
  templateUrl: 'periodictable/periodictable-list.component.html',
  styleUrls: ['periodictable/periodictable-list.component.css'],
  providers: [PeriodictableService],
  pipes: [DeltavaluePipe],
})

export class Periodictable implements OnInit, OnDestroy {
  @ViewChild('meth1') mymeth1;
  @ViewChild('meth2') mymeth2;
  title = 'Periodic Table';
  //  blank_pt: Object[];
  
  f_elements=  {'test':{'deltatest_H'  : -1.0, 'deltatest_He' : -1.0, 'deltatest_Li' : -1.0, 'deltatest_Be' : -1.0, 'deltatest_B'  : -1.0, 'deltatest_C'  : -1.0, 'deltatest_N'  : -1.0,
                        'deltatest_O'  : -1.0, 'deltatest_F'  : -1.0, 'deltatest_Ne' : -1.0, 'deltatest_Na' : -1.0, 'deltatest_Mg' : -1.0, 'deltatest_Al' : -1.0, 'deltatest_Si' : -1.0,
                        'deltatest_P'  : -1.0, 'deltatest_S'  : -1.0, 'deltatest_Cl' : -1.0, 'deltatest_Ar' : -1.0, 'deltatest_K'  : -1.0, 'deltatest_Ca' : -1.0, 'deltatest_Sc' : -1.0,
                        'deltatest_Ti' : -1.0, 'deltatest_V'  : -1.0, 'deltatest_Cr' : -1.0, 'deltatest_Mn' : -1.0, 'deltatest_Fe' : -1.0, 'deltatest_Co' : -1.0, 'deltatest_Ni' : -1.0,
                        'deltatest_Cu' : -1.0, 'deltatest_Zn' : -1.0, 'deltatest_Ga' : -1.0, 'deltatest_Ge' : -1.0, 'deltatest_As' : -1.0, 'deltatest_Se' : -1.0, 'deltatest_Br' : -1.0,
                        'deltatest_Kr' : -1.0, 'deltatest_Rb' : -1.0, 'deltatest_Sr' : -1.0, 'deltatest_Y'  : -1.0, 'deltatest_Zr' : -1.0, 'deltatest_Nb' : -1.0, 'deltatest_Mo' : -1.0,
                        'deltatest_Tc' : -1.0, 'deltatest_Ru' : -1.0, 'deltatest_Rh' : -1.0, 'deltatest_Pd' : -1.0, 'deltatest_Ag' : -1.0, 'deltatest_Cd' : -1.0, 'deltatest_In' : -1.0,
                        'deltatest_Sn' : -1.0, 'deltatest_Sb' : -1.0, 'deltatest_Te' : -1.0, 'deltatest_I'  : -1.0, 'deltatest_Xe' : -1.0, 'deltatest_Cs' : -1.0, 'deltatest_Ba' : -1.0,
                        'deltatest_La' : -1.0, 'deltatest_Hf' : -1.0, 'deltatest_Ta' : -1.0, 'deltatest_W'  : -1.0, 'deltatest_Re' : -1.0, 'deltatest_Os' : -1.0, 'deltatest_Ir' : -1.0,
                        'deltatest_Pt' : -1.0, 'deltatest_Au' : -1.0, 'deltatest_Hg' : -1.0, 'deltatest_Tl' : -1.0, 'deltatest_Pb' : -1.0, 'deltatest_Bi' : -1.0, 'deltatest_Po' : -1.0,
                        'deltatest_At' : -1.0, 'deltatest_Rn' : -1.0 },
                'methods': [0, 1],
                'summary': {'N':0, 'avg': 0., 'stdev':0. }};

  private _sub: any;

  errorMessage:  string;
  elements:      Object[];
  methods:       Object[];
  method1:       number;
  method2:       number;

  constructor(
    private _route: ActivatedRoute,
    private _router: Router,
    private _PeriodictableService: PeriodictableService) {}

  ngOnInit() {
    this._sub = this._route.params.subscribe(params => {
      this.method1 = +params['method1'];
      this.method2 = +params['method2'];

      this.getMethodlist();
      this.getPeriodictable(this.method1, this.method2);
    });
  }

  ngOnDestroy() {
    this._sub.unsubscribe();
  }

  myComplete() { 
     this.f_elements['methods'][0] = this.method1;
     this.f_elements['methods'][1] = this.method2;
     for (var t in this.elements['test']) {
        this.f_elements['test'][t] = this.elements['test'][t][6];
     }
     this.f_elements['summary']['N'] = this.elements['summary']['N'];
     this.f_elements['summary']['avg'] = this.elements['summary']['avg'];
     this.f_elements['summary']['stdev'] = this.elements['summary']['stdev'];

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

  getPeriodictable(method1, method2) {
    this._PeriodictableService.getPeriodictable(method1, method2).subscribe(
      elements => this.elements = elements,
      error => this.errorMessage = <any>error,
         () =>   this.myComplete());
  }

  getMethodlist() {
    this._PeriodictableService.getMethodlist().subscribe(
      methods => this.methods = methods,
      error => this.errorMessage = <any>error);
  };

  onSelect(method1, method2){
    this._router.navigate(['/periodictable', method1, method2]);
  }
}
