import {Component, OnInit, Pipe, PipeTransform, ViewChild, AfterViewChecked} from '@angular/core';
import {OnActivate, Router, RouteSegment, ROUTER_DIRECTIVES} from '@angular/router';

import {PeriodictableService} from './periodictable.service';

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

//<td>H <br /><a [routerLink]="['details/H/{{ f_elements.methods[0] }}/{{ f_elements.methods[1] }}']">{{ f_elements.test.deltatest_H  | deltavalue}}</a></td> 

@Component({
  selector: 'decimal-pipe',
  directives: [ROUTER_DIRECTIVES],
  template: `
  <div style="font-size: 16pt">
    Comparing Methods:
      <select #meth1 (change)="onSelect(meth1.value,meth2.value);">
          <option *ngFor="let method of methods">{{ method.id }}</option>
      </select> and 
      <select #meth2 (change)="onSelect(meth1.value,meth2.value);">
          <option *ngFor="let method of methods">{{ method.id }}</option>
      </select>.
  </div>
  <a [routerLink]="['/reports/comparison', method1, method2]">Go to list view</a>
  <table class="periodictable table table-bordered table-striped table-condensed" *ngIf="f_elements">
      <tr> 
          <td>H <br /><a href="../plot?test=deltatest_H&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_H  | deltavalue}}</a></td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td>He <br /><a href="../plot?test=deltatest_He&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_He | deltavalue}}</a></td></tr> 
      <tr> 
          <td>Li <br /><a href="../plot?test=deltatest_Li&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Li | deltavalue }}</a></td> 
          <td>Be <br /><a href="../plot?test=deltatest_Be&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Be | deltavalue}}</a></td> 
        <td style='border:0'></td>
        <td style='border:0'></td>
        <td style='border:0'></td>
        <td style='border:0'></td>
        <td style='border:0'></td>
        <td style='border:0'></td>
        <td style='border:0'></td>
        <td style='border:0'></td>
        <td style='border:0'></td>
        <td style='border:0'></td>
        <td>B <br /><a href="../plot?test=deltatest_B&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_B | deltavalue}}</a></td> 
        <td>C <br /><a href="../plot?test=deltatest_C&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_C | deltavalue}}</a></td> 
        <td>N <br /><a href="../plot?test=deltatest_N&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_N | deltavalue}}</a></td> 
        <td>O <br /><a href="../plot?test=deltatest_O&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_O | deltavalue}}</a></td> 
        <td>F <br /><a href="../plot?test=deltatest_F&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_F | deltavalue}}</a></td> 
        <td>Ne<br /><a href="../plot?test=deltatest_Ne&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Ne | deltavalue}}</a></td></tr> 
      <tr> 
          <td>Na<br /><a href="../plot?test=deltatest_Na&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Na | deltavalue}}</a></td> 
          <td>Mg<br /><a href="../plot?test=deltatest_Mg&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Mg | deltavalue}}</a></td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td style='border:0'> </td> 
        <td>Al<br /><a href="../plot?test=deltatest_Al&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Al | deltavalue}}</a></td> 
        <td>Si<br /><a href="../plot?test=deltatest_Si&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Si | deltavalue}}</a></td> 
        <td>P <br /><a href="../plot?test=deltatest_P&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_P  | deltavalue}}</a></td> 
        <td>S <br /><a href="../plot?test=deltatest_S&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_S  | deltavalue}}</a></td> 
        <td>Cl<br /><a href="../plot?test=deltatest_Cl&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Cl | deltavalue}}</a></td> 
        <td>Ar<br /><a href="../plot?test=deltatest_Ar&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Ar | deltavalue}}</a></td></tr> 
      <tr> 
          <td>K <br /><a href="../plot?test=deltatest_K&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_K  | deltavalue}}</a></td> 
          <td>Ca<br /><a href="../plot?test=deltatest_Ca&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Ca | deltavalue}}</a></td> 
          <td>Sc<br /><a href="../plot?test=deltatest_Sc&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Sc | deltavalue}}</a></td> 
          <td>Ti<br /><a href="../plot?test=deltatest_Ti&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Ti | deltavalue}}</a></td> 
          <td>V <br /><a href="../plot?test=deltatest_V&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_V  | deltavalue}}</a></td> 
          <td>Cr<br /><a href="../plot?test=deltatest_Cr&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Cr | deltavalue}}</a></td> 
          <td>Mn<br /><a href="../plot?test=deltatest_Mn&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Mn | deltavalue}}</a></td> 
          <td>Fe<br /><a href="../plot?test=deltatest_Fe&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Fe | deltavalue}}</a></td> 
          <td>Co<br /><a href="../plot?test=deltatest_Co&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Co | deltavalue}}</a></td> 
          <td>Ni<br /><a href="../plot?test=deltatest_Ni&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Ni | deltavalue}}</a></td> 
          <td>Cu<br /><a href="../plot?test=deltatest_Cu&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Cu | deltavalue}}</a></td> 
          <td>Zn<br /><a href="../plot?test=deltatest_Zn&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Zn | deltavalue}}</a></td> 
          <td>Ga<br /><a href="../plot?test=deltatest_Ga&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Ga | deltavalue}}</a></td> 
          <td>Ge<br /><a href="../plot?test=deltatest_Ge&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Ge | deltavalue}}</a></td> 
          <td>As<br /><a href="../plot?test=deltatest_As&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_As | deltavalue}}</a></td> 
          <td>Se<br /><a href="../plot?test=deltatest_Se&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Se | deltavalue}}</a></td> 
          <td>Br<br /><a href="../plot?test=deltatest_Br&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Br | deltavalue}}</a></td> 
          <td>Kr<br /><a href="../plot?test=deltatest_Kr&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Kr | deltavalue}}</a></td></tr> 
      <tr> 
          <td>Rb <br /><a href="../plot?test=deltatest_Rb&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Rb| deltavalue}}</a></td> 
          <td>Sr<br /><a href="../plot?test=deltatest_Sr&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Sr | deltavalue}}</a></td> 
          <td>Y <br /><a href="../plot?test=deltatest_Y&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Y  | deltavalue}}</a></td> 
          <td>Zr<br /><a href="../plot?test=deltatest_Zr&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Zr | deltavalue}}</a></td> 
          <td>Nb<br /><a href="../plot?test=deltatest_Nb&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Nb | deltavalue}}</a></td> 
          <td>Mo<br /><a href="../plot?test=deltatest_Mo&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Mo | deltavalue}}</a></td> 
          <td>Tc<br /><a href="../plot?test=deltatest_Tc&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Tc | deltavalue}}</a></td> 
          <td>Ru<br /><a href="../plot?test=deltatest_Ru&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Ru | deltavalue}}</a></td> 
          <td>Rh<br /><a href="../plot?test=deltatest_Rh&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Rh | deltavalue}}</a></td> 
          <td>Pd<br /><a href="../plot?test=deltatest_Pd&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Pd | deltavalue}}</a></td> 
          <td>Ag<br /><a href="../plot?test=deltatest_Ag&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Ag | deltavalue}}</a></td> 
          <td>Cd<br /><a href="../plot?test=deltatest_Cd&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Cd | deltavalue}}</a></td> 
          <td>In<br /><a href="../plot?test=deltatest_In&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_In | deltavalue}}</a></td> 
          <td>Sn<br /><a href="../plot?test=deltatest_Sn&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Sn | deltavalue}}</a></td> 
          <td>Sb<br /><a href="../plot?test=deltatest_Sb&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Sb | deltavalue}}</a></td> 
          <td>Te<br /><a href="../plot?test=deltatest_Te&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Te | deltavalue}}</a></td> 
          <td>I <br /><a href="../plot?test=deltatest_I&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_I  | deltavalue}}</a></td> 
          <td>Xe<br /><a href="../plot?test=deltatest_Xe&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Xe | deltavalue}}</a></td></tr> 
      <tr> 
          <td>Cs<br /><a href="../plot?test=deltatest_Cs&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Cs | deltavalue}}</a></td> 
          <td>Ba<br /><a href="../plot?test=deltatest_Ba&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Ba | deltavalue}}</a></td> 
          <td>La<br /><a href="../plot?test=deltatest_La&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_La | deltavalue}}</a></td> 
          <td>Hf<br /><a href="../plot?test=deltatest_Hf&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Hf | deltavalue}}</a></td> 
          <td>Ta<br /><a href="../plot?test=deltatest_Ta&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Ta | deltavalue}}</a></td> 
          <td>W <br /><a href="../plot?test=deltatest_W&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_W  | deltavalue}}</a></td> 
          <td>Re<br /><a href="../plot?test=deltatest_Re&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Re | deltavalue}}</a></td> 
          <td>Os<br /><a href="../plot?test=deltatest_Os&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Os | deltavalue}}</a></td> 
          <td>Ir<br /><a href="../plot?test=deltatest_Ir&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Ir | deltavalue}}</a></td> 
          <td>Pt<br /><a href="../plot?test=deltatest_Pt&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Pt | deltavalue}}</a></td> 
          <td>Au<br /><a href="../plot?test=deltatest_Au&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Au | deltavalue}}</a></td> 
          <td>Hg<br /><a href="../plot?test=deltatest_Hg&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Hg | deltavalue}}</a></td> 
          <td>Tl<br /><a href="../plot?test=deltatest_Tl&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Tl | deltavalue}}</a></td> 
          <td>Pb<br /><a href="../plot?test=deltatest_Pb&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Pb | deltavalue}}</a></td> 
          <td>Bi<br /><a href="../plot?test=deltatest_Bi&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Bi | deltavalue}}</a></td> 
          <td>Po<br /><a href="../plot?test=deltatest_Po&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Po | deltavalue}}</a></td> 
          <td>At<br /><a href="../plot?test=deltatest_At&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_At | deltavalue}}</a></td> 
          <td>Rn<br /><a href="../plot?test=deltatest_Rn&method={{ f_elements.methods[0] }}&method={{ f_elements.methods[1] }}">{{ f_elements.test.deltatest_Rn | deltavalue}}</a></td></tr> 
  </table>
  <h4>TODO: Clicking on an element should take the user to a 'Details' page.</h4>
  <h4>TODO: Proper description of method/settings.</h4>`,
  providers: [PeriodictableService],
  pipes: [DeltavaluePipe],
})

export class Periodictable {
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
                'methods': [0, 1]};

  errorMessage:  string;
  elements:      Object[];
  methods:       Object[];
  method1:       number;
  method2:       number;

  constructor(private _PeriodictableService: PeriodictableService,
              private _router: Router) { };

  myComplete() { 
     this.f_elements['methods'][0] = this.method1;
     this.f_elements['methods'][1] = this.method2;
     for (var t in this.elements['test']) {
        this.f_elements['test'][t] = this.elements['test'][t][6];
     }

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

  routerOnActivate(curr: RouteSegment): void {
    this.method1 = +curr.getParam('method1');
    this.method2 = +curr.getParam('method2');

    this.getMethodlist();
    this.getPeriodictable(this.method1, this.method2);
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
