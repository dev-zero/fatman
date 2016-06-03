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
  <h3 *ngIf="elements">Comparing methods {{ f_elements.methods[0] }} and  {{ f_elements.methods[1] }}.</h3>
  <table class="periodictable table table-bordered table-striped table-condensed" *ngIf="f_elements">
      <tr> 
          <td>H <br /><a [routerLink]="['../../details/H', f_elements.methods[0], f_elements.methods[1]]">{{ f_elements.test.deltatest_H  | deltavalue}}</a></td> 
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
        <td>He <br />{{ f_elements.test.deltatest_He | deltavalue}}</td> </tr> 
      <tr> 
        <td>Li <br />{{ f_elements.test.deltatest_Li | deltavalue }}</td> 
        <td>Be <br />{{ f_elements.test.deltatest_Be | deltavalue}}</td> 
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
        <td>B <br />{{ f_elements.test.deltatest_B | deltavalue}}</td> 
        <td>C <br />{{ f_elements.test.deltatest_C | deltavalue}}</td> 
        <td>N <br />{{ f_elements.test.deltatest_N | deltavalue}}</td> 
        <td>O <br />{{ f_elements.test.deltatest_O | deltavalue}}</td> 
        <td>F <br />{{ f_elements.test.deltatest_F | deltavalue}}</td> 
        <td>Ne<br />{{ f_elements.test.deltatest_Ne | deltavalue}}</td> </tr> 
      <tr> 
        <td>Na<br />{{ f_elements.test.deltatest_Na | deltavalue}}</td> 
        <td>Mg<br />{{ f_elements.test.deltatest_Mg | deltavalue}}</td> 
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
        <td>Al<br />{{ f_elements.test.deltatest_Al | deltavalue}}</td> 
        <td>Si<br />{{ f_elements.test.deltatest_Si | deltavalue}}</td> 
        <td>P <br />{{ f_elements.test.deltatest_P  | deltavalue}}</td> 
        <td>S <br />{{ f_elements.test.deltatest_S  | deltavalue}}</td> 
        <td>Cl<br />{{ f_elements.test.deltatest_Cl | deltavalue}}</td> 
        <td>Ar<br />{{ f_elements.test.deltatest_Ar | deltavalue}}</td> </tr> 
      <tr> 
        <td>K <br />{{ f_elements.test.deltatest_K  | deltavalue}}</td> 
        <td>Ca<br />{{ f_elements.test.deltatest_Ca | deltavalue}}</td> 
        <td>Sc<br />{{ f_elements.test.deltatest_Sc | deltavalue}}</td> 
        <td>Ti<br />{{ f_elements.test.deltatest_Ti | deltavalue}}</td> 
        <td>V <br />{{ f_elements.test.deltatest_V  | deltavalue}}</td> 
        <td>Cr<br />{{ f_elements.test.deltatest_Cr | deltavalue}}</td> 
        <td>Mn<br />{{ f_elements.test.deltatest_Mn | deltavalue}}</td> 
        <td>Fe<br />{{ f_elements.test.deltatest_Fe | deltavalue}}</td> 
        <td>Co<br />{{ f_elements.test.deltatest_Co | deltavalue}}</td> 
        <td>Ni<br />{{ f_elements.test.deltatest_Ni | deltavalue}}</td> 
        <td>Cu<br />{{ f_elements.test.deltatest_Cu | deltavalue}}</td> 
        <td>Zn<br />{{ f_elements.test.deltatest_Zn | deltavalue}}</td> 
        <td>Ga<br />{{ f_elements.test.deltatest_Ga | deltavalue}}</td> 
        <td>Ge<br />{{ f_elements.test.deltatest_Ge | deltavalue}}</td> 
        <td>As<br />{{ f_elements.test.deltatest_As | deltavalue}}</td> 
        <td>Se<br />{{ f_elements.test.deltatest_Se | deltavalue}}</td> 
        <td>Br<br />{{ f_elements.test.deltatest_Br | deltavalue}}</td> 
        <td>Kr<br />{{ f_elements.test.deltatest_Kr | deltavalue}}</td> </tr> 
      <tr> 
        <td>Rb <br />{{ f_elements.test.deltatest_Rb| deltavalue}}</td> 
        <td>Sr<br />{{ f_elements.test.deltatest_Sr | deltavalue}}</td> 
        <td>Y <br />{{ f_elements.test.deltatest_Y  | deltavalue}}</td> 
        <td>Zr<br />{{ f_elements.test.deltatest_Zr | deltavalue}}</td> 
        <td>Nb<br />{{ f_elements.test.deltatest_Nb | deltavalue}}</td> 
        <td>Mo<br />{{ f_elements.test.deltatest_Mo | deltavalue}}</td> 
        <td>Tc<br />{{ f_elements.test.deltatest_Tc | deltavalue}}</td> 
        <td>Ru<br />{{ f_elements.test.deltatest_Ru | deltavalue}}</td> 
        <td>Rh<br />{{ f_elements.test.deltatest_Rh | deltavalue}}</td> 
        <td>Pd<br />{{ f_elements.test.deltatest_Pd | deltavalue}}</td> 
        <td>Ag<br />{{ f_elements.test.deltatest_Ag | deltavalue}}</td> 
        <td>Cd<br />{{ f_elements.test.deltatest_Cd | deltavalue}}</td> 
        <td>In<br />{{ f_elements.test.deltatest_In | deltavalue}}</td> 
        <td>Sn<br />{{ f_elements.test.deltatest_Sn | deltavalue}}</td> 
        <td>Sb<br />{{ f_elements.test.deltatest_Sb | deltavalue}}</td> 
        <td>Te<br />{{ f_elements.test.deltatest_Te | deltavalue}}</td> 
        <td>I <br />{{ f_elements.test.deltatest_I  | deltavalue}}</td> 
        <td>Xe<br />{{ f_elements.test.deltatest_Xe | deltavalue}}</td> </tr> 
      <tr> 
        <td>Cs<br />{{ f_elements.test.deltatest_Cs | deltavalue}}</td> 
        <td>Ba<br />{{ f_elements.test.deltatest_Ba | deltavalue}}</td> 
        <td>La<br />{{ f_elements.test.deltatest_La | deltavalue}}</td> 
        <td>Hf<br />{{ f_elements.test.deltatest_Hf | deltavalue}}</td> 
        <td>Ta<br />{{ f_elements.test.deltatest_Ta | deltavalue}}</td> 
        <td>W <br />{{ f_elements.test.deltatest_W  | deltavalue}}</td> 
        <td>Re<br />{{ f_elements.test.deltatest_Re | deltavalue}}</td> 
        <td>Os<br />{{ f_elements.test.deltatest_Os | deltavalue}}</td> 
        <td>Ir<br />{{ f_elements.test.deltatest_Ir | deltavalue}}</td> 
        <td>Pt<br />{{ f_elements.test.deltatest_Pt | deltavalue}}</td> 
        <td>Au<br />{{ f_elements.test.deltatest_Au | deltavalue}}</td> 
        <td>Hg<br />{{ f_elements.test.deltatest_Hg | deltavalue}}</td> 
        <td>Tl<br />{{ f_elements.test.deltatest_Tl | deltavalue}}</td> 
        <td>Pb<br />{{ f_elements.test.deltatest_Pb | deltavalue}}</td> 
        <td>Bi<br />{{ f_elements.test.deltatest_Bi | deltavalue}}</td> 
        <td>Po<br />{{ f_elements.test.deltatest_Po | deltavalue}}</td> 
        <td>At<br />{{ f_elements.test.deltatest_At | deltavalue}}</td> 
        <td>Rn<br />{{ f_elements.test.deltatest_Rn | deltavalue}}</td> </tr> 
  </table>
  <div style="font-size: 16pt">
    Go to Comparison:
      <select #meth1 (change)="onSelect(meth1.value,meth2.value);">
          <option *ngFor="let method of methods">{{ method.id }}</option>
      </select>
      <select #meth2 (change)="onSelect(meth1.value,meth2.value);">
          <option *ngFor="let method of methods">{{ method.id }}</option>
      </select>
  </div>
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
