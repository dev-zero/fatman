import {Component, OnInit, Pipe, PipeTransform} from '@angular/core';

import {PeriodictableService} from './periodictable.service';

@Pipe({name: 'deltavalue'})
export class DeltavaluePipe implements PipeTransform {
    transform(value) {
        if (value === null) {
            return "";
        }
        else {
            return value.toFixed(3);
        }
    }
}

@Component({
  selector: 'decimal-pipe',
  template: `
  <h3 *ngIf="elements">Comparing methods {{ elements.methods[0] }} and  {{ elements.methods[1] }}.</h3>
  <table class="periodictable table table-bordered table-striped table-condensed" *ngIf="elements">
      <tr> 
        <td>H <br />{{ elements.test?.deltatest_H[6]  | deltavalue}}</td> 
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
        <td>He <br />{{ elements.test.deltatest_He[6] | number:'.3'}}</td> </tr> 
      <tr> 
        <td>Li <br />{{ elements.test.deltatest_Li[6] | number:'.3'}}</td> 
        <td>Be <br />{{ elements.test.deltatest_Be[6] | number:'.3'}}</td> 
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
        <td>B <br />{{elements.text.deltatest_B[6] | number:'.3'}}</td> 
        <td>C <br />{{elements.text.deltatest_C[6] | number:'.3'}}</td> 
        <td>N <br />{{elements.text.deltatest_N[6] | number:'.3'}}</td> 
        <td>O <br />{{elements.text.deltatest_O[6] | number:'.3'}}</td> 
        <td>F <br />{{elements.text.deltatest_F[6] | number:'.3'}}</td> 
        <td>Ne <br />{{elements.text.deltatest_Ne[6] | number:'.3'}}</td> </tr> 
      <tr> 
        <td>Na <br />{{elements.text.deltatest_Na[6] | number:'.3'}}</td> 
        <td>Mg <br />{{elements.text.deltatest_Mg[6] | number:'.3'}}</td> 
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
        <td>Al <br />{{elements.text.deltatest_Al[6] | number:'.3'}}</td> 
        <td>Si <br />{{elements.text.deltatest_Si[6] | number:'.3'}}</td> 
        <td>P <br />{{elements.text.deltatest_P[6] | number:'.3'}}</td> 
        <td>S <br />{{elements.text.deltatest_S[6] | number:'.3'}}</td> 
        <td>Cl<br />{{elements.text.deltatest_Cl[6] | number:'.3'}}</td> 
        <td>Ar<br />{{elements.text.deltatest_Ar[6] | number:'.3'}}</td> </tr> 
      <tr> 
        <td>K <br />{{elements.text.deltatest_K[6] | number:'.3'}}</td> 
        <td>Ca<br />{{elements.text.deltatest_Ca[6] | number:'.3'}}</td> 
        <td>Sc<br />{{elements.text.deltatest_Sc[6] | number:'.3'}}</td> 
        <td>Ti<br />{{elements.text.deltatest_Ti[6] | number:'.3'}}</td> 
        <td>V <br />{{elements.text.deltatest_V[6] | number:'.3'}}</td> 
        <td>Cr<br />{{elements.text.deltatest_Cr[6] | number:'.3'}}</td> 
        <td>Mn<br />{{elements.text.deltatest_Mn[6] | number:'.3'}}</td> 
        <td>Fe<br />{{elements.text.deltatest_Fe[6] | number:'.3'}}</td> 
        <td>Co<br />{{elements.text.deltatest_Co[6] | number:'.3'}}</td> 
        <td>Ni<br />{{elements.text.deltatest_Ni[6] | number:'.3'}}</td> 
        <td>Cu<br />{{elements.text.deltatest_Cu[6] | number:'.3'}}</td> 
        <td>Zn<br />{{elements.text.deltatest_Zn[6] | number:'.3'}}</td> 
        <td>Ga<br />{{elements.text.deltatest_Ga[6] | number:'.3'}}</td> 
        <td>Ge<br />{{elements.text.deltatest_Ge[6] | number:'.3'}}</td> 
        <td>As<br />{{elements.text.deltatest_As[6] | number:'.3'}}</td> 
        <td>Se<br />{{elements.text.deltatest_Se[6] | number:'.3'}}</td> 
        <td>Br<br />{{elements.text.deltatest_Br[6] | number:'.3'}}</td> 
        <td>Kr<br />{{elements.text.deltatest_Kr[6] | number:'.3'}}</td> </tr> 
      <tr> 
        <td>Rb <br />{{elements.text.deltatest_Rb[6] | number:'.3'}}</td> 
        <td>Sr<br />{{elements.text.deltatest_Sr[6] | number:'.3'}}</td> 
        <td>Y <br />{{elements.text.deltatest_Y[6] | number:'.3'}}</td> 
        <td>Zr<br />{{elements.text.deltatest_Zr[6] | number:'.3'}}</td> 
        <td>Nb<br />{{elements.text.deltatest_Nb[6] | number:'.3'}}</td> 
        <td>Mo<br />{{elements.text.deltatest_Mo[6] | number:'.3'}}</td> 
        <td>Tc<br />{{elements.text.deltatest_Tc[6] | number:'.3'}}</td> 
        <td>Ru<br />{{elements.text.deltatest_Ru[6] | number:'.3'}}</td> 
        <td>Rh<br />{{elements.text.deltatest_Rh[6] | number:'.3'}}</td> 
        <td>Pd<br />{{elements.text.deltatest_Pd[6] | number:'.3'}}</td> 
        <td>Ag<br />{{elements.text.deltatest_Ag[6] | number:'.3'}}</td> 
        <td>Cd<br />{{elements.text.deltatest_Cd[6] | number:'.3'}}</td> 
        <td>In<br />{{elements.text.deltatest_In[6] | number:'.3'}}</td> 
        <td>Sn<br />{{elements.text.deltatest_Sn[6] | number:'.3'}}</td> 
        <td>Sb<br />{{elements.text.deltatest_Sb[6] | number:'.3'}}</td> 
        <td>Te<br />{{elements.text.deltatest_Te[6] | number:'.3'}}</td> 
        <td>I <br />{{elements.text.deltatest_I[6] | number:'.3'}}</td> 
        <td>Xe<br />{{elements.text.deltatest_Xe[6] | number:'.3'}}</td> </tr> 
      <tr> 
        <td>Cs <br />{{elements.text.deltatest_Cs[6] | number:'.3'}}</td> 
        <td>Ba<br />{{elements.text.deltatest_Ba[6] | number:'.3'}}</td> 
        <td>La<br />{{elements.text.deltatest_La[6] | number:'.3'}}</td> 
        <td>Hf<br />{{elements.text.deltatest_Hf[6] | number:'.3'}}</td> 
        <td>Ta<br />{{elements.text.deltatest_Ta[6] | number:'.3'}}</td> 
        <td>W <br />{{elements.text.deltatest_W[6] | number:'.3'}}</td> 
        <td>Re<br />{{elements.text.deltatest_Re[6] | number:'.3'}}</td> 
        <td>Os<br />{{elements.text.deltatest_Os[6] | number:'.3'}}</td> 
        <td>Ir<br />{{elements.text.deltatest_Ir[6] | number:'.3'}}</td> 
        <td>Pt<br />{{elements.text.deltatest_Pt[6] | number:'.3'}}</td> 
        <td>Au<br />{{elements.text.deltatest_Au[6] | number:'.3'}}</td> 
        <td>Hg<br />{{elements.text.deltatest_Hg[6] | number:'.3'}}</td> 
        <td>Tl<br />{{elements.text.deltatest_Tl[6] | number:'.3'}}</td> 
        <td>Pb<br />{{elements.text.deltatest_Pb[6] | number:'.3'}}</td> 
        <td>Bi<br />{{elements.text.deltatest_Bi[6] | number:'.3'}}</td> 
        <td>Po<br />{{elements.text.deltatest_Po[6] | number:'.3'}}</td> 
        <td>At<br />{{elements.text.deltatest_At[6] | number:'.3'}}</td> 
        <td>Rn<br />{{elements.text.deltatest_Rn[6] | number:'.3'}}</td> </tr> 
  </table>
  <h4>TODO: Dropdown boxes to choose the two methods to compare; currently hardcoded :(</h4>
  <h4>TODO: Angular can't handle binding of non-existant dictionary items.</h4>
  <h4>TODO: Clicking on an element should take the user to a 'Details' page.</h4>`,
  providers: [PeriodictableService],
  pipes: [DeltavaluePipe],
})

export class Periodictable implements OnInit {
  title = 'Periodic Table';

  errorMessage: string;
  elements: Object[];

  constructor(private _PeriodictableService: PeriodictableService) { }

  getPeriodictable(method1, method2) {
    this._PeriodictableService.getPeriodictable(method1, method2).subscribe(
      elements => this.elements = elements,
      error => this.errorMessage = <any>error);
  }

  ngOnInit() {
    this.getPeriodictable(3,72);
  }
}
