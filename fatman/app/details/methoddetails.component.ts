import {Component, OnInit, ViewChild, AfterViewChecked, Pipe, PipeTransform} from '@angular/core';
import {OnActivate, Router, RouteSegment, ROUTER_DIRECTIVES} from '@angular/router';

import {DetailsService} from './details.service';

@Pipe({name: 'settingspipe'})
export class SettingsPipe implements PipeTransform {
    transform(value) {
      var ret = "";
      var val_s = "";
      for (var key in value) {
        var val = value[key];
        if (typeof val == 'number') {
          val_s = val.toFixed(0);
        }
        else {
          val_s = JSON.stringify (val);
        }
        ret += key + ": " + val_s + ", " ;  //todo: nicer formatting
      }
      return (ret);
      //if (value == -1.0) {
      //    return "-";
      //}
      //else {
      //    return value.toFixed(3);
      //}
    }
}

@Component({
  selector:'methoddetails', 
  inputs: ['method_id'],
  providers: [DetailsService],
  pipes: [SettingsPipe],
  template: `
     <div *ngIf='method' style="margin-bottom: 20px; margin-top:20px; padding:5px; border: 1px solid; width: 1000px; border-color: #777777">
       <table border=0 width="100%" style='padding-left: 5px'>
         <tr><td rowspan=2 width=60px style='text-align: center'>
             <div style='height:12px; font-size: 10pt'>Method</div><div style='color: #AAAAAA; font-weight:bold; font-size:32pt'>{{ method.id }}</div> </td>
         
           <td style='padding-left: 15px; width: 1px; white-space: nowrap' > Code: <b>{{ method.code }}</b> </td>
           <td style='padding-left: 25px; width: 1px; white-space: nowrap'> Pseudopotential: <b>{{ method.pseudopotential }}</b> </td>
           <td style='padding-left: 25px; width: 1px; white-space: nowrap'> Basis Set: <b>{{ method.basis_set }}</b></td>
           <td></td>
        </tr>
       <tr>
        <td  style='padding-left: 15px' colspan=4>{{ method.settings | settingspipe}}&nbsp;</td>
   </tr>
 </table>
     </div>
  `,
})
export class MethoddetailsComponent implements OnInit {
  constructor(private _service: DetailsService) { }

  errorMessage: string;
  method: Object[];

  public method_id: number;

  getMethodDetails(method_id) {
    this._service.getMethodDetails(method_id).subscribe(
      method => this.method = method,
      error => this.errorMessage = <any>error);
  };

  ngOnInit(){
    this.getMethodDetails(this.method_id);
  }

}
//  vim: set ts=2 sw=2 tw=0 :
