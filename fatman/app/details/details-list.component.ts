import {Component, OnInit, ViewChild, AfterViewChecked} from '@angular/core';
import {OnActivate, Router, RouteSegment, ROUTER_DIRECTIVES} from '@angular/router';

import {DetailsService} from './details.service';
import {MethoddetailsComponent} from '../details/methoddetails.component';

@Component({
  directives: [ROUTER_DIRECTIVES, MethoddetailsComponent],
  providers: [DetailsService],
   template: `
   <table><tr>
   <td width="1000px">
       <img width="1000px" src="{{ img_url }}">
   </td>
   <td width="200px">
     <div class="form-group container" style="border: 2px solid; width: 200px">
       <label>
      <select #testselect (change)="onSelect(testselect.value);">
          <option *ngFor="let test of tests" value="{{ test[1] }}">{{ test[1] }}</option>
      </select>
       </label>
     </div>
     <div class="form-group container" style="width: 200px; height: 600px; overflow-y: scroll; border: 2px solid">
       <div *ngFor="let method of methods">
         <label>
           <input type="checkbox" 
             name="selectedboxes"
             value="{{ method.id }}"
             [checked]="selectedmethods.indexOf(method.id)>=0"
             (change)="updateSelectedMethods(method,$event)" /> Method {{ method.id }} <br /> <font size=1>{{ method.pseudopotential }} </font>
         </label>
       </div>
     </div>
   </td>
 </tr></table>
    <methoddetails *ngFor="let method of selectedmethods" method_id="{{ method }}"></methoddetails>
    <div class="container">
      <router-outlet></router-outlet>
    </div>
  `,
})

export class Details implements OnActivate {
  @ViewChild('testselect') mytest1;
  constructor(
    private _service: DetailsService,
    private _router: Router) { }

  errorMessage: string;
  method1: number;
  method2: number;
  test1:   string;
  methods: Object[];
  tests:   Object[];
  selectedmethods = [];
  img_url: string;

  updateSelectedMethods(method, event) {
    if (event.target.checked) {
      this.selectedmethods.push(parseInt(event.target.value));
    }
    else {
      var idx = this.selectedmethods.indexOf(parseInt(event.target.value));
      this.selectedmethods.splice(idx,1);
    }
    this.img_url = this.getImgURL();
    //nothing for now;
  }

  routerOnActivate(curr: RouteSegment): void {
    this.method1 = +curr.getParam('method1');
    this.method2 = +curr.getParam('method2');
    this.test1 = curr.getParam('test1');

    this.selectedmethods.push(this.method1);
    this.selectedmethods.push(this.method2);
    this.getMethods(this.test1);
    this.getTests();

    this.img_url = this.getImgURL();
  }

  getMethods(test) {
    this._service.getMethods(test).subscribe(
      methods => this.methods = methods,
      error => this.errorMessage = <any>error);
  };

  getTests() {
    this._service.getTests().subscribe(
      tests => this.tests = tests,
      error => this.errorMessage = <any>error,
      () => this.myComplete());
  };

  onSelect(test){
    this.test1 = test; 
    this.img_url = this.getImgURL();
    this.getMethods(this.test1);
  }

  myComplete() {
     var opts = this.mytest1.nativeElement.options;
     for(var j = 0; j<opts.length; j++) {
        var opt = opts[j];
        if(opt.value == this.test1) {
            this.mytest1.nativeElement.selectedIndex = j;
        }
      }
  }

  ngAfterViewChecked(){
    this.myComplete();
  }

  getImgURL(): string {
    var url="../plot?test=" + this.test1;

    for (var i = 0; i<this.selectedmethods.length; i++) {
      url+="&method="+this.selectedmethods[i];
    }
    return url;
  }

}
//  vim: set ts=2 sw=2 tw=0 :
