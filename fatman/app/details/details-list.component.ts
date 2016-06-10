import {Component, OnInit, ViewChild, AfterViewChecked} from '@angular/core';
import {OnActivate, Router, RouteSegment, ROUTER_DIRECTIVES} from '@angular/router';

import {DetailsService} from './details.service';
import {MethoddetailsComponent} from '../details/methoddetails.component';
import {TestdetailsComponent} from '../details/testdetails.component';

@Component({
  directives: [ROUTER_DIRECTIVES, MethoddetailsComponent, TestdetailsComponent],
  providers: [DetailsService],
   template: `
   <div class="row row-centered">
     <div class="col-md-9 container">
       <img style="max-width:100%;max-height:100%" src="{{ img_url }}" />
     </div>
     <div class="col-md-3 container-fluid">
      <select #testselect (change)="onSelect(testselect.value);" class="form-control">
        <option *ngFor="let test of tests" value="{{ test[1] }}">{{ test[1] }}</option>
      </select>
      <div class="panel panel-default">
        <div *ngFor="let method of methods">
           <input type="checkbox" 
             name="selectedboxes"
             value="{{ method.id }}"
             [checked]="selectedmethods.indexOf(method.id)>=0"
             (change)="updateSelectedMethods(method,$event)" /><methoddetails method_id="{{ method.id }}" small=True></methoddetails>
         </div>
       </div>
     </div>
   </div>
<!--<div class="row">
     <div class="container">
       <methoddetails *ngFor="let method of selectedmethods" method_id="{{ method }}"></methoddetails>
     </div>
   </div>-->
   <div class="row">
     <div class="container">
       <testdetails *ngFor="let method of selectedmethods" method_id="{{ method }}" test="{{ test1 }}"></testdetails>
     </div>
   </div>
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
//           (change)="updateSelectedMethods(method,$event)" /> Method {{ method.id }} <font size=1>{{ method.pseudopotential }} </font>
