import {Component, OnInit, ViewChild, AfterViewChecked} from '@angular/core';
import {OnActivate, Router, RouteSegment, ROUTER_DIRECTIVES} from '@angular/router';
import {Location} from '@angular/common';
//had 

import {DetailsService} from './details.service';
import {MethoddetailsComponent} from '../details/methoddetails.component';
import {TestdetailsComponent} from '../details/testdetails.component';
import {ComparematrixComponent} from '../details/comparematrix.component';

@Component({
  directives: [ROUTER_DIRECTIVES, MethoddetailsComponent, TestdetailsComponent, ComparematrixComponent],
  inputs: ['test1'],
  providers: [DetailsService],
   template: `
   <p>
   <a [routerLink]="['/periodictable/', selectedmethods[0], selectedmethods[1]]"><button type="button" class="btn btn-sm btn-primary">Back to Periodic Table</button></a>
   <a [routerLink]="['/reports/comparison/', selectedmethods[0], selectedmethods[1]]"><button type="button" class="btn btn-sm btn-info">Back to List View</button></a>
   </p>
   <div class="row row-centered">
     <div class="col-md-9 container">
       <img style="max-width:100%;max-height:100%" src="{{ img_url }}" />
     </div>
     <div class="col-md-3 container-fluid">
      <select #testselect (change)="onSelect(testselect.value);" class="form-control">
        <option *ngFor="let test of tests" value="{{ test[1] }}">{{ test[1] }}</option>
      </select>
      <div class="panel panel-default">
        <div *ngFor="let method of unavailablemethods">
           <input type="checkbox" 
             name="selectedboxes"
             value="{{ method }}"
             [checked]="selectedmethods.indexOf(method)>=0"
             (change)="updateSelectedMethods(method,$event)" /><span style='text-decoration: line-through'> <methoddetails method_id="{{ method }}" small=True></methoddetails></span>
         </div>
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
   <div class="row">
     <div class="container">
       <testdetails *ngFor="let method of selectedmethods" method_id="{{ method }}" test="{{ test1 }}"></testdetails>
     </div>
   </div>
   <div class="container col-md-9">
     <!-- this iframe solution is obviously ugly, but it avoids intermixing 'external' javascript from the Chemdoodle WebGL stuff with angular2 -->
     <iframe seamless frameborder=0 marginheight=0 marginwidth=0 src="../structure?test={{ test1 }}&repeat=2&viewer=True&size=600" width=620 height=630></iframe>
   </div>
   <div class="container col-md-3">
     <comparematrix [methods]="selectedmethods" [test]="test1"></comparematrix>
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
    private _router: Router,
    private location: Location) { this.location = location; }


  errorMessage: string;
  method1: number;
  method2: number;

  test1:   string;
  methods: Object[];
  tests:   Object[];
  selectedmethods = [];
  img_url: string;
  unavailablemethods = [];

  updateSelectedMethods(method, event) {
    if (event.target.checked) {
      this.selectedmethods.push(parseInt(event.target.value));
    }
    else {
      var idx = this.selectedmethods.indexOf(parseInt(event.target.value));
      this.selectedmethods.splice(idx,1);
    }
    this.location.go(
      this._router.serializeUrl(
        this._router.createUrlTree(['details', this.test1, {methods: this.selectedmethods}])
      )
    );
    this.img_url = this.getImgURL();
    this.updateUnavailableMethods();
  }

  routerOnActivate(curr: RouteSegment): void {
    this.selectedmethods= [];
    
    var ms = curr.getParam('methods');

    if (typeof(ms)=='string') {
      var methods_str = ms.split(',');
      for (var i =0;i<methods_str.length;i++) {
        this.selectedmethods.push(+methods_str[i]);
      }
    }
    else {
      for (var i =0;i<ms.length;i++) {
        this.selectedmethods.push(+ms[i]);
      }
    }

    this.test1 = curr.getParam('test1');
    this.getMethods(this.test1);

    this.getTests();
    this.img_url = this.getImgURL();
  }

  getMethods(test) {
    this._service.getMethods(test).subscribe(
      methods => this.methods = methods,
      error => this.errorMessage = <any>error,
      () => this.updateUnavailableMethods());
  };

  updateUnavailableMethods() {
    this.unavailablemethods = [];
    for (var i=0; i<this.selectedmethods.length; i++) {
      var found=false;
      for (var j=0; j<this.methods.length; j++) {
        if (this.methods[j]['id'] == this.selectedmethods[i])  {
          found=true;
        }
      }
      if (!found) this.unavailablemethods.push(this.selectedmethods[i]);
    }
  }

  getTests() {
    this._service.getTests().subscribe(
      tests => this.tests = tests,
      error => this.errorMessage = <any>error,
      () => this.myComplete());
  };

  onSelect(test){
    this._router.navigate(['details', test, {methods: this.selectedmethods}])
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
//           (change)="updateSelectedMethods(method,$event)" /><span style='text-decoration: line-through'> <methoddetails method_id="{{ method.id }}" small=True></methoddetails></span>
