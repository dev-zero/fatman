import { Component, OnInit, OnDestroy, ViewChild, AfterViewChecked } from '@angular/core';
import { Router, ActivatedRoute } from '@angular/router';
import { Location } from '@angular/common';
import { SafeResourceUrl, DomSanitizer } from '@angular/platform-browser';

import { DetailsService } from './details.service';

@Component({
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
     <iframe seamless frameborder=0 marginheight=0 marginwidth=0 [src]="chemdoodle_url" width=620 height=630></iframe>
   </div>
   <div class="container col-md-3">
     <comparematrix [methods]="selectedmethods" [test]="test1"></comparematrix>
   </div>
    <div class="container">
      <router-outlet></router-outlet>
    </div>
  `,
})

export class Details implements OnInit, OnDestroy {
  @ViewChild('testselect') mytest1;

  constructor(
    private _route: ActivatedRoute,
    private _router: Router,
    private _location: Location,
    private _sanitizer: DomSanitizer,
    private _service: DetailsService) {}


  private _sub: any;

  errorMessage: string;
  method1: number;
  method2: number;

  test1:   string;
  methods: Object[];
  tests:   Object[];
  selectedmethods = [];
  img_url: string;
  chemdoodle_url: SafeResourceUrl;
  unavailablemethods = [];

  ngOnInit() {
    this._sub = this._route.params.subscribe(params => {
      this.selectedmethods = [];

      let ms = params['methods'];

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

      this.test1 = params['test1'];
      this.getMethods(this.test1);

      this.getTests();
      this.img_url = this.getImgURL();
      this.chemdoodle_url = this._sanitizer.bypassSecurityTrustResourceUrl(`../structure?test=${this.test1}&repeat=2&viewer=True&size=600`);
    });
  }

  ngOnDestroy() {
    this._sub.unsubscribe();
  }

  updateSelectedMethods(method, event) {
    if (event.target.checked) {
      this.selectedmethods.push(parseInt(event.target.value));
    }
    else {
      var idx = this.selectedmethods.indexOf(parseInt(event.target.value));
      this.selectedmethods.splice(idx,1);
    }
    this._location.go(
      this._router.serializeUrl(
        this._router.createUrlTree(['details', this.test1, {methods: this.selectedmethods}])
      )
    );
    this.img_url = this.getImgURL();
    this.updateUnavailableMethods();
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
