import {Component, OnInit, Pipe, PipeTransform, ViewChild, AfterViewChecked} from '@angular/core';
import {OnActivate, Router, RouteSegment} from '@angular/router';

import {DetailsService} from './details.service';

@Component({
  selector: 'decimal-pipe',
  template: `
  <h4>TODO: Proper description of method/settings.</h4>`,
  providers: [DetailsService],
})

export class Details {
  title = 'Periodic Table';
  //  blank_pt: Object[];
  

  errorMessage:  string;
  elements:      Object[];
  methods:       Object[];
  method1:       number;
  method2:       number;

  constructor(private _DetailsService: DetailsService,
              private _router: Router) { };

  myComplete() { 
  };

  getDetails(method1, method2) {
    this._DetailsService.getDetails(method1, method2).subscribe(
      elements => this.elements = elements,
      error => this.errorMessage = <any>error,
         () =>   this.myComplete());
  }

  routerOnActivate(curr: RouteSegment): void {
    this.method1 = +curr.getParam('method1');
    this.method2 = +curr.getParam('method2');
  }

}
