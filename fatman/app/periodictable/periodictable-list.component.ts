import { Component, OnInit, OnDestroy, Pipe, PipeTransform, ViewChild, AfterViewChecked } from '@angular/core';
import { ActivatedRoute, Router, ROUTER_DIRECTIVES } from '@angular/router';

import {PeriodictableService} from './periodictable.service';
import {MethoddetailsComponent} from '../details/methoddetails.component';

@Pipe({name: 'deltavalue'})
export class DeltavaluePipe implements PipeTransform {
    transform(value) {
      if (!!value)
        return value.toFixed(3);

      return "-";
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
  
  elements = {
    'test': {},
    'methods': [0, 1],
    'summary': {'N': 0, 'avg': 0., 'stdev': 0.}
  };

  private _sub: any;

  errorMessage: string;

  methods: Object[];
  method1: number;
  method2: number;

  constructor(
    private _route: ActivatedRoute,
    private _router: Router,
    private _PeriodictableService: PeriodictableService) {}

  ngOnInit() {
    this.getMethodlist();
    this._sub = this._route.params.subscribe(params => {
      this.method1 = +params['method1'];
      this.method2 = +params['method2'];

      this.getPeriodictable(this.method1, this.method2);
    });
  }

  ngOnDestroy() {
    this._sub.unsubscribe();
  }

  getMethodlist() {
    this._PeriodictableService.getMethodlist().subscribe(
      methods => this.methods = methods,
      error => this.errorMessage = <any>error);
  };

  getPeriodictable(method1, method2) {
    this.elements.test = {}

    this._PeriodictableService.getPeriodictable(method1, method2).subscribe(
      elements => {
        this.elements['methods'][0] = this.method1;
        this.elements['methods'][1] = this.method2;
        for (var t in elements['test']) {
           this.elements['test'][t] = elements['test'][t][6];
        }
        this.elements['summary']['N'] = elements['summary']['N'];
        this.elements['summary']['avg'] = elements['summary']['avg'];
        this.elements['summary']['stdev'] = elements['summary']['stdev'];
      },
      error => this.errorMessage = <any>error,
    );
  }

  onSelect(method1: number, method2: number) {
    this._router.navigate(['/periodictable', method1, method2]);
  }
}
