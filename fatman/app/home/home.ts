import { Component, OnInit, ViewChild } from '@angular/core';
import { HomeService } from './home.service';

@Component({
  providers: [HomeService],
  template: `
    <h3>Welcome to FATMAN</h3>

    <div style="font-size: 12pt" class="container" style="overflow-y:scroll; height:600px;">
      <p>Reference Method:
        <select #meth1 (change)="onSelect(meth1.value);">
            <option *ngFor="let method of methods" value="{{ method.id }}">{{ method.id }} ({{ method.pseudopotential }})</option>
        </select>
      </p>

      <table class="table table-bordered table-striped table-condensed">
        <thead>
          <tr>
            <th (click)="re_sort(0)">Code</th>
            <th (click)="re_sort(2)">Average Delta</th>
          </tr>
        </thead>
        <tbody>
          <tr *ngFor="let id of f_results">
            <td><methoddetails method_id="{{ id[0] }}" small=True></methoddetails></td>
            <td>
              <a [routerLink]="['/periodictable',method1, id[0]]">
                {{ id[1] | number:'.3'}} &plusmn; {{ id[2] | number:'.3'}} ({{ id[3] | number:'.0'}})
              </a>
            </td>
          </tr>
        </tbody>
      </table>
    </div>
  `,
})

export class Home implements OnInit {
  @ViewChild('meth1') mymeth1;

  constructor(
    private _service: HomeService) { }

  errorMessage: string;
  method1 = 3;

  methods: Object[];
  results: Object[];
  f_results = [];


  getMethods() {
    this._service.getMethods().subscribe(
      methods => this.methods = methods,
      error => this.errorMessage = <any>error,
          ()=>this.completeMethod() );
  };

  getResults(method) {
    this.f_results = [];
    if (this.methods) {
      for (var i=0; i<this.methods.length; i++) {
        this._service.getResults(method, this.methods[i]['id']).subscribe(
          results => {
            this.f_results.push([results['methods'][1]]
                                .concat(results['summary']['avg'])
                                .concat(results['summary']['stdev'])
                                .concat(results['summary']['N']))
          },
          error => this.errorMessage = <any>error,
          () => this.re_sort(1)
        );
      }
    }
  };

  ngOnInit(): void {
    this.getMethods();
    this.getResults(this.method1);
  }

  onSelect(method1){
    this.method1 = method1;
    this.getResults(method1);
  }

  completeMethod() { 
     var opts = this.mymeth1.nativeElement.options;
     for(var j = 0; j<opts.length; j++) {
        var opt = opts[j];
        if(opt.value == this.method1) {
            this.mymeth1.nativeElement.selectedIndex = j;
        }
     }
     this.getResults(this.method1);
  };

//completeResult() {
//  for (var i =0;i<this.results['methods'].length; i++) {
//      if (this.results['method'][this.results['methods'][i]][1]>=0){
//         this.f_results.push([this.results['methods'][i]].concat(this.results['method'][this.results['methods'][i]]));
//      }
//      }
//   this.re_sort(2);
//}

  re_sort(col) {
     this.f_results = this.f_results.sort(function(a,b) { return a[col] - b[col]});
  };

//ngAfterViewChecked() {
//  this.completeMethod();
//}
}
//  vim: set ts=2 sw=2 tw=0 :
