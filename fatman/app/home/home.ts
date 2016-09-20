import { Component, OnInit, ViewChild } from '@angular/core';
import { HomeService } from './home.service';

@Component({
  providers: [HomeService],
  template: `
    <h3>Welcome to FATMAN</h3>

    <div style="font-size: 12pt" class="container" style="overflow-y:scroll; height:600px;">
      <div class="form-group">
        <label for="methodselection">Reference method:</label>
        <select [(ngModel)]="reference_method" (ngModelChange)="onSelect()"
            class="form-control" id="methodselection">
          <option *ngFor="let method of methods" [ngValue]="method.id">
            {{ method.code }} / {{ method.pseudopotential }} / {{ method.basis_set }}
          </option>
        </select>
      </div>

      <table class="table table-bordered table-striped table-condensed">
        <thead>
          <tr>
            <th>Code</th>
            <th>Average Delta</th>
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
      methods => {
        this.methods = methods;
        this.getResults();
      },
      error => this.errorMessage = <any>error
    );
  }

  getResults() {
    this.f_results = [];
    for (var i=0; i < this.methods.length; i++) {
      this._service.getResults(this.reference_method, this.methods[i]['id']).subscribe(
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

  ngOnInit(): void {
    this.getMethods();
  }

  onSelect(reference_method){
    this.results = [];
    this.getResults();
  }

  re_sort(col) {
     this.f_results = this.f_results.sort(function(a,b) { return a[col] - b[col]});
  };

}
//  vim: set ts=2 sw=2 tw=0 :
