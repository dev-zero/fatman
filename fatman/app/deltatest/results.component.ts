import { Component, OnInit } from '@angular/core';

import { DeltatestResultsService } from './results.service';
import { TestResult } from './testresult';

@Component({
  template: `
  <div class="row">
    <div class="col-md-12">
      <div class="form-group">
        <label for="methodselection">Select reference method:</label>
        <select [(ngModel)]="reference_method" (ngModelChange)="onSelect()"
            class="form-control" id="methodselection">
          <option *ngFor="let method of methods" [ngValue]="method.id">
            {{ method.id }} ({{ method.pseudopotential }})
          </option>
        </select>
      </div>

  <div class="row">
	<div class="col-md-12">
	  <methoddetails method_id="{{ reference_method }}"></methoddetails>
	</div>
  </div>

  <div class="row">
    <div class="col-md-12">
      <table class="table table-bordered table-striped table-condensed">
        <thead>
          <tr>
            <th>#</th>
            <th>Added</th>
            <th>Test</th>
            <th>Method</th>
            <th>Δ</th>
          </tr>
        </thead>
        <tbody>
          <tr *ngFor="let testresult of testresults">
            <td>{{ testresult.id }}</td>
            <td>{{ testresult.ctime }}</td>
            <td>{{ testresult.test.name }}</td>
            <td>
              <a [routerLink]="['/periodictable', reference_method, testresult.method.id]">
                {{ testresult.method.id }} ({{ testresult.method.pseudopotential }})
              </a>
            </td>
            <td>{{ testresult.data.deltavalue }}</td>
          </tr>
        </tbody>
      </table>
    </div>
  </div>
  `,
})

export class DeltatestResultsComponent implements OnInit {
  errorMessage: string;
  testresults: TestResult[] = [];
  reference_method: number = 3;
  methods: Object[] = [];

  constructor(private _deltatestresultsServce: DeltatestResultsService) { }

  ngOnInit() {
    this.getTestResults();
    this.getMethods();
  }

  getTestResults() {
    this._deltatestresultsServce.getTestResults(this.reference_method).subscribe(
      testresults => this.testresults = testresults,
      error => this.errorMessage = <any>error,
    );
  }

  getMethods() {
    this._deltatestresultsServce.getMethods().subscribe(
      methods => this.methods = methods,
      error => this.errorMessage = <any>error,
    );
  }

  onSelect() {
    // TODO: only update the deltatest-value
    this.testresults = [];
    this.getTestResults();
  }
}
//  vim: set ts=2 sw=2 tw=0 :
