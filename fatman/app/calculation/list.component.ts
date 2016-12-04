import { Component, OnInit } from '@angular/core';

import { CalculationsService } from './calculations.service';
import { Calculation } from './calculation';

@Component({
  template: `
  <div class="row">
    <div class="col-md-12">
      <table class="table table-bordered table-striped table-condensed">
        <thead>
          <tr>
            <th>ID</th>
            <th>Test</th>
            <th>Structure</th>
            <th>Code</th>
            <th>Collection</th>
          </tr>
        </thead>
        <tbody>
          <tr *ngFor="let calculation of calculations">
            <td>{{ calculation.id }}</td>
            <td>{{ calculation.test }}</td>
            <td>{{ calculation.structure }}</td>
            <td>{{ calculation.code }}</td>
            <td>{{ calculation.collection }}</td>
          </tr>
        </tbody>
      </table>
    </div>
  </div>
  `,
})

export class CalculationListComponent implements OnInit {
  errorMessage: string;
  calculations: Calculation[] = [];

  constructor(private _calculationsService: CalculationsService) { }

  ngOnInit() {
    this.getCalculations();
  }

  getCalculations() {
    this._calculationsService.getCalculations().subscribe(
      calculations => this.calculations = calculations,
      error => this.errorMessage = <any>error,
    );
  }

  onSelect() {
    this.calculations = [];
    this.getCalculations();
  }
}
//  vim: set ts=2 sw=2 tw=0 :
