
import { Component, OnInit } from '@angular/core';
import { Router, ActivatedRoute, Params } from '@angular/router';

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
            <th>Created</th>
            <th>Modified</th>
            <th>Status</th>
            <th>Result</th>
          </tr>
        </thead>
        <tbody>
          <tr *ngFor="let calculation of calculations"
            (click)="onSelect(calculation)">
            <td>{{ calculation.id }}</td>
            <td>{{ calculation.test }}</td>
            <td>{{ calculation.structure }}</td>
            <td>{{ calculation.code }}</td>
            <td>{{ calculation.collection }}</td>
            <td>{{ calculation.current_task.ctime | date:'medium' }}</td>
            <td>{{ calculation.current_task.mtime | date:'medium' }}</td>
            <td>{{ calculation.current_task.status }}</td>
            <td><i class="fa" aria-hidden="true"
                 [class.fa-check]="calculation.result_available"
                 [class.fa-times]="!calculation.result_available"></i>
            </td>
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

  private _selectedId: string;

  constructor(
    private _calculationsService: CalculationsService,
    private _route: ActivatedRoute,
    private _router: Router
  ) {}

  ngOnInit() {
    this.getCalculations();
  }

  getCalculations() {
    this._calculationsService.getCalculations().subscribe(
      calculations => this.calculations = calculations,
      error => this.errorMessage = <any>error,
    );
  }

  onSelect(calculation: Calculation) {
    this._router.navigate(['/calculations', calculation.id]);
  }
}
