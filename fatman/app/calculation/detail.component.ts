
import { Component, OnInit, HostBinding,
         trigger, transition, animate,
         style, state } from '@angular/core';
import { Router, ActivatedRoute, Params } from '@angular/router';

import { Calculation }         from './calculation';
import { CalculationsService } from './calculations.service';

@Component({
  template: `
  <div *ngIf="calculation">
    <h3>Calculation {{ calculation.id }}</h3>
    <div>
      <label>Id: </label>{{ calculation.id }}
	</div>
    <p>
      <button (click)="gotoCalculations()">Back</button>
    </p>

    <div class="panel panel-default">
      <div class="panel-heading">Tasks</div>
      <table class="table">
        <thead>
          <tr>
            <th>Created</th>
            <th>Last update</th>
            <th>Machine</th>
            <th>Status</th>
          </tr>
        </thead>
        <tbody>
          <tr *ngFor="let task of calculation.tasks">
            <td>{{ task.ctime }}</td>
            <td>{{ task.mtime }}</td>
            <td>{{ task.machine }}</td>
            <td>{{ task.status }}</td>
          </tr>
        </tbody>
      </table>
    </div>
  </div>
  `,
  animations: [
    trigger('routeAnimation', [
      state('*',
        style({
          opacity: 1,
          transform: 'translateX(0)'
        })
      ),
      transition(':enter', [
        style({
          opacity: 0,
          transform: 'translateX(-100%)'
        }),
        animate('0.2s ease-in')
      ]),
      transition(':leave', [
        animate('0.5s ease-out', style({
          opacity: 0,
          transform: 'translateY(100%)'
        }))
      ])
    ])
  ]
})
export class CalculationDetailComponent implements OnInit {
  @HostBinding('@routeAnimation') get routeAnimation() {
    return true;
  }

  @HostBinding('style.display') get display() {
    return 'block';
  }

  @HostBinding('style.position') get position() {
    return 'absolute';
  }

  calculation: Calculation;

  constructor(
    private _route: ActivatedRoute,
    private _router: Router,
    private _service: CalculationsService
  ) {}

  ngOnInit() {
    this._route.params
      .switchMap((params: Params) => this._service.getCalculation(params['id']))
      .subscribe((calculation: Calculation) => this.calculation = calculation);
  }

  gotoCalculations() {
    let calculationId = this.calculation ? this.calculation.id : null;
    // Pass along the calculation id if available
    // so that the CalculationList component can select that calculation.
    this._router.navigate(['/calculations', { id: calculationId }]);
  }
}
