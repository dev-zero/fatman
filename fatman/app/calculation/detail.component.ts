
import { Component, OnInit, HostBinding,
         trigger, transition, animate,
         style, state } from '@angular/core';
import { Router, ActivatedRoute, Params } from '@angular/router';

import { Calculation }         from './calculation';
import { CalculationsService } from './calculations.service';

@Component({
  template: `
  <div class="row">
    <div class="col-md-12"
      *ngIf="calculation">

      <h3>Calculation {{ calculation.id }}</h3>
      <p>
        <button (click)="gotoCalculations()">Back</button>
      </p>

      <div class="panel-group" id="taskAccordion" role="tablist" aria-multiselectable="true">
        <div class="panel panel-default"
          *ngFor="let task of calculation.tasks"
          [ngClass]="{
            'panel-danger': task.status == 'error',
            'panel-warning': task.status == 'running',
            'panel-success': task.status == 'done'
            }">
          <div class="panel-heading" role="tab"
            [id]="'heading_' + task.id">
            <h4 class="panel-title">
              <a role="button" data-toggle="collapse" data-parent="#taskAccordion"
                [href]="'#collapse_' + task.id"
                [attr.aria-controls]="'collapse_' + task.id"
                [attr.aria-expanded]="isFirst">
                Task {{ task.id }} <span class="pull-right">{{ task.status | uppercase }}</span>
              </a>
            </h4>
          </div>
          <div class="panel-collapse collapse in" role="tabpanel" 
              [id]="'collapse_' + task.id" [attr.aria-labelledby]="'heading_' + task.id">
            <div class="panel-body">
              <dl class="dl-horizontal">
                <dt>Last Update</dt><dd>{{ task.mtime }}</dd>
                <dt>Created</dt><dd>{{ task.ctime }}</dd>
                <dt>Machine</dt><dd>{{ task.machine }}</dd>
              </dl>
            </div>
          </div>
        </div>
      </div>
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
