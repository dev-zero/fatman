
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

      <h3>
        <a class="btn btn-default" (click)="gotoCalculations()">
          <i class="fa fa-arrow-left" aria-hidden="true"></i>
        </a>
        Calculation {{ calculation.id }}
      </h3>

      <div class="panel-group" id="taskAccordion" role="tablist" aria-multiselectable="true">
        <div class="panel panel-default"
          *ngFor="let task of calculation.tasks; let isFirst = first"
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
          <div class="panel-collapse collapse" role="tabpanel" 
              [class.in]="isFirst"
              [id]="'collapse_' + task.id"
              [attr.aria-labelledby]="'heading_' + task.id">
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

      <div class="panel panel-default">
        <div class="panel-heading">
          <h4 class="panel-title">Results</h4>
        </div>
        <div class="panel-body" *ngIf="!calculation.results">
          <i class="fa fa-clock-o" aria-hidden="true"></i> not yet available
        </div>
        <div class="panel-body" *ngIf="calculation.results">
          <dl class="dl-horizontal">
            <template ngFor let-result [ngForOf]="calculation.results | mapToIterable">
              <dt>{{ result.key }}</dt>
              <dd>{{ result.value }}</dd>
            </template>
          </dl>
        </div>
        <div class="list-group" *ngIf="calculation.testresults">
          <template ngFor let-testresult [ngForOf]="calculation.testresults">
            <div class="list-group-item">
              <h4 class="list-group-item-heading">Test: {{ testresult.test }}</h4>
              <dl class="list-group-item-text dl-horizontal">
                <template ngFor let-result [ngForOf]="testresult.data | mapToIterable">
                  <dt>{{ result.key }}</dt>
                  <dd>{{ result.value | stringifydata }}</dd>
                </template>
              </dl>
            </div>
          </template>
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
