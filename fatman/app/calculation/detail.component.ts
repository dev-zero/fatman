
import { Component, OnInit, HostBinding,
         trigger, transition, animate,
         style, state } from '@angular/core';
import { Router, ActivatedRoute, Params } from '@angular/router';

import { Calculation }         from './calculation';
import { CalculationsService } from './calculations.service';
import { Task2Service }        from '../task2/task.service';
import { Task2 }               from '../task2/task';

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

      <div class="panel panel-default">
        <div class="panel-heading">
          <h4 class="panel-title">Settings</h4>
        </div>
        <div class="panel-body">
          <dl class="dl-horizontal">
            <dt>Structure</dt>
            <dd>{{ calculation.structure }}</dd>

            <dt>Code</dt>
            <dd>{{ calculation.code }}</dd>

            <dt>Collection</dt>
            <dd>{{ calculation.collection }}</dd>

            <template ngFor let-basis_set [ngForOf]="calculation.basis_sets">
              <dt>Basis Set for {{ basis_set.basis_set.element }} ({{ basis_set.type }})</dt>
              <dd>{{ basis_set.basis_set.family }}</dd>
            </template>

            <template ngFor let-pseudo [ngForOf]="calculation.pseudos">
              <dt>Pseudo</dt>
              <dd>{{ pseudo }}</dd>
            </template>
          </dl>
        </div>
      </div>

      <task2-detail [task]="calculation.current_task"></task2-detail>

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

              <div class="list-group-item-text"
               *ngIf="testresult.test == 'deltatest'">
                <dl class="dl-horizontal">
                  <dt>Status</dt>
                  <dd>{{ testresult.data.status}}</dd>
                  <template ngFor let-check [ngForOf]="testresult.data.checks | mapToIterable">
                    <dt>Check: {{ check.key }}</dt>
                    <dd>{{ check.value }}</dd>
                  </template>
                  <dt>E<sub>0</sub></dt>
                  <dd *ngIf="testresult.data.status == 'fitted'">{{ testresult.data.coefficients.E0 }}</dd>
                  <dt>B<sub>0</sub></dt>
                  <dd *ngIf="testresult.data.status == 'fitted'">{{ testresult.data.coefficients.B0 }}</dd>
                  <dt>B<sub>1</sub></dt>
                  <dd *ngIf="testresult.data.status == 'fitted'">{{ testresult.data.coefficients.B1 }}</dd>
                  <dt>V</dt>
                  <dd *ngIf="testresult.data.status == 'fitted'">{{ testresult.data.coefficients.V }}</dd>
                  <dt>R</dt>
                  <dd *ngIf="testresult.data.status == 'fitted'">{{ testresult.data.coefficients.R }}</dd>
                </dl>

                <table class="table table-bordered table-condensed">
                  <thead>
                    <tr>
                      <td>Volume</td>
                      <td>Energy</td>
                    </tr>
                  </thead>
                  <tbody>
                    <tr *ngFor="let volume of testresult.data.volumes; let idx = index">
                      <td>{{ volume }}</td>
                      <td>{{ testresult.data.energies[idx] }}</td>
                    </tr>
                  </tbody>
                </table>

                <div style="display: block;">
                  <canvas baseChart width="400" height="400"
                    [datasets]="lineChartData"
                    [labels]="lineChartLabels"
                    [chartType]="lineChartType"
                    [options]="lineChartOptions"
                    >
                  </canvas>
                </div>
              </div>

              <div class="list-group-item-text"
               *ngIf="testresult.test != 'deltatest'">
                <dl class="list-group-item-text dl-horizontal">
                  <template ngFor let-result [ngForOf]="testresult.data | mapToIterable">
                    <dt>{{ result.key }}</dt>
                    <dd>{{ result.value | stringifydata }}</dd>
                  </template>
                </dl>
              </div>

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

  public calculation: Calculation;
  public lineChartData: Array<any> = [];
  public lineChartType:string = 'line';
  public lineChartOptions:any = {
    animation: false,
    responsive: true,
	scales: {
	  xAxes: [{
		type: 'linear',
		position: 'bottom'
	  }]
	}
  };

  constructor(
    private _route: ActivatedRoute,
    private _router: Router,
    private _service: CalculationsService,
    private _tService: Task2Service
  ) {}

  ngOnInit() {
    this._route.params
      .switchMap((params: Params) => this._service.getCalculation(params['id']))
      .subscribe((calc: Calculation) => {
        this.calculation = calc;

        // fetch the full task object
        if (this.calculation.current_task)
          this._tService.getTask(this.calculation.current_task.id)
            .subscribe((full_task: Task2) => { this.calculation.current_task = full_task; });

        if (calc.testresults && calc.testresults[0].test == 'deltatest') {
          let dtres = calc.testresults[0].data;

          let E0 = 0.;
          if (dtres.status == 'fitted')
            E0 = dtres.coefficients.E0;

          this.lineChartData.push({
            data: dtres.volumes.map(function(x, i) {return {x: x, y: dtres.energies[i] - E0}; }),
            label: 'calculated',
            fill: false,
            showLine: false
          });

          if (dtres.status == 'fitted') {
            let Vmin = dtres.volumes[0];
            let Vmax = dtres.volumes[dtres.volumes.length-1];
            let V0 = dtres.coefficients.V;
            let B0 = dtres.coefficients.B0 * 1.e9 / 1.602176565e-19 / 1.e30;
            let B1 = dtres.coefficients.B1;
            let n = 100;
            let VE = Array<any>(n);
            for (let i=0; i < n; ++i) {
              let v = Math.min(0.93*V0, Vmin)*(1-i/n) + Math.max(1.07*V0, Vmax)*(i/n);

              // the calculated values are shifted by -E0 such that we get the minimum at 0
              let e = 9./16. * V0 * B0 * (
                  ((V0/v)**(2./3.) -1)**3 * B1 +
                  ((V0/v)**(2./3.) -1)**2 * (6-4*(V0/v)**(2./3.)) );

              VE[i] = {x: v, y: e};
            }

            this.lineChartData.push({
              data: VE,
              label: 'fitted',
              fill: false,
              pointRadius: 0
            });
          }
        }
      });
  }

  public gotoCalculations() {
    let calculationId = this.calculation ? this.calculation.id : null;
    // Pass along the calculation id if available
    // so that the CalculationList component can select that calculation.
    this._router.navigate(['/calculations', { id: calculationId }]);
  }
}
