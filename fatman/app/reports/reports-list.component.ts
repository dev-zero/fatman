import {Component, OnInit} from '@angular/core';
import {Router} from '@angular/router';

import {ReportService} from './reports.service';

@Component({
   template: `
    <table class="table table-bordered table-striped table-condensed table-hover">
      <thead>
        <th>&nbsp;</th>
          <th *ngFor="let method of methods">{{method['id']}}</th>
      </thead>
      <tbody>
        <tr *ngFor="let m_row of methods">
          <th scope="row">{{m_row['id']}}</th>
            <td *ngFor="let m_col of methods">
              <div *ngIf="m_row.id < m_col.id"
                (click)="onSelect(m_row, m_col)">
                  {{m_row['id']}}/{{m_col['id']}}
              </div>
            </td>
        </tr>
      </tbody>
    </table>
  `,
})

export class Reports implements OnInit {
  title = 'Reports';

  errorMessage: string;
  methods: Object[];
  mode = 'Observable';

  constructor(
    private _service: ReportService,
    private _router: Router) { }

  ngOnInit() {
    this._service.getMethods().subscribe(
      methods => this.methods = methods,
      error => this.errorMessage = <any>error);
  }

  onSelect(m_row, m_col) {
    this._router.navigate(['/reports/comparison', m_row.id, m_col.id]);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
