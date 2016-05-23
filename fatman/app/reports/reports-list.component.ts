import {Component} from '@angular/core';
import {OnActivate, Router, RouteSegment, RouteTree} from '@angular/router';

import {ReportService} from './reports.service';

@Component({
  selector: 'reports',
  template: `
    <table class="table table-bordered table-striped table-condensed">
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
    </table>`,
  providers: [ReportService],
})

export class Reports implements OnActivate {
  title = 'Reports';

  methods: Object[];

  private _selectedRowId: number;
  private _selectedColId: number;

  constructor(
    private _service: ReportService,
    private _router: Router) { }

  routerOnActivate(curr: RouteSegment,
                   prev?: RouteSegment,
                   currTree?: RouteTree,
                   prevTree?: RouteTree): void {
    this._selectedRowId = +curr.getParam('m_row_id');
    this._selectedColId = +curr.getParam('m_col_id');

    this._service.getMethods().then(methods => this.methods = methods);
  }

  isSelected(m_row, m_col) {
    return (m_row.id === this._selectedRowId) && (m_col.id === this._selectedColId);
  }

  onSelect(m_row, m_col) {
    this._router.navigate(['/report', m_row.id, m_col.id]);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
