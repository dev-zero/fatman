import {Component} from '@angular/core';
import {OnActivate, Router, RouteSegment} from '@angular/router';

import {ReportService} from './reports.service';

@Component({
  template: `
    <table class="table table-bordered table-striped table-condensed">
      <thead>
        <th>z</th>
        <th>Element</th>
        <th>V0</th>
        <th>B0</th>
        <th>B1</th>
        <th>V0,r</th>
        <th>B0,r</th>
        <th>B1,r</th>
      </thead>
      <tbody>
      </tbody>
    </table>
  `,
  providers: [ReportService],
})

export class ReportsComparison implements OnActivate {

  constructor(
    private _service: ReportService,
    private _router: Router) { }

  routerOnActivate(curr: RouteSegment): void {
    let id1 = +curr.getParam('id1');
    let id2 = +curr.getParam('id2');

    // TODO
  }

  gotoReports() {
    this._router.navigate(['/reports']);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
