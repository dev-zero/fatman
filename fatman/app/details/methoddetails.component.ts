import { Component, Input, OnChanges, SimpleChange } from '@angular/core';

import { DetailsService } from './details.service';

@Component({
  selector:'methoddetails', 
  providers: [DetailsService],
  template: `
    <div *ngIf="method && !small" class="panel panel-default">
      <div class="panel-body">
        <div class="text-center pull-left" style="margin-right: 1em;">
          Method <p class="lead" style="font-size: 200%;"><strong>{{ method.id }}</strong></p>
        </div>
        <ul class="list-inline">
          <li><strong>Code:</strong> {{ method.code }}</li>
          <li><strong>Pseudopotential:</strong> {{ method.pseudopotential }}</li>
          <li><strong>Basis Set:</strong> {{ method.basis_set }}</li>
        </ul>
        <ul class="list-inline">
          <li *ngFor="let settingpair of method.settings | mapToIterable">
            {{settingpair.key}}: {{ settingpair.value | stringifysetting }}
          </li>
        </ul>
      </div>
    </div>
    <span *ngIf="method && small" title="{{ method.settings | mapToIterable | concatkvpairlist }}">
      <span title="{{ method.settings | mapToIterable | concatkvpairlist }}">
        <span class="text-muted">{{ method.id }}</span>
        <span class="small">{{ method.basis_set }} / {{ method.pseudopotential }} ({{ method.code }})</span>
      </span>
    </span>
  `,
})
export class MethoddetailsComponent implements OnChanges {
  @Input()
  method_id: number;
  @Input()
  small: boolean = false;

  constructor(private _service: DetailsService) { }

  errorMessage: string;
  method: Object[];


  getMethodDetails(method_id) {
    this._service.getMethodDetails(method_id).subscribe(
      method => this.method = method,
      error => this.errorMessage = <any>error);
  };

  ngOnChanges(changes: {[propertyName: string]: SimpleChange}) {
    if (changes['method_id'])
      this.getMethodDetails(this.method_id);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
