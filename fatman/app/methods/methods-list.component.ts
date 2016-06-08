import {Component, OnInit} from '@angular/core';

import {PeriodictableService} from '../periodictable/periodictable.service';
import {MethoddetailsComponent} from '../details/methoddetails.component';

@Component({
  directives: [MethoddetailsComponent],
  template: `
     <div class="form-group container" style="width: 1100px; height: 900px; overflow-y: scroll; border: 2px solid">
        <methoddetails *ngFor="let method of methods" method_id="{{ method.id }}"></methoddetails>

     </div>`,
  providers: [PeriodictableService]
})

export class MethodList implements OnInit {
  title = 'Machines';

  errorMessage: string;
  methods: Object[];

  constructor(private _service: PeriodictableService) { }

  getMethods() {
    this._service.getMethodlist().subscribe(
      methods => this.methods = methods,
      error => this.errorMessage = <any>error);
  }

  ngOnInit() {
    this.getMethods();
  }
}
