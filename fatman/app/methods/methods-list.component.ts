import { Component, OnInit } from '@angular/core';

import { PeriodictableService } from '../periodictable/periodictable.service';

@Component({
  template: `
  <!-- <div class="form-group container" style="width: 1100px; height: 900px; overflow-y: scroll; border: 2px solid">-->
     <div class="container">
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
