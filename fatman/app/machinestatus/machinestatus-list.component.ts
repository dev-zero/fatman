import {Component, OnInit} from '@angular/core';

import {MachinestatusService} from './machinestatus.service';

@Component({
  template: `
  <table class="machinestatus table table-bordered table-striped table-condensed">
    <thead>
    <tr style="border-bottom: 2px solid; border-top: 2px solid;">
        <th> Machine </th>
        <th> Running </th>
        <th> Total </th>
    </tr>
    </thead>
    <tbody>
    <tr class="" *ngFor="let machine of machines">
        <td *ngFor="let item of machine">
            {{item}}
        </td>
    </tr>
    </tbody>
  </table>`,
  providers: [MachinestatusService]
})

export class MachineStatus implements OnInit {
  title = 'Machines';

  errorMessage: string;
  machines: Object[];

  constructor(private _MachinestatusService: MachinestatusService) { }

  getMachinestatus() {
    this._MachinestatusService.getMachinestatus().subscribe(
      machines => this.machines = machines,
      error => this.errorMessage = <any>error);
  }

  ngOnInit() {
    this.getMachinestatus();
  }
}
