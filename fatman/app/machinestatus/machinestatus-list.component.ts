import {Component, OnInit} from '@angular/core';

import {MachinestatusService} from './machinestatus.service';

@Component({
  template: `
  <div style="position: relative; left: 30px; /*border: 1px solid #73AD21; padding:5px; display: inline-block;*/">
  <table >
    <thead>
    <tr style="border-bottom: 2px solid; border-top: 2px solid;">
        <th style="padding-right:20px"> <b>Machine</b> </th>
        <th style="padding-right:20px"> <b>Running</b> </th>
        <th style="padding-right:20px"> <b>Total</b> </th>
    </tr>
    </thead>
    <tbody>
    <tr class="" *ngFor="let machine of machines" style="border-bottom: 1px solid;">
        <td *ngFor="let item of machine">
            {{item}}
        </td>
    </tr>
    </tbody>
  </table></div>`,
  providers: [MachinestatusService]
})

export class MachineStatus implements OnInit {
  title = 'Machines';
  machines: Object[];

  constructor(private _MachinestatusService: MachinestatusService) { }

  getMachinestatus() {
    this._MachinestatusService.getMachinestatus().then(machines => this.machines = machines);
  }

  ngOnInit() {
    this.getMachinestatus();
  }
}
