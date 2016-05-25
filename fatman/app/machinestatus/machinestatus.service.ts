import { Http, Response } from '@angular/http';
import { Injectable } from '@angular/core';

@Injectable()
export class MachinestatusService {

  constructor(private _http: Http) {
  }

  getMachinestatus() {
    return this._http.get('../machinestatus')
      .map((response) => response.json())
      .toPromise();
  }
}
//  vim: set ts=2 sw=2 tw=0 :
