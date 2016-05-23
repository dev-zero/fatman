import { Http, Response } from '@angular/http';
import { Injectable } from '@angular/core';

@Injectable()
export class ReportService {

  constructor(private _http: Http) {
  }

  getMethods() {
    return this._http.get('../methods')
      .map((response) => response.json())
      .toPromise();
  }
}
//  vim: set ts=2 sw=2 tw=0 :
