import { Http, Response } from '@angular/http';
import { Injectable } from '@angular/core';

import { Observable } from 'rxjs/Observable';

@Injectable()
export class PeriodictableService {

  constructor(private _http: Http) {
  }

  private _periodictableUrl = '../compare?method1=3&method2=72';

  getPeriodictable(): Observable<Object[]> {
    return this._http.get(this._periodictableUrl)
      .map((response) => response.json())
      .catch(this.handleError);
  }

  private handleError(error: any) {
    let errMsg = (error.message) ? error.message :
      error.status ? `${error.status} - ${error.statusText}` : 'Server error';
    console.error(errMsg);
    return Observable.throw(errMsg);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
