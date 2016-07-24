import { Http, Response, RequestOptions, URLSearchParams } from '@angular/http';
import { Injectable } from '@angular/core';

import { Observable } from 'rxjs/Observable';

@Injectable()
export class PeriodictableService {

  constructor(private _http: Http) {
  }

  private _periodictableUrl = '../compare';
  private _methodlistUrl = '../methods';

  elements: Object[];

  getPeriodictable(method1, method2): Observable<Object[]> {

    let myoptions = new RequestOptions ({
      search: new URLSearchParams(`method1=${method1}&method2=${method2}`)
    });

    return this._http.get(this._periodictableUrl, myoptions)
       .map((response) => response.json())
       .catch(this.handleError);
  }

  getMethodlist(): Observable<Object[]> {
    return this._http.get(this._methodlistUrl)
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
