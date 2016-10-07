import { Injectable } from '@angular/core';
import { Http, Response, RequestOptions, URLSearchParams } from '@angular/http';

import { Observable } from 'rxjs/Observable';

@Injectable()
export class ReportService {
  constructor(private _http: Http) {
  }

  private _methodsUrl = '../api/v1/methods';
  private _compareUrl = '../api/v1/compare';
  private _testlistUrl = '../api/v1/tests';

  getMethods(): Observable<Object[]> {
    return this._http.get(this._methodsUrl)
      .map((response) => response.json())
      .catch(this.handleError);
  }

  getTests(): Observable<Object[]> {
    return this._http.get(this._testlistUrl)
      .map((response) => response.json())
      .catch(this.handleError);
  }

  getComparison(method1, method2): Observable<Object[]> {
    let myoptions = new RequestOptions ({
      search: new URLSearchParams(`method1=${method1}&method2=${method2}`)
    });

    return this._http.get(this._compareUrl, myoptions)
       .map((response) => response.json())
       .catch(this.handleError);
  }

  getElementComparison(method1, test1): Observable<Object[]> {
    let myoptions = new RequestOptions ({
      search: new URLSearchParams(`method1=${method1}&test=${test1}`)
    });

    return this._http.get(this._compareUrl, myoptions)
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
