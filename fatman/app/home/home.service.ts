import { Injectable } from '@angular/core';
import { Http, Response, RequestOptions, URLSearchParams } from '@angular/http';

import { Observable } from 'rxjs/Observable';

@Injectable()
export class HomeService {
  constructor(private _http: Http) {
  }

  private _methodsUrl = '../methods';
  private _compareUrl = '../compare';

  getMethods(): Observable<Object[]> {
    return this._http.get(this._methodsUrl)
      .map((response) => response.json())
      .catch(this.handleError);
  }

  getResults(method1, method2): Observable<Object[]> {
    let myoptions = new RequestOptions ({
      search: new URLSearchParams(`method1=${method1}&method2=${method2}`)
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
