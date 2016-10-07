import { Injectable } from '@angular/core';
import { Http, Response, RequestOptions, URLSearchParams } from '@angular/http';

import { Observable } from 'rxjs/Observable';

@Injectable()
export class DetailsService {
  constructor(private _http: Http) {
  }

  private _methodsUrl = '../api/v1/methods';
  private _testlistUrl = '../api/v1/tests';
  private _resultsUrl = '../api/v1/results';
  private _testresultUrl = '../api/v1/testresults';
  private _compareUrl = '../api/v1/compare';

  getComparisonOneTest(method1, method2, test): Observable<Object[]> {
    let myoptions = new RequestOptions ({
      search: new URLSearchParams(`method1=${method1}&method2=${method2}&test=${test}`)
    });

    return this._http.get(this._compareUrl, myoptions)
       .map((response) => response.json())
       .catch(this.handleError);
  }

  getMethods(test): Observable<Object[]> {
    let myoptions = new RequestOptions ({
      search: new URLSearchParams(`test=${test}`)
    });
    return this._http.get(this._methodsUrl, myoptions)
      .map((response) => response.json())
      .catch(this.handleError);
  }

  getMethodDetails(method_id): Observable<Object[]> {
    return this._http.get(this._methodsUrl + '/' + method_id)
      .map((response) => response.json())
      .catch(this.handleError);
  }

  getTests(): Observable<Object[]> {
    return this._http.get(this._testlistUrl)
      .map((response) => response.json())
      .catch(this.handleError);
  }

  getResults(method, test): Observable<Object[]> {
    let myoptions = new RequestOptions ({
      search: new URLSearchParams(`test=${test}&method=${method}`)
    });
    return this._http.get(this._resultsUrl, myoptions)
      .map((response) => response.json())
      .catch(this.handleError);
  }

  getTestresults(method, test): Observable<Object[]> {
    let myoptions = new RequestOptions ({
      search: new URLSearchParams(`test=${test}&method=${method}`)
    });
    return this._http.get(this._testresultUrl, myoptions)
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
