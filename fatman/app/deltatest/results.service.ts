import { Http, Response, URLSearchParams } from '@angular/http';
import { Injectable } from '@angular/core';

import { Observable } from 'rxjs/Observable';

import { TestResult, toTestResult } from './testresult';

@Injectable()
export class DeltatestResultsService {

  constructor(private _http: Http) {
  }

  private _testresultsUrl = '../api/v1/testresults';
  private _methodsUrl = '../api/v1/methods';

  getTestResults(ref_method: string = null) : Observable<TestResult[]> {
    let params: URLSearchParams = new URLSearchParams();

    params.set('test', 'deltatest')

    if (!!ref_method)
      params.set('deltaref', ref_method.toString())

    return this._http.get(this._testresultsUrl, {search: params})
      .map(this._mapTestResults)
      .catch(this.handleError);
  }

  private _mapTestResults(response: Response) : TestResult[] {
    return response.json().map(toTestResult);
  }

  private handleError(error: any) {
    let errMsg = (error.message) ? error.message :
      error.status ? `${error.status} - ${error.statusText}` : 'Server error';
    console.error(errMsg);
    return Observable.throw(errMsg);
  }


  getMethods(): Observable<Object[]> {
    return this._http.get(this._methodsUrl)
      .map((response) => response.json())
      .catch(this.handleError);
  }

}
//  vim: set ts=2 sw=2 tw=0 :
