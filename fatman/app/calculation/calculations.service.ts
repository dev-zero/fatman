import { Http, Response, URLSearchParams } from '@angular/http';
import { Injectable } from '@angular/core';

import { Observable } from 'rxjs/Observable';

import { Calculation, toCalculation } from './calculation';

@Injectable()
export class CalculationsService {

  constructor(private _http: Http) {
  }

  private _calculationsUrl = '../api/v2/calculations';

  getCalculations(ref_method: string = null) : Observable<Calculation[]> {
    return this._http.get(this._calculationsUrl)
      .map(this._mapCalculations)
      .catch(this.handleError);
  }

  private _mapCalculations(response: Response) : Calculation[] {
    return response.json().map(toCalculation);
  }

  private handleError(error: any) {
    let errMsg = (error.message) ? error.message :
      error.status ? `${error.status} - ${error.statusText}` : 'Server error';
    console.error(errMsg);
    return Observable.throw(errMsg);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
