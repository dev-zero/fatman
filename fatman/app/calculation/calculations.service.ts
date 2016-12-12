import { Http, Headers, Response, URLSearchParams } from '@angular/http';
import { Injectable } from '@angular/core';

import { Observable } from 'rxjs/Observable';

import { Calculation, toCalculation } from './calculation';

@Injectable()
export class CalculationsService {
  constructor(private _http: Http) {
  }

  private _calculationsUrl = '../api/v2/calculations';
  private _calculations: Calculation[];
  private _observable: Observable<any>;

  getCalculations() : Observable<Calculation[]> {
    if (this._calculations) {
      return Observable.of(this._calculations);
    } else if (this._observable) {
      return this._observable;
    } else {
      let headers = new Headers();
      headers.append('Content-Type', 'application/json');
      this._observable = this._http.get(this._calculationsUrl, {headers: headers})
        .map(response => this._mapCalculations(response))
        .catch(this.handleError)
        .share();

      return this._observable;
    }
  }

  public getCalculation(id: string) : Observable<Calculation> {
    return Observable.of(this._calculations.find(calc => calc.id == id));
  }

  private _mapCalculations(response: Response) : Calculation[] {
    this._observable = null;
    this._calculations = response.json().map(toCalculation);
    return this._calculations;
  }

  private handleError(error: any) {
    let errMsg = (error.message) ? error.message :
      error.status ? `${error.status} - ${error.statusText}` : 'Server error';
    console.error(errMsg);
    return Observable.throw(errMsg);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
