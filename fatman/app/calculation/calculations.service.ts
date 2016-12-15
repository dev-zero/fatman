import { Http, Headers, Response, URLSearchParams } from '@angular/http';
import { Injectable } from '@angular/core';

import { Observable } from 'rxjs/Observable';

import { Calculation } from './calculation';

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
    let headers = new Headers();
    headers.append('Content-Type', 'application/json');
    return this._http.get(this._calculationsUrl + '/' + id, {headers: headers})
      .map(response => this._mapCalculation(response))
      .catch(this.handleError);
  }

  private _mapCalculations(response: Response) : Calculation[] {
    this._observable = null;
    this._calculations = response.json().map(calc => calc as Calculation);
    return this._calculations;
  }

  private _mapCalculation(response: Response) : Calculation {
    let calculation = response.json() as Calculation;
    // sort by inverse mtime:
    calculation.tasks = calculation.tasks.sort((t1,t2) => Date.parse(t2['mtime']) - Date.parse(t1['mtime']));
    return calculation;
  }

  private handleError(error: any) {
    let errMsg = (error.message) ? error.message :
      error.status ? `${error.status} - ${error.statusText}` : 'Server error';
    console.error(errMsg);
    return Observable.throw(errMsg);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
