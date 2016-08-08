import { Http, Response, URLSearchParams } from '@angular/http';
import { Injectable } from '@angular/core';

import { Observable } from 'rxjs/Observable';

@Injectable()
export class TaskService {

  constructor(private _http: Http) {
  }

  private _tasksUrl = '../tasks';
  private _statsUrl = '../stats/tasks';

  getTasks(statusName: string = ""): Observable<Object[]> {
    let params: URLSearchParams = new URLSearchParams();
    params.set('limit', '100');
    params.set('timeorder', 'True');

    if (statusName != "")
      params.set('status', statusName)

    return this._http.get(this._tasksUrl, {search: params})
      .map((response) => response.json())
      .catch(this.handleError);
  }

  getStats(): Observable<Object[]> {
    return this._http.get(this._statsUrl)
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
