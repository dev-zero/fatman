import { Http, Headers, Response, URLSearchParams } from '@angular/http';
import { Injectable } from '@angular/core';

import { Observable } from 'rxjs/Observable';

import { Task2 } from './task';

@Injectable()
export class Task2Service {
  constructor(private _http: Http) {
  }

  private _tasksUrl = '../api/v2/tasks';

  public getTask(id: string) : Observable<Task2> {
    let headers = new Headers();
    headers.append('Content-Type', 'application/json');
    return this._http.get(this._tasksUrl + '/' + id, {headers: headers})
      .map(response => response.json() as Task2)
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
