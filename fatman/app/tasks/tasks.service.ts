import { Http, Response } from '@angular/http';
import { Injectable } from '@angular/core';

@Injectable()
export class TaskService {

  constructor(private _http: Http) {
  }

  getTasks() {
    return this._http.get('../tasks')
      .map((response) => response.json())
      .toPromise();
  }
}
//  vim: set ts=2 sw=2 tw=0 :
