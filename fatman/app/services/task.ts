import { Http, Response, Request } from 'angular2/http';
import { Injectable } from 'angular2/core';
import { Observable } from 'rxjs/Observable';
import { Observer } from 'rxjs/Observer';

// @Injectable() - normally for Dart but since we just have to have metadata
// on our object to make it inject Http into the constructor
// we use this annotation. We also could have done an
// @Inject(http) in the parameter as an alternative, but this
// is a bit drier and works for multiple constructor DI args.
// see http://blog.thoughtram.io/angular/2015/09/17/resolve-service-dependencies-in-angular-2.html
@Injectable()
export class TaskService {
  tasks$: Observable<Array<Object>>;
  private _tasksObserver: Observer<Array<Object>>;
  private _dataStore: {
    tasks: Array<Object>
  };

  constructor(private _http: Http) {
    this.tasks$ = new Observable(observer => this._tasksObserver = observer)
                                 .share();
    this._dataStore = { tasks: [] };

    console.log('Task Service created.', _http);
  }

  load() {
    this._http.get('../tasks')
      .map((response: Response) => response.json())
      .subscribe(data => {
        this._dataStore.tasks = data;
        this._tasksObserver.next(this._dataStore.tasks);
      }, error => console.log('Could not load tasks.'));
  }
}
//  vim: set ts=2 sw=2 tw=0 :
