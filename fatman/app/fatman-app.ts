import {Component} from 'angular2/core';
import {RouteConfig, ROUTER_DIRECTIVES, ROUTER_PROVIDERS} from 'angular2/router';

import {TaskService} from './services/task';

import {Home} from './components/home/home';
import {Tasks} from './components/tasks/tasks';
import {About} from './components/about/about';

@Component({
  selector: 'fatman-app',
  providers: [
    ROUTER_PROVIDERS,
    TaskService
  ],
  templateUrl: 'fatman-app.html',
  directives: [ROUTER_DIRECTIVES],
})
@RouteConfig([
  { path: '/home', component: Home, name: 'Home', useAsDefault: true},
  { path: '/tasks', component: Tasks, name: 'Tasks'},
  { path: '/about', component: About, name: 'About'},
])
export class FatmanApp {

  constructor() {}

}
