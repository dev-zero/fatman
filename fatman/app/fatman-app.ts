import {Component, OnInit} from '@angular/core';
import {Routes, Router, ROUTER_DIRECTIVES, ROUTER_PROVIDERS} from '@angular/router';

import {Home} from './components/home/home';
import {Tasks} from './components/tasks/tasks';
import {About} from './components/about/about';

@Component({
  selector: 'fatman-app',
  providers: [ROUTER_PROVIDERS],
  directives: [ROUTER_DIRECTIVES],
  templateUrl: 'fatman-app.html',
})

@Routes([
  { path: '/home',  component: Home },
  { path: '/tasks', component: Tasks },
  { path: '/about', component: About },
])

export class FatmanApp implements OnInit {
  constructor(private _router: Router) {}

  ngOnInit() {
    this._router.navigate(['/home']);
  }
}
