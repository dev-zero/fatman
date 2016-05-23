import {Component} from '@angular/core';
import {Routes, Router, ROUTER_DIRECTIVES} from '@angular/router';

import {Home} from './home/home';
import {Tasks} from './tasks/tasks-list.component';
import {Reports} from './reports/reports-list.component';

@Component({
  selector: 'fatman-app',
  directives: [ROUTER_DIRECTIVES],
  template: `
    <nav class="navbar navbar-inverse navbar-static-top">
      <div class="container">
	<div class="navbar-header">
	  <a class="navbar-brand" href="#">FATMAN</a>
	</div>
	<div id="navbar" class="collapse navbar-collapse">
	  <ul class="nav navbar-nav">
	    <li><a [routerLink]="['/home']">Home</a></li>
	    <li><a [routerLink]="['/tasks']">Tasks</a></li>
	    <li><a [routerLink]="['/reports']">Reports</a></li>
	  </ul>
	</div>
      </div>
    </nav>
    <div>
      <router-outlet></router-outlet>
    </div>
  `,
})

@Routes([
  { path: '/', component: Home },
  { path: '/home', component: Home },
  { path: '/tasks', component: Tasks },
  { path: '/reports', component: Reports },
])

export class FatmanApp {}
