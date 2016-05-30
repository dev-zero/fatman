import {Component} from '@angular/core';

// Add the RxJS Observable operators we need in this app.
import './rxjs-operators';

import {Routes, Router, ROUTER_DIRECTIVES} from '@angular/router';

import {Home} from './home/home';
import {Tasks} from './tasks/tasks-list.component';
import {MachineStatus} from './machinestatus/machinestatus-list.component';
import {Reports} from './reports/reports-list.component';
import {ReportsComparison} from './reports/reports-comparison.component';
import {Periodictable} from './periodictable/periodictable-list.component';

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
            <li><a [routerLink]="['/machinestatus']">Machine Status</a></li>
            <li><a [routerLink]="['/periodictable']">Periodic Table</a></li>
          </ul>
        </div> <!--/.navbar-collapse -->
      </div>
    </nav>
    <div class="container">
      <router-outlet></router-outlet>
    </div>
  `,
})

@Routes([
  { path: '/', component: Home },
  { path: '/home', component: Home },
  { path: '/tasks', component: Tasks },
  { path: '/reports/comparison/:id1/:id2', component: ReportsComparison },
  { path: '/reports', component: Reports },
  { path: '/machinestatus', component: MachineStatus },
  { path: '/periodictable', component: Periodictable },
])

export class FatmanApp {}
