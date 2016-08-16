import {Component} from '@angular/core';
import { ROUTER_DIRECTIVES } from '@angular/router';

// Add the RxJS Observable operators we need in this app.
import './rxjs-operators';

@Component({
  selector: 'fatman-app',
  directives: [ROUTER_DIRECTIVES],
  template: `
    <nav class="navbar navbar-inverse navbar-static-top">
      <div class="container-fluid">
        <div class="navbar-header">
          <a class="navbar-brand" href="#">FATMAN</a>
        </div>
        <div id="navbar" class="collapse navbar-collapse">
          <ul class="nav navbar-nav">
            <li><a [routerLink]="['/home']">Home</a></li>
            <li><a [routerLink]="['/tasks']">Tasks</a></li>
            <li><a [routerLink]="['/reports']">Reports</a></li>
            <li><a [routerLink]="['/machinestatus']">Machine Status</a></li>
            <li><a [routerLink]="['/periodictable/3/72']">Periodic Table</a></li>
            <li><a [routerLink]="['/methodlist']">Method Overview</a></li>
            <li><a [routerLink]="['/pseudos']">Pseudos</a></li>
          </ul>
        </div> <!--/.navbar-collapse -->
      </div>
    </nav>
    <div class="container">
      <router-outlet></router-outlet>
    </div>
  `,
})

export class FatmanApp {
  title = 'F.A.T.M.A.N.';
}
