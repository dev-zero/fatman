import { NgModule }       from '@angular/core';
import { BrowserModule }  from '@angular/platform-browser';
import { FormsModule }    from '@angular/forms';
import { HTTP_PROVIDERS } from '@angular/http'; // make this provider global

import { AppComponent }   from './app.component';
import { routing }        from './app.routing';

import {Home}                     from './home/home';
import {TasksStats}               from './tasks/tasks-stats.component';
import {MachineStatus}            from './machinestatus/machinestatus-list.component';
import {Reports}                  from './reports/reports-list.component';
import {ReportsComparison}        from './reports/reports-comparison.component';
import {ReportsElementComparison} from './reports/reports-elementcomparison.component';
import {Periodictable}            from './periodictable/periodictable-list.component';
import {Details}                  from './details/details-list.component';
import {MethodList}               from './methods/methods-list.component';
import {PseudoList}               from './pseudos/pseudo-list.component';

import { ReportService } from './reports/reports.service';
import { TaskService }   from './tasks/tasks.service';
import { PseudoService } from './pseudos/pseudo.service';

@NgModule({
  imports: [
    BrowserModule,
    FormsModule,
    routing
  ],
  declarations: [
    AppComponent,
    Home,
    TasksStats,
    MachineStatus,
    Reports,
    ReportsComparison,
    ReportsElementComparison,
    Periodictable,
    Details,
    MethodList,
    PseudoList
  ],
  providers: [
    ReportService,
    TaskService,
    PseudoService,
    HTTP_PROVIDERS
  ],
  bootstrap: [ AppComponent ]
})
export class AppModule {
}
//  vim: set ts=2 sw=2 tw=0 :
