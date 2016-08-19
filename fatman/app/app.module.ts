import { NgModule }       from '@angular/core';
import { BrowserModule }  from '@angular/platform-browser';
import { FormsModule }    from '@angular/forms';
import { HttpModule }     from '@angular/http';

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
import {PseudoFamilyList}         from './pseudos/pseudofamily-list.component';

import { ReportService } from './reports/reports.service';
import { TaskService }   from './tasks/tasks.service';
import { PseudoService } from './pseudos/pseudo.service';

@NgModule({
  imports: [
    BrowserModule,
    FormsModule,
    routing,
    HttpModule
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
    PseudoFamilyList
  ],
  providers: [
    ReportService,
    TaskService,
    PseudoService,
  ],
  bootstrap: [ AppComponent ]
})
export class AppModule {
}
//  vim: set ts=2 sw=2 tw=0 :
