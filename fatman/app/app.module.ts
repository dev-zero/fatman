import { NgModule }       from '@angular/core';
import { BrowserModule }  from '@angular/platform-browser';
import { FormsModule }    from '@angular/forms';
import { HttpModule }     from '@angular/http';

import { AppComponent }   from './app.component';
import { routing }        from './app.routing';

import { Home }                     from './home/home';
import { TasksStats }               from './tasks/tasks-stats.component';
import { TasksList }                from './tasks/tasks-list.component';
import { MachineStatus }            from './machinestatus/machinestatus-list.component';
import { Reports }                  from './reports/reports-list.component';
import { ReportsComparison }        from './reports/reports-comparison.component';
import { ReportsElementComparison } from './reports/reports-elementcomparison.component';
import { Periodictable }            from './periodictable/periodictable-list.component';
import { MethodList }               from './methods/methods-list.component';
import { PseudoFamilyList }         from './pseudos/pseudofamily-list.component';
import { PseudoTable }              from './pseudos/pseudo-table.component';
import { PseudoTableElement }       from './pseudos/pseudo-table-element.component';
import { Details }                  from './details/details-list.component';
import { MethoddetailsComponent }   from './details/methoddetails.component';
import { TestdetailsComponent }     from './details/testdetails.component';
import { ComparematrixComponent }   from './details/comparematrix.component';

import { ReportService } from './reports/reports.service';
import { TaskService }   from './tasks/tasks.service';
import { PseudoService } from './pseudos/pseudo.service';

import { TruncatePipe }                           from './details/testdetails.component';
import { DeltavaluePipe }                         from './periodictable/periodictable-list.component';
import { IterablePipe, SettingsPipe, ConcatPipe } from './common/pipes';

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
    TasksList,
    MachineStatus,
    Reports,
    ReportsComparison,
    ReportsElementComparison,
    Periodictable,
    Details,
    MethoddetailsComponent,
    TestdetailsComponent,
    ComparematrixComponent,
    MethodList,
    PseudoFamilyList,
    PseudoTable,
    PseudoTableElement,
    TruncatePipe,
    IterablePipe,
    SettingsPipe,
    ConcatPipe,
    DeltavaluePipe
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
