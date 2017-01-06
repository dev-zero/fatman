import { NgModule }       from '@angular/core';
import { BrowserModule }  from '@angular/platform-browser';
import { FormsModule }    from '@angular/forms';
import { HttpModule }     from '@angular/http';

import { ChartsModule }   from 'ng2-charts/ng2-charts';

import { AppComponent }   from './app.component';
import { routing }        from './app.routing';

import { Home }                       from './home/home';
import { TasksStats }                 from './tasks/tasks-stats.component';
import { TasksList }                  from './tasks/tasks-list.component';
import { MachineStatus }              from './machinestatus/machinestatus-list.component';
import { Reports }                    from './reports/reports-list.component';
import { ReportsComparison }          from './reports/reports-comparison.component';
import { ReportsElementComparison }   from './reports/reports-elementcomparison.component';
import { Periodictable }              from './periodictable/periodictable-list.component';
import { MethodList }                 from './methods/methods-list.component';
import { PseudoFamilyList }           from './pseudos/pseudofamily-list.component';
import { PseudoTable }                from './pseudos/pseudo-table.component';
import { PseudoTableElement }         from './pseudos/pseudo-table-element.component';
import { Details }                    from './details/details-list.component';
import { MethoddetailsComponent }     from './details/methoddetails.component';
import { TestdetailsComponent }       from './details/testdetails.component';
import { ComparematrixComponent }     from './details/comparematrix.component';
import { DeltatestResultsComponent }  from './deltatest/results.component';
import { CalculationListComponent }   from './calculation/list.component';
import { CalculationDetailComponent } from './calculation/detail.component';
import { Task2DetailComponent }       from './task2/detail.component';

import { ReportService }           from './reports/reports.service';
import { TaskService }             from './tasks/tasks.service';
import { PseudoService }           from './pseudos/pseudo.service';
import { DeltatestResultsService } from './deltatest/results.service';
import { CalculationsService }     from './calculation/calculations.service';
import { Task2Service }            from './task2/task.service';

import { TruncatePipe }                           from './details/testdetails.component';
import { DeltavaluePipe }                         from './periodictable/periodictable-list.component';
import { IterablePipe, SettingsPipe, ConcatPipe, StringifyPipe } from './common/pipes';

@NgModule({
  imports: [
    BrowserModule,
    FormsModule,
    routing,
    HttpModule,
    ChartsModule
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
    DeltatestResultsComponent,
    CalculationListComponent,
    CalculationDetailComponent,
    Task2DetailComponent,
    TruncatePipe,
    IterablePipe,
    SettingsPipe,
    ConcatPipe,
    DeltavaluePipe,
    StringifyPipe
  ],
  providers: [
    ReportService,
    TaskService,
    PseudoService,
    DeltatestResultsService,
    CalculationsService,
    Task2Service
  ],
  bootstrap: [ AppComponent ]
})
export class AppModule {
}
//  vim: set ts=2 sw=2 tw=0 :
