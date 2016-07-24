
import {bootstrap} from '@angular/platform-browser-dynamic';
//import { enableProdMode } from '@angular/core';
import {HTTP_PROVIDERS} from '@angular/http';

import {FatmanApp} from './app.component';
import {ReportService} from './reports/reports.service';
import {TaskService} from './tasks/tasks.service';
import {PseudoService} from './pseudos/pseudo.service';
import {appRouterProviders} from './app.routes';

//enableProdMode();

// Since services are simple classes, we have to specify providers here
// instead of in some @Component section.
bootstrap(FatmanApp, [
  [ReportService, TaskService, PseudoService, HTTP_PROVIDERS, appRouterProviders],
]);

//  vim: set ts=2 sw=2 tw=0 :
