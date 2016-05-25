
import {bootstrap} from '@angular/platform-browser-dynamic';
//import { enableProdMode } from '@angular/core';
import {HTTP_PROVIDERS} from '@angular/http';
import {ROUTER_PROVIDERS} from '@angular/router';

import {FatmanApp} from './app/app.component';
import {ReportService} from './app/reports/reports.service';

//enableProdMode();

// Since services are simple classes, we have to specify providers here
// instead of in some @Component section.
bootstrap(FatmanApp, [
         [ReportService, HTTP_PROVIDERS, ROUTER_PROVIDERS],
]);

//  vim: set ts=4 sw=2 tw=0 :
