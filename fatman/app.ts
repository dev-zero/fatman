
import {bootstrap} from 'angular2/platform/browser';
import {HTTP_PROVIDERS} from 'angular2/http';
import {ROUTER_PROVIDERS} from 'angular2/router';
//import { enableProdMode } from 'angular2/core';
//enableProdMode();

import {FatmanApp} from './app/fatman-app';

import 'rxjs/Rx'

bootstrap(FatmanApp, [HTTP_PROVIDERS, ROUTER_PROVIDERS])
  .catch(err => console.error(err));
