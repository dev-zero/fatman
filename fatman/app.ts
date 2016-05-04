
import {bootstrap} from '@angular/platform-browser-dynamic';
import {HTTP_PROVIDERS} from '@angular/http';
import {ROUTER_PROVIDERS} from '@angular/router';
//import { enableProdMode } from '@angular/core';
//enableProdMode();

import {FatmanApp} from './app/fatman-app';

import 'rxjs/Rx'

bootstrap(FatmanApp, [HTTP_PROVIDERS, ROUTER_PROVIDERS])
  .catch(err => console.error(err));
