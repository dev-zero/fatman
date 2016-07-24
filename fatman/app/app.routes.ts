
import {provideRouter, RouterConfig} from '@angular/router';

import {Home} from './home/home';
import {Tasks} from './tasks/tasks-list.component';
import {MachineStatus} from './machinestatus/machinestatus-list.component';
import {Reports} from './reports/reports-list.component';
import {ReportsComparison} from './reports/reports-comparison.component';
import {ReportsElementComparison} from './reports/reports-elementcomparison.component';
import {Periodictable} from './periodictable/periodictable-list.component';
import {Details} from './details/details-list.component';
import {MethodList} from './methods/methods-list.component';
import {PseudoList} from './pseudos/pseudo-list.component';


const routes: RouterConfig = [
  { path: '', redirectTo: '/home', pathMatch: 'full' },
  { path: 'home', component: Home },
  { path: 'tasks', component: Tasks },
  { path: 'reports/comparison/:id1/:id2', component: ReportsComparison },
  { path: 'reports/elementcomparison/:id1/:test1', component: ReportsElementComparison },
  { path: 'reports', component: Reports },
  { path: 'machinestatus', component: MachineStatus },
  { path: 'periodictable/:method1/:method2', component: Periodictable },
  { path: 'periodictable/', component: Home },
  { path: 'details/', component: Home },
  { path: 'details/:test1', component: Details },
  { path: 'methodlist', component: MethodList },
  { path: 'pseudos', component: PseudoList }
];

export const appRouterProviders = [
  provideRouter(routes)
];

//  vim: set ts=2 sw=2 tw=0 :
