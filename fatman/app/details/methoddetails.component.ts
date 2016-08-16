import {Component, OnInit, Pipe, PipeTransform} from '@angular/core';
import {Router} from '@angular/router';

import {DetailsService} from './details.service';

/**
 * Iterable Pipe
 *
 * It accepts Objects and [Maps](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Map)
 * From here:
 *   https://gist.github.com/amcdnl/202596c5b85cc66d7002d10bde3ab514
 *
 * Example:
 *
 *  <div *ngFor="let keyValuePair of someObject | mapToIterable">
 *    key {{keyValuePair.key}} and value {{keyValuePair.value}}
 *  </div>
 *
 */
@Pipe({name: 'mapToIterable'})
export class IterablePipe implements PipeTransform {
  transform(iterable: any, args: any[]): any {
    let result = [];

    if (iterable.entries) {
      iterable.forEach((key, value) => {
        result.push({key, value});
      });
    } else {
      for (let key in iterable) {
        if (iterable.hasOwnProperty(key)) {
          result.push({key, value: iterable[key]});
        }
      }
    }

    return result;
  }
}

@Pipe({name: 'stringifysetting'})
export class SettingsPipe implements PipeTransform {
  transform(value: any) {
    if (typeof value == 'number')
      return value.toFixed(0);

    return JSON.stringify(value);
  }
}

@Pipe({name: 'concatkvpairlist'})
export class ConcatPipe implements PipeTransform {
  transform(iterable: any[]) {
    return iterable.map(entry => entry.key + ": " + JSON.stringify(entry.value)).join(", ");
  }
}

@Component({
  selector:'methoddetails', 
  inputs: ['method_id', 'small'],
  providers: [DetailsService],
  pipes: [SettingsPipe, IterablePipe, ConcatPipe],
  template: `
    <div *ngIf="method && !small" class="panel panel-default">
      <div class="panel-body">
        <div class="text-center pull-left" style="margin-right: 1em;">
          Method <p class="lead" style="font-size: 200%;"><strong>{{ method.id }}</strong></p>
        </div>
        <ul class="list-inline">
          <li><strong>Code:</strong> {{ method.code }}</li>
          <li><strong>Pseudopotential:</strong> {{ method.pseudopotential }}</li>
          <li><strong>Basis Set:</strong> {{ method.basis_set }}</li>
        </ul>
        <ul class="list-inline">
          <li *ngFor="let settingpair of method.settings | mapToIterable">
            {{settingpair.key}}: {{ settingpair.value | stringifysetting }}
          </li>
        </ul>
      </div>
    </div>
    <span *ngIf="method && small" title="{{ method.settings | mapToIterable | concatkvpairlist }}">
      <span title="{{ method.settings | mapToIterable | concatkvpairlist }}">
        <span class="text-muted">{{ method.id }}</span>
        <span class="small">{{ method.basis_set }} / {{ method.pseudopotential }} ({{ method.code }})</span>
      </span>
    </span>
  `,
})
export class MethoddetailsComponent implements OnInit {
  constructor(private _service: DetailsService) { }

  errorMessage: string;
  method: Object[];

  public method_id: number;

  getMethodDetails(method_id) {
    this._service.getMethodDetails(method_id).subscribe(
      method => this.method = method,
      error => this.errorMessage = <any>error);
  };

  ngOnInit(){
    this.getMethodDetails(this.method_id);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
