import { Pipe, PipeTransform } from '@angular/core';

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

@Pipe({name: 'stringifydata'})
export class StringifyPipe implements PipeTransform {
  transform(value: any) {
    if (typeof value == 'number')
      return value.toFixed(0);

    if (typeof value == 'string')
      return value;

    // if the passed object is not an array (e.g. an arbitrary object)
    // figure out its keys => values and put them into an array
    if (!(value instanceof Array)) {
      let value_array: any[] = [];
      for (let key in value)
        if (value.hasOwnProperty(key))
          value_array.push(`${key}=` + JSON.stringify(value[key]));
      value = value_array;
    }

    if (value instanceof Array)
      return value.join(', ');

    return JSON.stringify(value);
  }
}

//  vim: set ts=2 sw=2 tw=0 :
