import { Component, Input, SimpleChange,  Pipe, PipeTransform } from '@angular/core';

import { Pseudo } from './pseudo';

@Component({
  selector: 'pseudotable-element',
  template: `
    <ul *ngIf="pseudo" class="list-inline">
      <li>
        <a [href]="pseudo.links['self'] + '.' + pseudo.format | lowercase" target="_blank">{{ pseudo.format }}</a>
      </li>
      <li *ngFor="let pseudo of pseudo.converted_pseudos">
        <a [href]="pseudo.links['self'] + '.' + pseudo.format | lowercase" target="_blank">{{ pseudo.format }}</a>
      </li>
    </ul>
  `,
})

export class PseudoTableElement {
  @Input()
  pseudo: Pseudo;
}
