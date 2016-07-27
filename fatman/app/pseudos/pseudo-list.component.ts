import {Component, OnInit} from '@angular/core';

import {PseudoService} from './pseudo.service';

@Component({
  template: `
  <table class="tasks table table-bordered table-striped table-condensed">
    <thead>
      <tr>
          <th> Element </th>
          <th> Family </th>
          <th> Formats </th>
      </tr>
    </thead>
    <tbody>
        <tr *ngFor="let pseudo of pseudos">
            <td> {{pseudo.element}} </td>
            <td> {{pseudo.family}} </td>
            <td> {{pseudo.getAllFormats()}} </td>
        </tr>
    </tbody>
  </table>
  `,
})

export class PseudoList implements OnInit {
  title = 'Pseudos';

  errorMessage: string;
  pseudos: Object[];

  constructor(private _pseudoService: PseudoService) { }

  getPseudos() {
    this._pseudoService.getPseudos().subscribe(
      pseudos => this.pseudos = pseudos,
      error   => this.errorMessage = <any>error);
  }

  ngOnInit() {
    this.getPseudos();
  }
}
