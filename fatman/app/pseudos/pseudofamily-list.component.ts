import { Component, OnInit } from '@angular/core';

import { Pseudo, PseudoFamily } from './pseudo';
import { PseudoService }        from './pseudo.service';
import { PseudoTable }          from './pseudo-table.component';

@Component({
  template: `
  <div class="row">
    <div class="col-md-12">
      <div class="form-group">
        <label for="familyselection">Select Pseudopotential Family:</label>
        <select [(ngModel)]="selectedFamily" class="form-control" id="familyselection">
          <option *ngFor="let family of pseudofamilies" [ngValue]="family">
            {{ family.name }}
          </option>
        </select>
      </div>
    </div>
  </div>
  <div class="row">
    <div class="col-md-12">
      <pseudotable [family]=selectedFamily></pseudotable>
    </div>
  </div>
  `,
  directives: [PseudoTable],
})

export class PseudoFamilyList implements OnInit {
  errorMessage: string;
  pseudofamilies: PseudoFamily[];
  selectedFamily: PseudoFamily;

  constructor(private _pseudoService: PseudoService) { }

  getPseudoFamilies() {
    this._pseudoService.getPseudoFamilies().subscribe(
      pseudofamilies => this.pseudofamilies = pseudofamilies,
      error          => this.errorMessage = <any>error);
  }

  ngOnInit() {
    this.getPseudoFamilies();
  }
}
