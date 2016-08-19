import {Component, OnInit} from '@angular/core';

import {Pseudo, PseudoFamily} from './pseudo';
import {PseudoService} from './pseudo.service';

@Component({
  template: `
  <div class="row">
    <div class="col-md-12">
      <div class="form-group">
        <label for="familyselection">Select Pseudopotential Family:</label>
        <select class="form-control" id="familyselection">
          <option></option>
          <option *ngFor="let family of pseudofamilies">
            {{ family.name }}
          </option>
        </select>
      </div>
    </div>
  </div>
  <div class="row">
    <div class="col-md-12">
    </div>
  </div>
  `,
})

export class PseudoFamilyList implements OnInit {
  errorMessage: string;
  pseudofamilies: PseudoFamily[];

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
