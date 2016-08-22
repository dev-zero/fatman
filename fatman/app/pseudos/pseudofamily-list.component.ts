import { Component, OnInit, OnDestroy } from '@angular/core';
import { Router, ActivatedRoute } from '@angular/router';
import { Location } from '@angular/common';

import { Pseudo, PseudoFamily } from './pseudo';
import { PseudoService }        from './pseudo.service';
import { PseudoTable }          from './pseudo-table.component';

@Component({
  template: `
  <div class="row">
    <div class="col-md-12">
      <div class="form-group">
        <label for="familyselection">Select Pseudopotential Family:</label>
        <select [(ngModel)]="selectedFamily" (ngModelChange)="updateUrl($event)"
            class="form-control" id="familyselection">
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

export class PseudoFamilyList implements OnInit, OnDestroy {
  errorMessage: string;
  pseudofamilies: PseudoFamily[];
  selectedFamily: PseudoFamily;
  selectedFamilyString: string;

  private _sub: any;

  constructor(
      private _route: ActivatedRoute,
      private _router: Router,
      private _location: Location,
      private _pseudoService: PseudoService) { }

  ngOnInit() {
    this._sub = this._route.params.subscribe(params => {
      this.selectedFamilyString = params['family'];
      this.getPseudoFamilies();
    });
  }

  ngOnDestroy() {
    this._sub.unsubscribe();
  }

  getPseudoFamilies() {
    this._pseudoService.getPseudoFamilies().subscribe(
      pseudofamilies => {
        this.pseudofamilies = pseudofamilies;

        if (this.selectedFamilyString)
          this.selectedFamily = this.pseudofamilies.find(family => family.name === this.selectedFamilyString);
      },
      error          => this.errorMessage = <any>error);
  }

  updateUrl(family: PseudoFamily) {
    this._location.replaceState(`/pseudos/${family.name}`)
  }
}
