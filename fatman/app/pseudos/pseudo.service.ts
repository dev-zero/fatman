import { Http, Response, URLSearchParams } from '@angular/http';
import { Injectable } from '@angular/core';

import { Observable } from 'rxjs/Observable';

import { Pseudo, toPseudo, PseudoFamily, toPseudoFamily } from './pseudo';

@Injectable()
export class PseudoService {

  constructor(private _http: Http) {
  }

  private _pseudosUrl = '../api/v1/pseudos';
  private _pseudoFamiliesUrl = '../api/v1/pseudofamilies';

  getPseudoFamilies() : Observable<PseudoFamily[]> {
    return this._http.get(this._pseudoFamiliesUrl)
      .map(this._mapPseudoFamilies)
      .catch(this.handleError);
  }

  private _mapPseudoFamilies(response: Response) : PseudoFamily[] {
    return response.json().map(toPseudoFamily);
  }

  getPseudos(family: PseudoFamily) : Observable<Pseudo[]> {
    let params: URLSearchParams = new URLSearchParams();
    
    if (family)
      params.set('family', family.name)

    return this._http.get(this._pseudosUrl, {search: params})
      .map(this._mapPseudos)
      .catch(this.handleError);
  }

  private _mapPseudos(res: Response) : Pseudo[] {
    if (!res.json())
      return [];

    let plist = res.json().map(toPseudo);

    let pseudos : Pseudo[] = [];
    let unmatched : Pseudo[] = [];

    for (let pseudo of plist) {
      // add non-converted pseudos directly to the list of pseudos
      if (!pseudo.converted_from || pseudo.converted_from == null) {
        pseudos.push(pseudo);
        continue;
      }

      // we are using the fact that converted pseudos
      // are usually added shortly/directly after the original one
      let idx = -1;
      for (let i = pseudos.length-1; i >= 0; i--) {
        if (pseudos[i].id == pseudo.converted_from.id) {
          idx = i;
          break;
        }
      }

      // nevertheless, it may be that a converted pseudo was added before the original one
      if (idx >= 0) {
        pseudos[idx].converted_pseudos.push(pseudo);
      } else {
        unmatched.push(pseudo);
      }
    }

    // it is guaranteed that the converted_from id exists in the database if non-zero
    for (let pseudo of unmatched) {
      let spseudo = pseudos.find(function(p) { return p.id == pseudo.converted_from.id; })
      spseudo.converted_pseudos.push(pseudo);
    }

    return pseudos;
  }

  private handleError(error: any) {
    let errMsg = (error.message) ? error.message :
      error.status ? `${error.status} - ${error.statusText}` : 'Server error';
    console.error(errMsg);
    return Observable.throw(errMsg);
  }
}
//  vim: set ts=2 sw=2 tw=0 :
