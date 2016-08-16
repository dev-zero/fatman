export class Pseudo {

  public id: number = -1;
  public element: string = '';
  public family: string = '';
  public format: string = '';
  public converted_from: Pseudo;
  converted_pseudos: Pseudo[] = [];

  constructor(json : Object) {
    this.id = json['id'];
    this.element = json['element'];
    this.family = json['family'];
    this.format = json['format'];
    this.converted_from = json['converted_from'];
  }

  public addConverted(cp: Pseudo) {
    this.converted_pseudos.push(cp);
  }

  public getAllFormats() : string[] {
    return [this.format].concat(
      this.converted_pseudos.map(function(cp) { return cp.format; })
    );
  }
}
//  vim: set ts=2 sw=2 tw=0 :
