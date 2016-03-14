import {Component} from 'angular2/core';
import {Alert} from 'ng2-bootstrap/ng2-bootstrap';

@Component({
  selector: 'home',
  templateUrl: 'components/home/home.html',
  styleUrls: ['components/home/home.css'],
  providers: [],
  directives: [
    Alert
  ],
  pipes: []
})
export class Home {

  constructor() {}

}
