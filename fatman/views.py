
from os import path
from werkzeug.exceptions import NotFound as WerkzeugNotFound
from flask import render_template, send_from_directory, redirect, url_for

from fatman import app

@app.route('/app/')
def webapp():
    """Serve the initial Angular page"""
    return render_template('index.html', application_root=app.config['APPLICATION_ROOT'])

@app.route('/app/<path:filename>')
def webappfile(filename):
    """If a corresponding file exists, return that,
    otherwise return the initial Angular page again which will
    read the route from the URL.

    This is mainly when running only the Flask webserver. Otherwise the
    front/proxy webserver should do this directly.
    """
    try:
        return send_from_directory(path.join(app.root_path, 'app'), filename)
    except WerkzeugNotFound:
        return render_template('index.html', application_root=app.config['APPLICATION_ROOT'])
