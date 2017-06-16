"""Defines CLI commands for FATMAN"""

import click
import flask_security.cli  # import commands from flask-security

from . import app


@app.cli.command()
def createconfig():
    """Create initial configuration file"""

    click.echo("Creating a minimal configuration for this FATMAN instance")

    from os import urandom

    with open('fatman.cfg', 'a') as cfg:
        cfg.write("SECRET_KEY = {}\n".format(urandom(24)))
