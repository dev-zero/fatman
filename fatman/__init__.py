
import logging

from flask import Flask, redirect, url_for
from playhouse.postgres_ext import PostgresqlExtDatabase
from flask_uploads import UploadSet, configure_uploads
from flask_security import Security, PeeweeUserDatastore
from flask_caching import Cache

app = Flask(__name__)
app.config.from_object('fatman.default_settings')
app.config.from_pyfile('fatman.cfg', silent=True)
app.config.from_pyfile('../fatman.cfg', silent=True)
app.config.from_envvar('FATMAN_SETTINGS', silent=True)

if 'SECURITY_POST_LOGIN_VIEW' not in app.config:
    app.config['SECURITY_POST_LOGIN_VIEW'] = app.config['APPLICATION_ROOT']
if 'SECURITY_POST_LOGOUT_VIEW' not in app.config:
    app.config['SECURITY_POST_LOGOUT_VIEW'] = app.config['APPLICATION_ROOT']

if app.config.get('DATABASE_LOG_QUERIES', False):
    logger = logging.getLogger('peewee')
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())

db = PostgresqlExtDatabase(
    app.config['DATABASE'],
    autorollback=True, # automatically rollback when a query failed to avoid dead threads
    register_hstore=False,
    user=app.config.get('DATABASE_USER'),
    password=app.config.get('DATABASE_PASSWORD'),
    host=app.config.get('DATABASE_HOST'),
    port=app.config.get('DATABASE_PORT')
    )

# Explicitly connect
@app.before_request
def _db_connect():
    db.connect()

# Explicitly disconnect
@app.teardown_request
def _db_close(_):
    if not db.is_closed():
        db.close()

resultfiles = UploadSet('results')
configure_uploads(app, (resultfiles,))

# initialize Flask-Caching
cache = Cache(app, config=app.config)

# The imports are deliberately at this place.
# They import this file itself, but need the app and db objects to be ready.
# On the other hand we import them here to finish initialization of the app and db objects.
import fatman.models

from fatman.models import User, Role, UserRole
user_datastore = PeeweeUserDatastore(db, User, Role, UserRole)
security = Security(app, user_datastore)

import fatman.admin
import fatman.api
import fatman.views
