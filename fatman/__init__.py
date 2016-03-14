
from flask import Flask, redirect, url_for
from playhouse.postgres_ext import PostgresqlExtDatabase
from flask.ext.uploads import UploadSet, configure_uploads
from flask.ext.security import Security, PeeweeUserDatastore

app = Flask(__name__)
app.config.from_object('fatman.default_settings')
app.config.from_pyfile('fatman.cfg', silent=True)
app.config.from_pyfile('../fatman.cfg', silent=True)
app.config.from_envvar('FATMAN_SETTINGS', silent=True)

if 'SECURITY_POST_LOGIN_VIEW' not in app.config:
    app.config['SECURITY_POST_LOGIN_VIEW'] = app.config['APPLICATION_ROOT']
if 'SECURITY_POST_LOGOUT_VIEW' not in app.config:
    app.config['SECURITY_POST_LOGOUT_VIEW'] = app.config['APPLICATION_ROOT']

db = PostgresqlExtDatabase(
    app.config['DATABASE'],
    autorollback=True, # automatically rollback when a query failed to avoid dead threads
    register_hstore=False)

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
