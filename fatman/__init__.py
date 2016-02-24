
from flask import Flask
from playhouse.postgres_ext import PostgresqlExtDatabase
from flask.ext.uploads import UploadSet, configure_uploads

app = Flask(__name__)
app.config.from_object('fatman.default_settings')
app.config.from_envvar('FATMAN_SETTINGS', silent=True)
app.debug = True

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
def _db_close(exc):
    if not db.is_closed():
        db.close()


resultfiles = UploadSet('results')
configure_uploads(app, (resultfiles,))

import fatman.admin
import fatman.api
import fatman.models
