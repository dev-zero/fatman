
from flask import Flask
from playhouse.postgres_ext import PostgresqlExtDatabase
from flask.ext.uploads import UploadSet, configure_uploads

app = Flask(__name__)
app.config.from_object('fatman.default_settings')
app.config.from_envvar('FATMAN_SETTINGS', silent=True)

db = PostgresqlExtDatabase(app.config['DATABASE'], register_hstore=False)

resultfiles = UploadSet('results')
configure_uploads(app, (resultfiles,))

import fatman.admin
import fatman.api
import fatman.models
