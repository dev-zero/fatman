
from flask import Flask
from playhouse.postgres_ext import PostgresqlExtDatabase

app = Flask(__name__)
app.config.from_object('fatman.default_settings')
app.config.from_envvar('FATMAN_SETTINGS', silent=True)

db = PostgresqlExtDatabase(app.config['DATABASE'], register_hstore=False)

import fatman.admin
import fatman.api
import fatman.models
