
from flask import Flask
from playhouse.postgres_ext import PostgresqlExtDatabase

DATABASE = 'pseudotests'

app = Flask(__name__)
app.config.from_object(__name__)
db = PostgresqlExtDatabase(app.config['DATABASE'], register_hstore=False)
