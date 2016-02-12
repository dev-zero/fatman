
from flask import Flask
from playhouse.postgres_ext import PostgresqlExtDatabase

DATABASE = 'pseudotests'
SECRET_KEY = '06ff31da-fc90-4dc6-b177-26dd0018d2d7'

app = Flask(__name__)
app.config.from_object(__name__)
db = PostgresqlExtDatabase(app.config['DATABASE'], register_hstore=False)
