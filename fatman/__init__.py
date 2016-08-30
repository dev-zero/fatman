
import logging

from flask import Flask
from playhouse.postgres_ext import PostgresqlExtDatabase
from flask_uploads import UploadSet, configure_uploads
from flask_security import Security, PeeweeUserDatastore
from flask_caching import Cache
from celery import Celery

app = Flask(__name__)
app.config.from_object('fatman.default_settings')
app.config.from_pyfile('fatman.cfg', silent=True)
app.config.from_pyfile('../fatman.cfg', silent=True)
app.config.from_envvar('FATMAN_SETTINGS', silent=True)

if 'SECURITY_POST_LOGIN_VIEW' not in app.config:
    app.config['SECURITY_POST_LOGIN_VIEW'] = app.config['APPLICATION_ROOT']
if 'SECURITY_POST_LOGOUT_VIEW' not in app.config:
    app.config['SECURITY_POST_LOGOUT_VIEW'] = app.config['APPLICATION_ROOT']


def configure_db_logger():
    logger = logging.getLogger('peewee')
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())

if app.config.get('DATABASE_LOG_QUERIES', False):
    configure_db_logger()

db = PostgresqlExtDatabase(
    app.config['DATABASE'],
    autorollback=True,  # automatically rollback when a query failed to avoid dead threads
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
cache = Cache(app)


# initialize Celery
def setup_celery():
    # we need pickle support to serialize exceptions from
    # tasks, otherwise we would need a catch-all handler in them.
    app.config['CELERY_ACCEPT_CONTENT'] = ['pickle', 'json']

    capp = Celery(app.import_name,
                    backend=app.config['CELERY_BACKEND'],
                    broker=app.config['CELERY_BROKER_URL'])
    capp.conf.update(app.config)

    TaskBase = capp.Task

    # Inject the Flask context in Celery tasks
    class ContextTask(TaskBase):
        abstract = True

        def __call__(self, *args, **kwargs):
            with app.app_context():
                return TaskBase.__call__(self, *args, **kwargs)
    capp.Task = ContextTask

    return capp

capp = setup_celery()


from .models import User, Role, UserRole
user_datastore = PeeweeUserDatastore(db, User, Role, UserRole)
security = Security(app, user_datastore)


# The imports are deliberately at this place.
# They import this file itself, but need all other global objects to be ready.
# On the other hand we import them here to hook them up into Flask.
from . import models, views, admin, api
