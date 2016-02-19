#!/usr/bin/env python

from flask.ext.script import Manager

from fatman import app

manager = Manager(app)

@manager.command
def initdb():
    """Initialize the database (structure and initial data)"""

    from fatman import db
    from fatman.models import Structure, Method, Test, TestStructure, TaskStatus, Task, Result
    db.create_tables([Structure, Method, Test, TestStructure, TaskStatus, Task, Result], safe=True)

    for name in ['new', 'pending', 'running', 'done', 'error']:
        TaskStatus.create_or_get(name=name)


@manager.command
def cleardb():
    """Delete all the data in the Structure, Test, TestStructure, and Task databases. 
       The Result and Method tables remain intact."""

    from fatman import db
    from fatman.models import Structure, Test, TestStructure, Task
    for table in [Structure, TestStructure, Test, Task]:
        q = table.delete()
        q.execute()


@manager.command
def migration():
    """Do schema migration"""

    # the following is merely an example since it isn't idempotent
    # peewee does not have automated schema versioning/introspection
    from fatman import db
    from playhouse.migrate import PostgresqlMigrator, migrate
    from peewee import IntegerField
    # see http://docs.peewee-orm.com/en/latest/peewee/playhouse.html#schema-migrations
    migrator = PostgresqlMigrator(db)
    with db.transaction():
        migrate(
            # if a new field is to be not-NULL, a default must be provided,
            # the uniqueness constraint will get ignored here
            migrator.add_column('structure', 'ase_structure',
                                IntegerField(null=False, unique=True, default=0)),
            # ... and must be added manually via an index
            migrator.add_index('structure', ('ase_structure',), True),
        )


@manager.command
def createconfig():
    """Create initial configuration file"""

    from os import urandom

    with open('fatman.cfg', 'w') as cfg:
        cfg.write("SECRET_KEY = {}\n".format(urandom(24)))

@manager.shell
def make_shell_context():
    from fatman import db
    import fatman.models
    return dict(app=app, db=db, models=fatman.models)

if __name__ == '__main__':
    manager.run()
