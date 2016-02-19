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
def createconfig():
    """Create initial configuration file"""

    from os import urandom

    with open('fatman.cfg', 'w') as cfg:
        cfg.write("SECRET_KEY = {}\n".format(urandom(24)))

if __name__ == '__main__':
    manager.run()
