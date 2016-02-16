#!/usr/bin/env python

from flask.ext.script import Manager

from fatman import app, db

manager = Manager(app)

@manager.command
def initdb():
    """Initialize the database (structure and initial data)"""
    db.create_tables([Structure, Method, Test, TestStructure, TaskStatus, Task, Result], safe=True)

    for name in ['new', 'pending', 'running', 'done', 'error']:
        TaskStatus.create_or_get(name=name)

if __name__ == '__main__':
    manager.run()
