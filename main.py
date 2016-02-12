#!/usr/bin/env python3

from app import app, db

from models import *

def create_tables():
    db.create_tables([Structure, Method, Test, TestStructure, TaskStatus, Task, Results], safe=True)

def initialize_tables():
    for name in ['new', 'pending', 'running', 'done', 'error']:
        TaskStatus.create_or_get(name=name)

if __name__ == '__main__':
    create_tables()
    initialize_tables()
    app.run()
