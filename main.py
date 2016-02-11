#!/usr/bin/env python3

from app import app, db

from models import *

def create_tables():
    db.create_tables([Structure, Method, Test, TestStructure, TaskStatus, Task, Results], safe=True)

if __name__ == '__main__':
    create_tables()
    app.run()
