
from peewee import *
from playhouse.postgres_ext import BinaryJSONField

from app import db

class BaseModel(Model):
    class Meta:
        database = db

class Structure(BaseModel):
    name = CharField(unique=True, null=False)
    class Meta:
        order_by = ('name',)

class Method(BaseModel):
    # proper database design would demand introduction of
    # separate entities for the following fields, but let's make it easy
    code = CharField()
    pseudopotential = CharField()
    basis_set = CharField()
    settings = BinaryJSONField(null=True)

    class Meta:
        order_by = ('code',)

class Test(BaseModel):
    name = CharField(unique=True, null=False)
    description = TextField(null=True)

    class Meta:
        order_by = ('name',)

class TestStructure(BaseModel):
    structure = ForeignKeyField(Structure, related_name='tests')
    test = ForeignKeyField(Test, related_name='structures')

class TaskStatus(BaseModel):
    name = CharField()

    class Meta:
        order_by = ('name',)

class Task(BaseModel):
    structure = ForeignKeyField(Structure, related_name='tasks')
    method = ForeignKeyField(Method, related_name='tasks')
    status = ForeignKeyField(TaskStatus)
    ctime = DateTimeField()
    mtime = DateTimeField()
    machine = CharField()

class Results(BaseModel):
    energy = DoubleField()
    task = ForeignKeyField(Task, related_name='results')

